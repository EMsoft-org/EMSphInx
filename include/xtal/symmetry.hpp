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

#ifndef _symmetry_h_
#define _symmetry_h_

#include <string>
#include <vector>
#include <functional>

namespace xtal {

	//@brief: a class to abstract point group operations
	//@note : I've opted for an enumeration type approach over polymorphism since there are finite crystallographic point groups
	//        This choice has benefits and drawbacks and either option can work
	//        I choose an enumeration style primarily to make the class easier to instantiate and to prevent the potential for confusion from pointer casting
	//        now that there is a general position class a point group could also be defined as a vector of general positions but that seems a bit cumbersome
	struct PointGroup {
		////////////////////////////////////////////////////////////////////////
		//               Constructors / Basic Attribute Queries               //
		////////////////////////////////////////////////////////////////////////

		//@brief   : construct a point group from a space group
		//@param sg: space group number [1,230]
		//@param as: should the alternate setting 2 group be used (e.g. P 1 1 2 instead of P 1 2 1 for group 3)
		//@note    : alternate setting is only supported for monoclinic groups [3,15]
		PointGroup(const uint_fast8_t sg = 1, const bool as = false);

		//@brief   : construct a point group from a string
		//@param pg: name of point group (e.g. m3m)
		//@note    : setting can be selected by using the full name e.g. "121" or "112" instead of just "2"
		//@note    : laue group and schonflies names can also be used
		PointGroup(std::string pg);

		//@brief   : a factory method to produce rotated point groups
		//@param pg: name of point group (e.g. 222r for 2.22)
		//@return  : special rotated point group
		//@note    : these exist for e.g. the rotational group of -4m2, only use them if you know what you're doing
		static PointGroup BuildOrtho45(std::string pg);

		//@brief : get a list of acceptable names to construct a point group from
		//@return: list of acceptable names
		static std::vector<std::string> Names();

		//@brief : get the point group number in IUCr order (1-32)
		//@return: point group number
		uint_fast8_t number() const;

		//@brief : determine if there are multiple axis choices for this point group (e.g. 112 vs 121 or -4m2 vs -42m)
		//@return: true/false for multiple/single axis choice(s)
		bool hasMultAxis() const;

		//@brief    : get the short Hermann-Mauguin name for the point group
		//@param lng: true/false to get a longer unambiguous name for multi setting groups (e.g. 1m1 or 11m instead of just m)
		//@return   : the name of the point group (e.g. "4/mmm")
		//@note     : rotoinversion axis are denoted as -n instead of \bar{n}
		std::string name(const bool lng = true) const;

		//@brief    : get the full Hermann-Mauguin name for the point group
		//@param lng: true/false to get a longer unambiguous name for multi setting groups (e.g. 1m1 or 11m instead of just m)
		//@return   : the name of the point group (e.g. "\frac{4}{m}\frac{2}{m}\frac{2}{m}")
		//@note     : formatted for latex typesetting
		std::string fullName(const bool lng = true) const;

		//@brief    : get the Schonflies name for the point group
		//@param alt: true/false to get the alternate symbol (e.g. S2 instead of Ci for -1)
		//@return   : the name of the point group (e.g. "D4h")
		//@note     : all characters after the first are subscripts
		std::string schonflies(const bool alt = false) const;

		//@brief    : get the Groth (1921) name for the point group
		//@param alt: true/false to get the alternate symbol (e.g. S2 instead of Ci for -1)
		//@return   : the name of the point group (e.g. "Ditetragonal-dipyramidal")
		//@note     : all characters after the first are subscripts
		std::string groth() const;

		//@brief    : get the Friedel (1926) name for the point group
		//@param alt: true/false to get the alternate symbol (e.g. S2 instead of Ci for -1)
		//@return   : the name of the point group (e.g. "Holohedry")
		//@note     : all characters after the first are subscripts
		std::string friedel() const;

		//@brief : get the TSL 'numbering' for this point group
		//@return: TSL number
		uint_fast8_t tslNum() const;

		//@brief : get the HKL 'numbering' for this point group (laue group number [1,11])
		//@return: HKL number
		uint_fast8_t hklNum() const;

		//@brief    : construct a point group from a TSL number
		//@param tsl: TSL Laue group number (e.g. 62)
		static PointGroup FromTSL(const uint_fast8_t tsl);

		//@brief    : construct a point group from a HKL number
		//@param hkl: HKL Laue group number (e.g. 62)
		static PointGroup FromHKL(const uint_fast8_t hkl);

		//@brief : get the order of the group
		//@return: order
		uint_fast8_t order() const;

		////////////////////////////////////////////////////////////////////////
		//                     Point Group Relationships                      //
		////////////////////////////////////////////////////////////////////////

		//@brief : get the name of the laue group this point group belongs to
		//@return: name of laue group e.g. "CubicLow"
		std::string laueName() const;

		//@brief : get the laue group the point group belongs to
		//@return: point group for the the laue group this point group belongs to
		PointGroup laueGroup() const;

		//@brief : get the purely rotational group the point group belongs to
		//@return: point group for the the purely rotational group this point group belongs to
		PointGroup rotationGroup() const;

		//@brief    : get the symmorphic space group of this point group
		//@param lat: lattice type (must be one of p, c, a, f, i, or r)
		//@return   : space group number (or 0 if this point group doesn't have a space group for the choosen lattice type)
		uint_fast8_t symmorphic(const char lat = 'p') const;

		//@brief    : get the transformation matrix from the default symmorphic setting to this space group
		//@return   : 3x3 transformation matrix A such that A^T * M * A is the symmetry operation m in the new reference frame (or null if in standard frame)
		//@note     : returns 45@z for 222r type groups, 90@x for 112 type groups, null otherwise
		template <typename Real> Real const * symmorphicTrns() const;

		////////////////////////////////////////////////////////////////////////
		//                        Symmetry Attributes                         //
		////////////////////////////////////////////////////////////////////////

		//@brief : check if this point group has inversion symmetry
		//@return: true/false if this point group does/doesn't have an inversino center
		//@note  : true/false also corresponds to centrosymmetric/non-centrosymmetric
		bool inversion() const;

		//@brief : check if this point group has enantiomorphism
		//@return: true/false if this crystal does / doesn't have enantiomorphism
		//@note  : true means no mirror planes or inversion symmetry
		bool enantiomorphism() const;

		//@brief : check if this point group has a mirror plane perpendicular to the z axis
		//@return: true/false if there is/isn't a mirror perpendicular to the z axis
		bool zMirror() const;

		//@brief : check if this point group has a zRot() planes with normals in the equator (e.g. Nmm where N == zRot() )
		//@return: 0 - there are no mirrors with normals in the equatorial plane
		//       : 1 - there are are mirrors with normals in the equatorial plane with normal  alignment (e.g. 31m and -4m2)
		//       : 2 - there are are mirrors with normals in the equatorial plane with rotated alignment (e.g. 31m and -42m)
		uint_fast8_t mmType() const;

		//@brief : get rotational symmetry about z axis
		//@return: degree of rotational symmetry about z
		uint_fast8_t zRot() const;

		//@brief : get the number of rotational symmetry operators
		//@return: number of rotational symmetry operators
		uint_fast8_t numRotOps() const;

		//@brief : get the closed set of rotational symmetry operators
		//@return: pointer to rotational symmetry operators (as w,x,y,z quaternions)
		template <typename Real> Real const * rotOps() const;

		//@brief : get the number of mirror planes
		//@return: number of mirror planes
		uint_fast8_t numMirror() const;

		//@brief : get mirror plane normals
		//@return: pointer to list of mirror plane normals (as x,y,z unit vectors)
		template <typename Real> Real const * mirrors() const;

		//@brief : get the number of rotational symmetry operators
		//@return: number of rotational symmetry operators
		uint_fast8_t numRotAxis() const;

		//@brief : get the rotational symmetry axis
		//@return: pointer to rotational symmetry axis (as n,x,y,z where xyz is a unit axis and n is the order (negative for rotoinversion))
		template <typename Real> Real const * rotAxis() const;

		////////////////////////////////////////////////////////////////////////
		//                        Symmetry Operations                         //
		////////////////////////////////////////////////////////////////////////

		//@brief   : check if a rodrigues vector is in the fundamental zone
		//@param ro: orientation to check (x, y, z, tan(w/2))
		//@return  : true/false if ro is/isn't in the fundamental zone
		template <typename Real> bool roInFz(Real const * const ro) const;

		//@brief   : compute the symmetric equivalent orientation in the fundamental zone
		//@param qu: orientation to compute symmetric equivalent of
		//@param fz: location to write symmetric equivalent
		template <typename Real> void fzQu(Real const * const qu, Real * const fz) const;

		//@brief    : compute the disorientation between 2 orientations
		//@param qu1: first  orientation to compute disorientation of (as w,x,y,z quaternion)
		//@param qu2: second orientation to compute disorientation of (as w,x,y,z quaternion)
		//@param dis: location to write disorientation (as w,x,y,z quaternion)
		//@note     : qu2 is assumed to belong to the same point group
		template <typename Real> void disoQu(Real const * const qu1, Real const * const qu2, Real * const dis) const;

		//@brief: compute the symmetric equivalent orientation of qu2 closest to qu1
		//@param qu1: orientation to search near
		//@param qu2: orientation to compute symmetric equivalent of
		//@param qu3: location to write symmetric equivalent of qu2 closest to qu1
		template <typename Real> void nearbyQu(Real const * const qu1, Real const * const qu2, Real * const equiv) const;

		//@brief   : reduce a unit direction to the fundamental sector
		//@param n : unit direction to reduce (magnitude is assumed to be 1)
		//@param fs: location to write symmetric equivalent of n in the fundamental sector (can be the same as n)
		//@return  : true if n was already in the fundamental sector
		template <typename Real> bool fsDir(Real const * const n, Real * const fs) const;

		//@param n  : crystal direction to color to color (magnitude is assumed to be 1)
		//@param rgb: location to write color [0,1]
		//@param h2r: hsl2rgb like coloring function to use with h, s, and l in [0,1] and output as [0,1] rgb: void(Real const * const hsl, Real * const rgb)
		template <typename Real> void ipfColor(Real const * const n, Real * const rgb, std::function<void(Real const*const, Real*const)> h2r) const;

		//@brief    : check if two point groups are the same
		//@param rhs: other point group to compare against
		//@return   : true if rhs is the same point group, false otherwise
		bool operator==(const PointGroup& rhs) const {return type == rhs.type;}

		//@brief    : check if two point groups are different
		//@param rhs: other point group to compare against
		//@return   : false if rhs is the same point group, true otherwise
		bool operator!=(const PointGroup& rhs) const {return type != rhs.type;}

		//@brief    : comparison (for sorting)
		//@param rhs: other point group to compare against
		//@return   : true if this points groups type is higher in the enumeration the ths
		bool operator<(const PointGroup& rhs) const {return type < rhs.type;}

		private:
			enum class PG {
				//triclinic
				   _1,// 1             "1"
				  _b1,// \bar{1}      "-1"

				//monoclinic (b unique preferred)
				 _121,// 121         "121" (b unique)
				 _112,// 112         "112" (c unique)
				 _1m1,// 1m1         "1m1" (b unique)
				 _11m,// 11m         "11m" (c unique)
				_12m1,// 12/m1     "12/m1" (b unique)
				_112m,// 112/m     "112/m" (c unique)

				//orthorhombic
				 _222,// 222         "222"
				_222r,// 222 rotated 45 degrees (2 fold axis at z, xy, and -xy) [this is 2.22 in the international tables] {see 10.1.3 for details}
				 _mm2,// mm2         "mm2"
				_mm2r,// mm2 rotated 45 degrees (as 222r, this is a maximal subgroup of -42m where -4m2 has mm2) [this would be m.m2 in the international tables]
				 _mmm,// mmm         "mmm"
				_mmmr,// mmm rotated 45 degrees (laue group of _222r) [this would be m.mm in the international tables]

				//tetragonal
				   _4,// 4             "4"
				  _b4,// \bar{4}      "-4"
				  _4m,// 4/m         "4/m"
				 _422,// 422         "422"
				 _4mm,// 4mm         "4mm"
				_b42m,// \bar{4}2m  "-42m" (2m vs m2)
				_b4m2,// \bar{4}m2  "-4m2" (2m vs m2)
				_4mmm,// 4/mmm     "4/mmm"
	 
				//trigonal
				   _3,// 3             "3"
				  _b3,// \bar{3}      "-3"
				 _321,// 321         "321" (x secondary)
				 _312,// 312         "312" (y secondary)
				 _3m1,// 3m1         "3m1" (x secondary)
				 _31m,// 31m         "31m" (y secondary)
				_b3m1,// \bar{3}m1  "-3m1" (x secondary)
				_b31m,// \bar{3}1m  "-31m" (y secondary)

				//hexagonal
				   _6,// 6             "6"
				  _b6,// \bar{6}      "-6"
				  _6m,// 6/m         "6/m"
				 _622,// 622         "622"
				 _6mm,// 6mm         "6mm"
				_b6m2,// \bar{6}m2  "-6m2" (2m vs m2)
				_b62m,// \bar{6}2m  "-62m" (2m vs m2)
				_6mmm,// 6/mmm     "6/mmm"

				//cubic
				  _23,// 23           "23"
				 _mb3,// m3           "m3"
				 _432,// 432         "432"
				_b43m,// \bar{4}3m  "-43m"
				_mb3m,// m\bar{3}m   "m3m"
			};

			PG type;//poing group type

			//@brief   : construct a point group from a point group enumeration
			//@param pg: point group enumeration
			PointGroup(const PG& pg) : type(pg) {}

			////////////////////////////////////////////////////////////////////////
			//            Static Functions for Symmetry Specific Code             //
			////////////////////////////////////////////////////////////////////////

			//@brief   : check if an orientation is in the fundamental zone
			//@param ro: orientation to check as [x, y, z, tan(w/2)]
			//@return  : true/false if the orientation is/isn't in the fundamental zone
			//@note    : these are specializations for cyclic, dihedral, and cubic symmetries
			template <typename Real> static bool FZ1   (Real const * const ro) {return ro[3] >= 0;}//negative tan(theta/2) -> outside of [0,pi] (negate all)}
			template <typename Real> static bool FZ121 (Real const * const ro);
			template <typename Real> static bool FZ112 (Real const * const ro);
			template <typename Real> static bool FZ222 (Real const * const ro);
			template <typename Real> static bool FZ222r(Real const * const ro);
			template <typename Real> static bool FZ3   (Real const * const ro);
			template <typename Real> static bool FZ321 (Real const * const ro);
			template <typename Real> static bool FZ312 (Real const * const ro);
			template <typename Real> static bool FZ4   (Real const * const ro);
			template <typename Real> static bool FZ422 (Real const * const ro);
			template <typename Real> static bool FZ6   (Real const * const ro);
			template <typename Real> static bool FZ622 (Real const * const ro);
			template <typename Real> static bool FZ23  (Real const * const ro);// 32
			template <typename Real> static bool FZ432 (Real const * const ro);// 432

			//@brief   : test x and y against 60 degree ro fz cutting plane
			//@param ro: orientation to check as [x, y, z, tan(w/2)]
			//@param yx: true to swap x and y
			//@return  : true if inside fundamental zone
			template <typename Real> static bool FZ3xy (Real const * const ro, const bool yx);

			//@brief    : compute the disorientation between 2 orientations
			//@param qu1: first  orientation to compute disorientation of
			//@param qu2: second orientation to compute disorientation of
			//@param dis: location to write disorientation
			//@param op1: rotational symmetry operators as quaternions (wxyz) for the first  orientation
			//@param no1: number of rotational symmetry operators for the first  orientation
			//@param op2: rotational symmetry operators as quaternions (wxyz) for the second orientation
			//@param no2: number of rotational symmetry operators for the second orientation
			template <typename Real> static void DisoQu(Real const * const qu1, Real const * const qu2, Real * const dis, Real const * const op1, const uint_fast8_t no1, Real const * const op2, const uint_fast8_t no2);

			//@brief    : compute the disorientation between 2 cubic (432) orientations
			//@param qu1: first  orientation to compute disorientation of
			//@param qu2: second orientation to compute disorientation of
			//@param dis: location to write disorientation
			//@note     : this shortcut is significantly faster than the general algorithm
			template <typename Real> static void Diso432(Real const * const qu1, Real const * const qu2, Real * const dis);
	};
}

////////////////////////////////////////////////////////////////////////
//                          Implementations                           //
////////////////////////////////////////////////////////////////////////

#include <algorithm>
#include <cctype>
#include <cmath>

#include "quaternion.hpp"
#include "rotations.hpp"
#include "sphere_sector.hpp"
#include "constants.hpp"

namespace xtal {

	////////////////////////////////////////////////////////////////////////
	//               Constructors / Basic Attribute Queries               //
	////////////////////////////////////////////////////////////////////////

	//@brief   : construct a point group from a space group
	//@param sg: space group number [1,230]
	//@param as: should the alternate setting 2 group be used (e.g. P 1 1 2 instead of P 1 2 1 for group 3)
	//@note    : alternate setting is only supported for monoclinic groups [3,15]
	PointGroup::PointGroup(const uint_fast8_t sg, const bool as) {
		//look up table to convert from a space group to a point group
		static const PG SG2PG[230]  = {
			//triclinic
			PG::   _1, PG::  _b1,                                            //  1-  2

			//monoclinic
			PG:: _121, PG:: _121, PG:: _121,                                 //  3-  5
			PG:: _1m1, PG:: _1m1, PG:: _1m1, PG:: _1m1,                      //  6-  9
			PG::_12m1, PG::_12m1, PG::_12m1, PG::_12m1, PG::_12m1, PG::_12m1,// 10- 15

			//orthorhombic
			PG:: _222, PG:: _222, PG:: _222, PG:: _222, PG:: _222, PG:: _222,// 16- 24
			PG:: _222, PG:: _222, PG:: _222,                                 
			PG:: _mm2, PG:: _mm2, PG:: _mm2, PG:: _mm2, PG:: _mm2, PG:: _mm2,// 25- 46
			PG:: _mm2, PG:: _mm2, PG:: _mm2, PG:: _mm2, PG:: _mm2, PG:: _mm2,
			PG:: _mm2, PG:: _mm2, PG:: _mm2, PG:: _mm2, PG:: _mm2, PG:: _mm2,
			PG:: _mm2, PG:: _mm2, PG:: _mm2, PG:: _mm2,                      
			PG:: _mmm, PG:: _mmm, PG:: _mmm, PG:: _mmm, PG:: _mmm, PG:: _mmm,// 47- 74
			PG:: _mmm, PG:: _mmm, PG:: _mmm, PG:: _mmm, PG:: _mmm, PG:: _mmm,
			PG:: _mmm, PG:: _mmm, PG:: _mmm, PG:: _mmm, PG:: _mmm, PG:: _mmm,
			PG:: _mmm, PG:: _mmm, PG:: _mmm, PG:: _mmm, PG:: _mmm, PG:: _mmm,
			PG:: _mmm, PG:: _mmm, PG:: _mmm, PG:: _mmm,

			//tetragonal
			PG::   _4, PG::   _4, PG::   _4, PG::   _4, PG::   _4, PG::   _4,// 75- 80
			PG::  _b4, PG::  _b4,                                            // 81- 82
			PG::  _4m, PG::  _4m, PG::  _4m, PG::  _4m, PG::  _4m, PG::  _4m,// 83- 88
			PG:: _422, PG:: _422, PG:: _422, PG:: _422, PG:: _422, PG:: _422,// 89- 98
			PG:: _422, PG:: _422, PG:: _422, PG:: _422,
			PG:: _4mm, PG:: _4mm, PG:: _4mm, PG:: _4mm, PG:: _4mm, PG:: _4mm,// 99-110
			PG:: _4mm, PG:: _4mm, PG:: _4mm, PG:: _4mm, PG:: _4mm, PG:: _4mm,
			PG::_b42m, PG::_b42m, PG::_b42m, PG::_b42m,                      //111-114
			PG::_b4m2, PG::_b4m2, PG::_b4m2, PG::_b4m2, PG::_b4m2, PG::_b4m2,//115-120
			PG::_b42m, PG::_b42m,                                            //121-122
			PG::_4mmm, PG::_4mmm, PG::_4mmm, PG::_4mmm, PG::_4mmm, PG::_4mmm,//123-142
			PG::_4mmm, PG::_4mmm, PG::_4mmm, PG::_4mmm, PG::_4mmm, PG::_4mmm,
			PG::_4mmm, PG::_4mmm, PG::_4mmm, PG::_4mmm, PG::_4mmm, PG::_4mmm,
			PG::_4mmm, PG::_4mmm,

			//trigonal
			PG::   _3, PG::   _3, PG::   _3, PG::   _3,                      //143-146
			PG::  _b3, PG::  _b3,                                            //147-148
			PG:: _312, PG:: _321, PG:: _312, PG:: _321, PG:: _312, PG:: _321,//149-155
			PG:: _321,
			PG:: _3m1, PG:: _31m, PG:: _3m1, PG:: _31m, PG:: _3m1, PG:: _3m1,//156-161
			PG::_b31m, PG::_b31m, PG::_b3m1, PG::_b3m1, PG::_b3m1, PG::_b3m1,//162-167

			//hexagonal
			PG::   _6, PG::   _6, PG::   _6, PG::   _6, PG::   _6, PG::   _6,//168-173
			PG::  _b6,                                                       //174
			PG::  _6m, PG::  _6m,                                            //175-176
			PG:: _622, PG:: _622, PG:: _622, PG:: _622, PG:: _622, PG:: _622,//177-182
			PG:: _6mm, PG:: _6mm, PG:: _6mm, PG:: _6mm,                      //183-186
			PG::_b6m2, PG::_b6m2, PG::_b62m, PG::_b62m,                      //187-190
			PG::_6mmm, PG::_6mmm, PG::_6mmm, PG::_6mmm,                      //191-194

			//cubic
			PG::  _23, PG::  _23, PG::  _23, PG::  _23, PG::  _23,           //195-199
			PG:: _mb3, PG:: _mb3, PG:: _mb3, PG:: _mb3, PG:: _mb3, PG:: _mb3,//200-206
			PG:: _mb3,
			PG:: _432, PG:: _432, PG:: _432, PG:: _432, PG:: _432, PG:: _432,//207-214
			PG:: _432, PG:: _432,
			PG::_b43m, PG::_b43m, PG::_b43m, PG::_b43m, PG::_b43m, PG::_b43m,//215-220
			PG::_mb3m, PG::_mb3m, PG::_mb3m, PG::_mb3m, PG::_mb3m, PG::_mb3m,//221-230
			PG::_mb3m, PG::_mb3m, PG::_mb3m, PG::_mb3m
		};

		if(sg < 1 || sg > 230) throw std::runtime_error("point group number must be [1,230]");
		type = SG2PG[sg-1];
	}

	//@brief   : construct a point group from a string
	//@param pg: name of point group (e.g. m3m)
	//@note    : setting can be selected by using the full name e.g. "121" or "112" instead of just "2"
	//@note    : laue group and schonflies names can also be used
	PointGroup::PointGroup(std::string pg) {
		//first make a lowercase copy of the point group name with all white space removed
		pg.erase(remove_if(pg.begin(), pg.end(), ::isspace), pg.end());//remove white space
		std::transform(pg.begin(), pg.end(), pg.begin(), ::tolower);//convert to lower case

		//now make a list of all possible point groups
		//note that _222r, _mm2r, and _mmmr are excluded since they should be for internal use only
		// 222r can be constructed as PointGroup("-4m2").rotationGroup();
		// mmmr can be constructed as PointGroup("-4m2").rotationGroup().laueGroup();
		// mm2r can't currently be constructed publicly
		const PG groups[40] = {
			PG::   _1, PG::  _b1, PG:: _121, PG:: _112, PG:: _1m1, PG:: _11m,
			PG::_12m1, PG::_112m, PG:: _222, PG:: _mm2, PG:: _mmm, PG::   _4,
			PG::  _b4, PG::  _4m, PG:: _422, PG:: _4mm, PG::_b42m, PG::_b4m2,
			PG::_4mmm, PG::   _3, PG::  _b3, PG:: _312, PG:: _321, PG:: _31m,
			PG:: _3m1, PG::_b31m, PG::_b3m1, PG::   _6, PG::  _b6, PG::  _6m,
			PG:: _622, PG:: _6mm, PG::_b6m2, PG::_b62m, PG::_6mmm, PG::  _23,
			PG:: _mb3, PG:: _432, PG::_b43m, PG::_mb3m
		};

		//Schonflies start with a capital letter
		std::string pgUp = pg;
		if(!pg.empty()) pgUp.front() = ::toupper(pg.front());

		//loop over point groups checking names
		for(size_t i = 0; i < 40; i++) {
			PointGroup g(groups[i]);
			bool match = false;
			if     (0 == pg.compare(g.name(true ) )) match = true;//unambigous HM name
			else if(0 == pg.compare(g.name(false) )) match = true;//ambigous HM name (e.g. "2" instead of "121")
			else {
				if(0 == pgUp.compare(g.schonflies())) match = true;//schonflies name
			}
			if(match) {
				type = g.type;
				return;
			}
		}

		//if a point group didn't work, try laue group names
		if     (0 == pg.compare("triclinic"     )) {type = PG::  _b1; return;}
		else if(0 == pg.compare("monoclinic"    )) {type = PG::_12m1; return;}
		else if(0 == pg.compare("orthorhombic"  )) {type = PG:: _mmm; return;}
		else if(0 == pg.compare("trigonallow"   )) {type = PG::  _b3; return;}
		else if(0 == pg.compare("trigonalhigh"  )) {type = PG::_b3m1; return;}
		else if(0 == pg.compare("trigonal"      )) {type = PG::_b3m1; return;}//assume high symmetry if not specified
		else if(0 == pg.compare("tetragonallow" )) {type = PG::  _4m; return;}
		else if(0 == pg.compare("tetragonalhigh")) {type = PG::_4mmm; return;}
		else if(0 == pg.compare("tetragonal"    )) {type = PG::_4mmm; return;}//assume high symmetry if not specified
		else if(0 == pg.compare("hexagonallow"  )) {type = PG::  _6m; return;}
		else if(0 == pg.compare("hexagonalhigh" )) {type = PG::_6mmm; return;}
		else if(0 == pg.compare("hexagonal"     )) {type = PG::_6mmm; return;}//assume high symmetry if not specified
		else if(0 == pg.compare("cubiclow"      )) {type = PG:: _mb3; return;}
		else if(0 == pg.compare("cubichigh"     )) {type = PG::_mb3m; return;}
		else if(0 == pg.compare("cubic"         )) {type = PG::_mb3m; return;}//assume high symmetry if not specified

		//if we didn't find a match the string isn't valid
		throw std::runtime_error("couldn't parse a point group type from '" + pg + "'");
	}

	//@brief   : a factory method to produce rotated point groups
	//@param pg: name of point group (e.g. 222r for 2.22)
	//@return  : special rotated point group
	//@note    : these exist for e.g. the rotational group of -4m2, only use them if you know what you're doing
	PointGroup PointGroup::BuildOrtho45(std::string pg) {
		if     ("222r" == pg) return PointGroup(PG::_222r);
		else if("mmmr" == pg) return PointGroup(PG::_mmmr);
		else if("mm2r" == pg) return PointGroup(PG::_mm2r);
		else throw std::runtime_error("invalid name for 45 degree rotated orthorhombic group '" + pg + "'");
	}


	//@brief : get a list of acceptable names to construct a point group from
	//@return: list of acceptable names
	std::vector<std::string> PointGroup::Names() {
		//make a list of all possible point groups
		const PG groups[40] = {
			PG::   _1, PG::  _b1, PG:: _121, PG:: _112, PG:: _1m1, PG:: _11m,
			PG::_12m1, PG::_112m, PG:: _222, PG:: _mm2, PG:: _mmm, PG::   _4,
			PG::  _b4, PG::  _4m, PG:: _422, PG:: _4mm, PG::_b42m, PG::_b4m2,
			PG::_4mmm, PG::   _3, PG::  _b3, PG:: _312, PG:: _321, PG:: _31m,
			PG:: _3m1, PG::_b31m, PG::_b3m1, PG::   _6, PG::  _b6, PG::  _6m,
			PG:: _622, PG:: _6mm, PG::_b6m2, PG::_b62m, PG::_6mmm, PG::  _23,
			PG:: _mb3, PG:: _432, PG::_b43m, PG::_mb3m
		};

		//accumulate HM names
		std::vector<std::string> names;
		for(size_t i = 0; i < 40; i++) {
			PointGroup g(groups[i]);
			names.push_back(g.name(true));
			const std::string nm = g.name(false);//get ambigous name as well (e.g. "32" instead of "321")
			if(0 != names.back().compare(nm)) names.push_back(nm);//add ambigous name if it is different
		}

		//accmulate schonflies names
		for(size_t i = 0; i < 40; i++) {
			PointGroup g(groups[i]);
			const std::string nm = g.schonflies();
			if(0 != names.back().compare(nm)) names.push_back(nm);
		}

		//accumulate laue group names
		names.push_back("Triclinic"     );// -1
		names.push_back("Monoclinic"    );// 2/m
		names.push_back("Orthorhombic"  );// mmm
		names.push_back("TrigonalLow"   );// -3
		names.push_back("TrigonalHigh"  );// -3m
		names.push_back("Trigonal"      );// -3m
		names.push_back("TetragonalLow" );// 4/m
		names.push_back("TetragonalHigh");// 4/mmm
		names.push_back("Tetragonal"    );// 4/mmm
		names.push_back("HexagonalLow"  );// 6/m
		names.push_back("HexagonalHigh" );// 6/mmm
		names.push_back("Hexagonal"     );// 6/mmm
		names.push_back("CubicLow"      );// m3
		names.push_back("CubicHigh"     );// m3m
		names.push_back("Cubic"         );// m3m
		return names;
	}

	//@brief : get the point group number in IUCr order (1-32)
	//@return: point group number
	uint_fast8_t PointGroup::number() const {
		switch(type) {
			case PG::   _1: return  1;
			case PG::  _b1: return  2;
			case PG:: _121://intentional fall through
			case PG:: _112: return  3;
			case PG:: _1m1://intentional fall through
			case PG:: _11m: return  4;
			case PG::_12m1://intentional fall through
			case PG::_112m: return  5;
			case PG:: _222://intentional fall through
			case PG::_222r: return  6;
			case PG:: _mm2://intentional fall through
			case PG::_mm2r: return  7;
			case PG:: _mmm://intentional fall through
			case PG::_mmmr: return  8;
			case PG::   _4: return  9;
			case PG::  _b4: return 10;
			case PG::  _4m: return 11;
			case PG:: _422: return 12;
			case PG:: _4mm: return 13;
			case PG::_b42m://intentional fall through
			case PG::_b4m2: return 14;
			case PG::_4mmm: return 15;
			case PG::   _3: return 16;
			case PG::  _b3: return 17;
			case PG:: _321://intentional fall through
			case PG:: _312: return 18;
			case PG:: _3m1://intentional fall through
			case PG:: _31m: return 19;
			case PG::_b3m1://intentional fall through
			case PG::_b31m: return 20;
			case PG::   _6: return 21;
			case PG::  _b6: return 22;
			case PG::  _6m: return 23;
			case PG:: _622: return 24;
			case PG:: _6mm: return 25;
			case PG::_b62m://intentional fall through
			case PG::_b6m2: return 26;
			case PG::_6mmm: return 27;
			case PG::  _23: return 28;
			case PG:: _mb3: return 29;
			case PG:: _432: return 30;
			case PG::_b43m: return 31;
			case PG::_mb3m: return 32;
		}
		return -1;
	}

	//@brief : determine if there are multiple axis choices for this point group (e.g. 112 vs 121 or -4m2 vs -42m)
	//@return: true/false for multiple/single axis choice(s)
	bool PointGroup::hasMultAxis() const {
		switch(type) {
			case PG:: _121:
			case PG:: _112:
			case PG:: _1m1:
			case PG:: _11m:
			case PG::_12m1:
			case PG::_112m:
			case PG::_b42m:
			case PG::_b4m2:
			case PG:: _312:
			case PG:: _321:
			case PG:: _31m:
			case PG:: _3m1:
			case PG::_b31m:
			case PG::_b3m1:
			case PG::_b6m2:
			case PG::_b62m: return true;

			case PG::   _1:
			case PG::  _b1:
			case PG:: _222://these could be considered as multiple conventions
			case PG::_222r://but 222r is extremely unusual (and shouldn't be used)
			case PG:: _mm2://these could be considered as multiple conventions
			case PG::_mm2r://but mm2r is extremely unusual (and shouldn't be used)
			case PG:: _mmm://these could be considered as multiple conventions
			case PG::_mmmr://but mmmr is extremely unusual (and shouldn't be used)
			case PG::   _4:
			case PG::  _b4:
			case PG::  _4m:
			case PG:: _422:
			case PG:: _4mm:
			case PG::   _3:
			case PG::  _b3:
			case PG::   _6:
			case PG::  _b6:
			case PG::  _6m:
			case PG:: _622:
			case PG:: _6mm:
			case PG::_4mmm:
			case PG::_6mmm:
			case PG::  _23:
			case PG:: _mb3:
			case PG:: _432:
			case PG::_b43m:
			case PG::_mb3m: return false;
		}
		return false;
	}

	//@brief    : get the short Hermann-Mauguin name for the point group
	//@param lng: true/false to get a longer unambiguous name for multi setting groups (e.g. 1m1 or 11m instead of just m)
	//@return   : the name of the point group (e.g. "4/mmm")
	//@note     : rotoinversion axis are denoted as -n instead of \bar{n}
	std::string PointGroup::name(const bool lng) const {
		switch(type) {
			case PG::   _1: return     "1";
			case PG::  _b1: return    "-1";
			case PG:: _121: return lng ?   "121" :   "2";
			case PG:: _112: return   "112";//always use full name for non default ambigous
			case PG:: _1m1: return lng ?   "1m1" :   "m";
			case PG:: _11m: return   "11m";//always use full name for non default ambigous
			case PG::_12m1: return lng ? "12/m1" : "2/m";
			case PG::_112m: return "112/m";//always use full name for non default ambigous
			case PG:: _222: return   "222";
			case PG::_222r: return  "222r";
			case PG:: _mm2: return   "mm2";
			case PG::_mm2r: return  "mm2r";
			case PG:: _mmm: return   "mmm";
			case PG::_mmmr: return  "mmmr";
			case PG::   _4: return     "4";
			case PG::  _b4: return    "-4";
			case PG::  _4m: return   "4/m";
			case PG:: _422: return   "422";
			case PG:: _4mm: return   "4mm";
			case PG::_b42m: return  "-42m";
			case PG::_b4m2: return  "-4m2";
			case PG::_4mmm: return "4/mmm";
			case PG::   _3: return     "3";
			case PG::  _b3: return    "-3";
			case PG:: _321: return lng ?  "321" :  "32";
			case PG:: _312: return   "312";//always use full name for non default ambigous
			case PG:: _3m1: return lng ?  "3m1" :  "3m";
			case PG:: _31m: return   "31m";//always use full name for non default ambigous
			case PG::_b3m1: return lng ? "-3m1" : "-3m";
			case PG::_b31m: return  "-31m";//always use full name for non default ambigous
			case PG::   _6: return     "6";
			case PG::  _b6: return    "-6";
			case PG::  _6m: return   "6/m";
			case PG:: _622: return   "622";
			case PG:: _6mm: return   "6mm";
			case PG::_b62m: return  "-62m";
			case PG::_b6m2: return  "-6m2";
			case PG::_6mmm: return "6/mmm";
			case PG::  _23: return    "23";
			case PG:: _mb3: return    "m3";
			case PG:: _432: return   "432";
			case PG::_b43m: return  "-43m";
			case PG::_mb3m: return   "m3m";
		}
		return "";
	}

	//@brief    : get the full Hermann-Mauguin name for the point group
	//@param lng: true/false to get a longer unambiguous name for multi setting groups (e.g. 1m1 or 11m instead of just m)
	//@return   : the name of the point group (e.g. "\frac{4}{m}\frac{2}{m}\frac{2}{m}")
	//@note     : formatted for latex typesetting
	std::string PointGroup::fullName(const bool lng) const {
		switch(type) {
			case PG::   _1: return "1"                                                   ;
			case PG::  _b1: return "\\bar{1}"                                            ;
			case PG:: _121: return lng ? "121" : "2"                                     ;
			case PG:: _112: return "112"                                                 ;//always use full name for non default ambigous
			case PG:: _1m1: return lng ? "1m1" : "m"                                     ;
			case PG:: _11m: return "11m"                                                 ;//always use full name for non default ambigous
			case PG::_12m1: return lng ? "1\\frac{2}{m}1" : "\\frac{2}{m}"               ;
			case PG::_112m: return "11\\frac{2}{m}"                                      ;//always use full name for non default ambigous
			case PG:: _222: return "222"                                                 ;
			case PG::_222r: return "222r"                                                ;
			case PG:: _mm2: return "mm2"                                                 ;
			case PG::_mm2r: return "mm2r"                                                ;
			case PG:: _mmm: return "\\frac{2}{m}\\frac{2}{m}\\frac{2}{m}"                ;
			case PG::_mmmr: return "\\frac{2}{m}\\frac{2}{m}\\frac{2}{m}r"               ;
			case PG::   _4: return "4"                                                   ;
			case PG::  _b4: return "\\bar{4}"                                            ;
			case PG::  _4m: return "\\frac{4}{m}"                                        ;
			case PG:: _422: return "422"                                                 ;
			case PG:: _4mm: return "4mm"                                                 ;
			case PG::_b42m: return "\\bar{4}2m"                                          ;
			case PG::_b4m2: return "\\bar{4}m2"                                          ;
			case PG::_4mmm: return "\\frac{4}{m}\\frac{2}{m}\\frac{2}{m}"                ;
			case PG::   _3: return "3"                                                   ;
			case PG::  _b3: return "\\bar{3}"                                            ;
			case PG:: _321: return lng ? "321" : "32"                                    ;
			case PG:: _312: return "312"                                                 ;//always use full name for non default ambigous
			case PG:: _3m1: return lng ? "3m1" : "3m"                                    ;
			case PG:: _31m: return "31m"                                                 ;//always use full name for non default ambigous
			case PG::_b3m1: return lng ? "\\bar{3}\\frac{2}{m}1" : "\\bar{3}\\frac{2}{m}";
			case PG::_b31m: return "\\bar{3}1\\frac{2}{m}"                               ;//always use full name for non default ambigous
			case PG::   _6: return "6"                                                   ;
			case PG::  _b6: return "\\bar{6}"                                            ;
			case PG::  _6m: return "\\frac{6}{m}"                                        ;
			case PG:: _622: return "622"                                                 ;
			case PG:: _6mm: return "6mm"                                                 ;
			case PG::_b62m: return "\\bar{6}2m"                                          ;
			case PG::_b6m2: return "\\bar{6}m2"                                          ;
			case PG::_6mmm: return "\\frac{6}{m}\\frac{2}{m}\\frac{2}{m}"                ;
			case PG::  _23: return "23"                                                  ;
			case PG:: _mb3: return "\\frac{2}{m}\\bar{3}"                                ;
			case PG:: _432: return "432"                                                 ;
			case PG::_b43m: return "\\bar{4}3m"                                          ;
			case PG::_mb3m: return "\\frac{4}{m}\\bar{3}\\frac{2}{m}"                    ;
		}
		return "";
	}

	//@brief    : get the Schonflies name for the point group
	//@param alt: true/false to get the alternate symbol (e.g. S2 instead of Ci for -1)
	//@return   : the name of the point group (e.g. "D4h")
	//@note     : all characters after the first are subscripts
	std::string PointGroup::schonflies(const bool alt) const {
		switch(type) {
			case PG::   _1: return "C1" ;
			case PG::  _b1: return alt ? "S2"  : "Ci" ;
			case PG:: _121:
			case PG:: _112: return "C2" ;
			case PG:: _1m1:
			case PG:: _11m: return alt ? "C1h" : "Cs" ;
			case PG::_12m1:
			case PG::_112m: return "C2h";
			case PG:: _222:
			case PG::_222r: return alt ? "V"   : "D2" ;
			case PG:: _mm2:
			case PG::_mm2r: return "C2v";
			case PG:: _mmm:
			case PG::_mmmr: return alt ? "Vh"  : "D2h";
			case PG::   _4: return "C4" ;
			case PG::  _b4: return "S4" ;
			case PG::  _4m: return "C4h";
			case PG:: _422: return "D4" ;
			case PG:: _4mm: return "C4v";
			case PG::_b42m:
			case PG::_b4m2: return alt ? "Vd"  : "D2d";
			case PG::_4mmm: return "D4h";
			case PG::   _3: return "C3" ;
			case PG::  _b3: return alt ? "S6"  : "C3i";
			case PG:: _321:
			case PG:: _312: return "D3" ;
			case PG:: _3m1:
			case PG:: _31m: return "C3v";
			case PG::_b3m1:
			case PG::_b31m: return "D3d";
			case PG::   _6: return "C6" ;
			case PG::  _b6: return "C3h";
			case PG::  _6m: return "C6h";
			case PG:: _622: return "D6" ;
			case PG:: _6mm: return "C6v";
			case PG::_b62m:
			case PG::_b6m2: return "D3h";
			case PG::_6mmm: return "D6h";
			case PG::  _23: return "T"  ;
			case PG:: _mb3: return "Th" ;
			case PG:: _432: return "O"  ;
			case PG::_b43m: return "Td" ;
			case PG::_mb3m: return "Oh" ;
		}
		return "";
	}

	//@brief    : get the Groth (1921) name for the point group
	//@param alt: true/false to get the alternate symbol (e.g. S2 instead of Ci for -1)
	//@return   : the name of the point group (e.g. "Ditetragonal-dipyramidal")
	//@note     : all characters after the first are subscripts
	std::string PointGroup::groth() const {
		switch(type) {
			case PG::   _1: return "pedial (asymmetric)"                           ;
			case PG::  _b1: return "pinacoidal"                                    ;
			case PG:: _121:
			case PG:: _112: return "sphenoidal"                                    ;
			case PG:: _1m1:
			case PG:: _11m: return "domatic"                                       ;
			case PG::_12m1:
			case PG::_112m: return "prismatic"                                     ;
			case PG:: _222:
			case PG::_222r: return "disphenoidal"                                  ;
			case PG:: _mm2:
			case PG::_mm2r: return "pyramidal"                                     ;
			case PG:: _mmm:
			case PG::_mmmr: return "dipyramidal"                                   ;
			case PG::   _4: return "pryamidal"                                     ;
			case PG::  _b4: return "disphenoidal"                                  ;
			case PG::  _4m: return "dipyramidal"                                   ;
			case PG:: _422: return "trapezoihedral"                                ;
			case PG:: _4mm: return "ditetragonal-pyramidal"                        ;
			case PG::_b42m:
			case PG::_b4m2: return "scalenohedral"                                 ;
			case PG::_4mmm: return "ditetragonal-dipyramidal"                      ;
			case PG::   _3: return "pyramidal"                                     ;
			case PG::  _b3: return "rhombohedral"                                  ;
			case PG:: _321:
			case PG:: _312: return "trapezohedral"                                 ;
			case PG:: _3m1:
			case PG:: _31m: return "ditrigonal-pyramidal"                          ;
			case PG::_b3m1:
			case PG::_b31m: return "ditrigonal-scalenohedral"                      ;
			case PG::   _6: return "pyramidal"                                     ;
			case PG::  _b6: return "trigonal-dipyramidal"                          ;
			case PG::  _6m: return "dipyramidal"                                   ;
			case PG:: _622: return "trapezohedral"                                 ;
			case PG:: _6mm: return "dihexagonal-pyramidal"                         ;
			case PG::_b62m:
			case PG::_b6m2: return "ditrigonal-dipyramidal"                        ;
			case PG::_6mmm: return "dihexagonal-dipyramidal"                       ;
			case PG::  _23: return "tetrahedral-pentagondodecahedral (tetartoidal)";
			case PG:: _mb3: return "disdodecahderal (diploidal)"                   ;
			case PG:: _432: return "pentagon-icositetrahedral (gyroidal)"          ;
			case PG::_b43m: return "hexakistetrahedral (hextetrahedral)"           ;
			case PG::_mb3m: return "hexakisoctahedral (hexoctahedral)"             ;
		}
		return "";
	}

	//@brief    : get the Friedel (1926) name for the point group
	//@param alt: true/false to get the alternate symbol (e.g. S2 instead of Ci for -1)
	//@return   : the name of the point group (e.g. "Holohedry")
	//@note     : all characters after the first are subscripts
	std::string PointGroup::friedel() const {
		switch(type) {
			case PG::   _1: return "hemihedry"                         ;
			case PG::  _b1: return "holohedry"                         ;
			case PG:: _121:
			case PG:: _112: return "holoaxial hemihedry"               ;
			case PG:: _1m1:
			case PG:: _11m: return "antihemihedry"                     ;
			case PG::_12m1:
			case PG::_112m: return "holohedry"                         ;
			case PG:: _222:
			case PG::_222r: return "holoaxial hemihedry"               ;
			case PG:: _mm2:
			case PG::_mm2r: return "antihemihedry"                     ;
			case PG:: _mmm:
			case PG::_mmmr: return "holohedry"                         ;
			case PG::   _4: return "tetartohedry with 4-axis"          ;
			case PG::  _b4: return "sphenohedral tartohedry"           ;
			case PG::  _4m: return "parahemihedry"                     ;
			case PG:: _422: return "holoaxial hemihedry"               ;
			case PG:: _4mm: return "antihemihedry with 4-axis"         ;
			case PG::_b42m:
			case PG::_b4m2: return "sphenohderal antihemihedry"        ;
			case PG::_4mmm: return "holohedry"                         ;
			case PG::   _3: return "ogdohedry"                         ;//rhombohedral: tetartohedry
			case PG::  _b3: return "paratetartohedry"                  ;//rhombohedral: parahemihedry
			case PG:: _321:
			case PG:: _312: return "holoaxial tetartohedry with 3-axis";//rhombohedral: holoaxial hemihedry
			case PG:: _3m1:
			case PG:: _31m: return "hemimorphic antitetartohedry"      ;//rhombohedral: antihemihedry
			case PG::_b3m1:
			case PG::_b31m: return "parahemihedry with 3-axis"         ;//rhombohedral: holohedry
			case PG::   _6: return "tetartohedry with 6-axis"          ;
			case PG::  _b6: return "trigonohedral antitetartohedry"    ;
			case PG::  _6m: return "parahemihedry with 6-axis"         ;
			case PG:: _622: return "holoaxial hemihedry"               ;
			case PG:: _6mm: return "antihemihedry with 6-axis"         ;
			case PG::_b62m:
			case PG::_b6m2: return "trigonohedral antihemihedry"       ;
			case PG::_6mmm: return "holohedry"                         ;
			case PG::  _23: return "tetartohedry"                      ;
			case PG:: _mb3: return "parahemihedry"                     ;
			case PG:: _432: return "holoaxial hemihedry"               ;
			case PG::_b43m: return "antihemihedry"                     ;
			case PG::_mb3m: return "holohedry"                         ;
		}
		return "";
	}

	//@brief : get the TSL 'numbering' for this point group
	//@return: TSL number
	uint_fast8_t PointGroup::tslNum() const {
		switch(laueGroup().type) {
			case PG::  _b1: return  1;
			case PG::_12m1:
			case PG::_112m: return  2;
			case PG:: _mmm:
			case PG::_mmmr: return 22;
			case PG::  _4m: return  4;
			case PG::_4mmm: return 42;
			case PG::  _b3: return  3;
			case PG::_b3m1:
			case PG::_b31m: return 32;
			case PG::  _6m: return  6;
			case PG::_6mmm: return 62;
			case PG:: _mb3: return 23;
			case PG::_mb3m: return 43;
			default       : throw std::logic_error("laueGroup() '" + laueGroup().name() + "' not handled for tslNum()");
		}
	}

	//@brief : get the HKL 'numbering' for this point group (laue group number [1,11])
	//@return: HKL number
	uint_fast8_t PointGroup::hklNum() const {
		switch(laueGroup().type) {
			case PG::  _b1: return  1;
			case PG::_12m1:
			case PG::_112m: return  2;
			case PG:: _mmm:
			case PG::_mmmr: return  3;
			case PG::  _4m: return  4;
			case PG::_4mmm: return  5;
			case PG::  _b3: return  6;
			case PG::_b3m1:
			case PG::_b31m: return  7;
			case PG::  _6m: return  8;
			case PG::_6mmm: return  9;
			case PG:: _mb3: return 10;
			case PG::_mb3m: return 11;
			default       : throw std::logic_error("laueGroup() '" + laueGroup().name() + "' not handled for hklNum()");
		}
	}

	//@brief    : construct a point group from a TSL number
	//@param tsl: TSL Laue group number (e.g. 62)
	PointGroup PointGroup::FromTSL(const uint_fast8_t tsl) {
		switch(tsl) {
			case  1: return PointGroup(   "-1");
			case  2: return PointGroup(  "2/m");
			case 22: return PointGroup(  "mmm");
			case  4: return PointGroup(  "4/m");
			case 42: return PointGroup("4/mmm");
			case  3: return PointGroup(   "-3");
			case 32: return PointGroup(  "-3m");
			case  6: return PointGroup(  "6/m");
			case 62: return PointGroup("6/mmm");
			case 23: return PointGroup(   "m3");
			case 43: return PointGroup(  "m3m");
			default: throw std::domain_error("couldn't get point group from TSL number");
		}
	}

	//@brief    : construct a point group from a HKL number
	//@param hkl: HKL Laue group number (e.g. 11)
	PointGroup PointGroup::FromHKL(const uint_fast8_t hkl) {
		switch(hkl) {
			case  1: return PointGroup(   "-1");
			case  2: return PointGroup(  "2/m");
			case  3: return PointGroup(  "mmm");
			case  4: return PointGroup(  "4/m");
			case  5: return PointGroup("4/mmm");
			case  6: return PointGroup(   "-3");
			case  7: return PointGroup(  "-3m");
			case  8: return PointGroup(  "6/m");
			case  9: return PointGroup("6/mmm");
			case 10: return PointGroup(   "m3");
			case 11: return PointGroup(  "m3m");
			default: throw std::domain_error("couldn't get point group from HKL number");
		}
	}

	//@brief : get the order of the group
	//@return: order
	uint_fast8_t PointGroup::order() const {
		switch(type) {
			case PG::   _1: return  1;
			case PG::  _b1:
			case PG:: _121:
			case PG:: _112:
			case PG:: _1m1:
			case PG:: _11m: return  2;
			case PG::_12m1:
			case PG::_112m:
			case PG:: _222:
			case PG::_222r:
			case PG:: _mm2:
			case PG::_mm2r:
			case PG::   _4:
			case PG::  _b4: return  4;
			case PG:: _mmm:
			case PG::_mmmr:
			case PG::  _4m:
			case PG:: _422:
			case PG:: _4mm:
			case PG::_b42m:
			case PG::_b4m2: return  8;
			case PG::_4mmm: return 16;
			case PG::   _3: return  3;
			case PG::  _b3:
			case PG:: _321:
			case PG:: _312:
			case PG:: _3m1:
			case PG:: _31m:
			case PG::   _6:
			case PG::  _b6: return  6;
			case PG::_b3m1:
			case PG::_b31m:
			case PG::  _6m:
			case PG:: _622:
			case PG:: _6mm:
			case PG::_b62m:
			case PG::_b6m2:
			case PG::  _23: return 12;
			case PG::_6mmm:
			case PG:: _mb3:
			case PG:: _432:
			case PG::_b43m: return 24;
			case PG::_mb3m: return 48;
		}
		return 0;
	}

	////////////////////////////////////////////////////////////////////////
	//                     Point Group Relationships                      //
	////////////////////////////////////////////////////////////////////////

	//@brief : get the name of the laue group this point group belongs to
	//@return: name of laue group e.g. "CubicLow"
	std::string PointGroup::laueName() const {
		switch(laueGroup().type) {
			case PG::  _b1: return "Triclinic"    ;
			case PG::_12m1:
			case PG::_112m: return "Monoclinic"   ;
			case PG:: _mmm:
			case PG::_mmmr: return "Orthorhombic" ;
			case PG::  _4m: return "TetragonalLow";
			case PG::_4mmm: return "Tetragonal"   ;
			case PG::  _b3: return "TrigonalLow"  ;
			case PG::_b3m1:
			case PG::_b31m: return "Trigonal"     ;
			case PG::  _6m: return "HexagonalLow" ;
			case PG::_6mmm: return "Hexagonal"    ;
			case PG:: _mb3: return "CubicLow"     ;
			case PG::_mb3m: return "Cubic"        ;
			default       : throw std::logic_error("laueGroup() '" + laueGroup().name() + "' not handled for laueName()");
		}
	}

	//@brief : get the laue group the point group belongs to
	//@return: point group for the the laue group this point group belongs to
	PointGroup PointGroup::laueGroup() const {
		switch(type) {
			case PG::   _1:
			case PG::  _b1: return PointGroup(PG::  _b1);
			case PG:: _121:
			case PG:: _1m1:
			case PG::_12m1: return PointGroup(PG::_12m1);
			case PG:: _112:
			case PG:: _11m:
			case PG::_112m: return PointGroup(PG::_112m);
			case PG:: _222:
			case PG:: _mm2:
			case PG:: _mmm: return PointGroup(PG:: _mmm);
			case PG::_222r: 
			case PG::_mm2r:
			case PG::_mmmr: return PointGroup(PG::_mmmr);
			case PG::   _4:
			case PG::  _b4:
			case PG::  _4m: return PointGroup(PG::  _4m);
			case PG:: _422:
			case PG:: _4mm:
			case PG::_b42m:
			case PG::_b4m2:
			case PG::_4mmm: return PointGroup(PG::_4mmm);
			case PG::   _3:
			case PG::  _b3: return PointGroup(PG::  _b3);
			case PG:: _321:
			case PG:: _3m1:
			case PG::_b3m1: return PointGroup(PG::_b3m1);
			case PG:: _312:
			case PG:: _31m:
			case PG::_b31m: return PointGroup(PG::_b31m);
			case PG::   _6:
			case PG::  _b6:
			case PG::  _6m: return PointGroup(PG::  _6m);
			case PG:: _622:
			case PG:: _6mm:
			case PG::_b6m2:
			case PG::_b62m:
			case PG::_6mmm: return PointGroup(PG::_6mmm);
			case PG::  _23:
			case PG:: _mb3: return PointGroup(PG:: _mb3);
			case PG:: _432:
			case PG::_b43m:
			case PG::_mb3m: return PointGroup(PG::_mb3m);
		}
		return PointGroup(PG::_b1);
	}

	//@brief : get the purely rotational group the point group belongs to
	//@return: point group for the the purely rotational group this point group belongs to
	PointGroup PointGroup::rotationGroup() const {
		switch(type) {
			case PG::   _1:
			case PG::  _b1:
			case PG:: _1m1:
			case PG:: _11m: return PointGroup(PG::   _1);
			case PG:: _121:
			case PG::_12m1: return PointGroup(PG:: _121);
			case PG:: _112:
			case PG::_112m:
			case PG:: _mm2:
			case PG::_mm2r:
			case PG::  _b4: return PointGroup(PG:: _112);
			case PG:: _222:
			case PG:: _mmm:
			case PG::_b42m: return PointGroup(PG:: _222);
			case PG::_222r:
			case PG::_mmmr:
			case PG::_b4m2: return PointGroup(PG::_222r);//222 rotated 45 degrees
			case PG::   _4:
			case PG::  _4m:
			case PG:: _4mm: return PointGroup(PG::   _4);
			case PG:: _422:
			case PG::_4mmm: return PointGroup(PG:: _422);
			case PG::   _3:
			case PG::  _b3:
			case PG:: _3m1:
			case PG:: _31m:
			case PG::  _b6: return PointGroup(PG::   _3);
			case PG:: _321:
			case PG::_b3m1:
			case PG::_b62m: return PointGroup(PG:: _321);
			case PG:: _312:
			case PG::_b31m:
			case PG::_b6m2: return PointGroup(PG:: _312);
			case PG::   _6:
			case PG::  _6m:
			case PG:: _6mm: return PointGroup(PG::   _6);
			case PG:: _622:
			case PG::_6mmm: return PointGroup(PG:: _622);
			case PG::  _23:
			case PG:: _mb3:
			case PG::_b43m: return PointGroup(PG::  _23);
			case PG:: _432:
			case PG::_mb3m: return PointGroup(PG:: _432);
		}
		return PointGroup(PG::_1);
	}

	//@brief    : get the symmorphic space group of this point group
	//@param lat: lattice type (must be one of p, c, a, f, i, or r)
	//@return   : space group number (or 0 if this point group doesn't have a space group for the choosen lattice type)
	uint_fast8_t PointGroup::symmorphic(const char lat) const {
		const char x = std::toupper(lat);// e.g. 'p' -> 'P'
		switch(type) {
			case PG::   _1: return 'P' == x ? 1 : 0;
			case PG::  _b1: return 'P' == x ? 2 : 0;

			//monoclinic b (alternate cell choices currently commented out)
			case PG:: _121: return 'P' == x ?  3 : ( ('C' == x /*|| 'A' == x || 'I' == x*/) ?  5 : 0);
			case PG:: _1m1: return 'P' == x ?  6 : ( ('C' == x /*|| 'A' == x || 'I' == x*/) ?  8 : 0);
			case PG::_12m1: return 'P' == x ? 10 : ( ('C' == x /*|| 'A' == x || 'I' == x*/) ? 12 : 0);

			//monoclinic c (alternate cell choices currently commented out)
			case PG:: _112: return 'P' == x ?  3 : ( ('A' == x /*|| 'B' == x || 'I' == x*/) ?  5 : 0);
			case PG:: _11m: return 'P' == x ?  6 : ( ('A' == x /*|| 'B' == x || 'I' == x*/) ?  8 : 0);
			case PG::_112m: return 'P' == x ? 10 : ( ('A' == x /*|| 'B' == x || 'I' == x*/) ? 12 : 0);

			//orthorhombic (currently only abc setting)
			case PG::_222: return 'P' == x ? 16 : ( 'C' == x ? 21 : ( 'F' == x ? 22 : ( 'I' == x ? 23 :                   0   ) ) );
			case PG::_mm2: return 'P' == x ? 25 : ( 'C' == x ? 35 : ( 'F' == x ? 42 : ( 'I' == x ? 44 : ( 'A' == x ? 38 : 0 ) ) ) );
			case PG::_mmm: return 'P' == x ? 47 : ( 'C' == x ? 65 : ( 'F' == x ? 69 : ( 'I' == x ? 71 :                   0   ) ) );

			//rotated orthorhombic
			case PG::_222r: return 'P' == x ? 16 : ( 'C' == x ? 21 : ( 'F' == x ? 22 : ( 'I' == x ? 23 :                   0   ) ) );
			case PG::_mm2r: return 'P' == x ? 25 : ( 'C' == x ? 35 : ( 'F' == x ? 42 : ( 'I' == x ? 44 : ( 'A' == x ? 38 : 0 ) ) ) );
			case PG::_mmmr: return 'P' == x ? 47 : ( 'C' == x ? 65 : ( 'F' == x ? 69 : ( 'I' == x ? 71 :                   0   ) ) );

			//tetragonal
			case PG::   _4: return 'P' == x ?  75 : ('I' == x ?  79 : 0);
			case PG::  _b4: return 'P' == x ?  81 : ('I' == x ?  82 : 0);
			case PG::  _4m: return 'P' == x ?  83 : ('I' == x ?  87 : 0);
			case PG:: _422: return 'P' == x ?  89 : ('I' == x ?  97 : 0);
			case PG:: _4mm: return 'P' == x ?  99 : ('I' == x ? 107 : 0);
			case PG::_b42m: return 'P' == x ? 111 : ('I' == x ? 121 : 0);
			case PG::_b4m2: return 'P' == x ? 115 : ('I' == x ? 119 : 0);
			case PG::_4mmm: return 'P' == x ? 123 : ('I' == x ? 139 : 0);

			//trigonal
			case PG::   _3: return 'P' == x ? 143 : ('R' == x ? 146 : 0);
			case PG::  _b3: return 'P' == x ? 147 : ('R' == x ? 148 : 0);
			case PG:: _321: return 'P' == x ? 150 : ('R' == x ? 155 : 0);
			case PG:: _312: return 'P' == x ? 149 : ('R' == x ? 155 : 0);
			case PG:: _3m1: return 'P' == x ? 156 : ('R' == x ? 160 : 0);
			case PG:: _31m: return 'P' == x ? 157 : ('R' == x ? 160 : 0);
			case PG::_b3m1: return 'P' == x ? 164 : ('R' == x ? 166 : 0);
			case PG::_b31m: return 'P' == x ? 162 : ('R' == x ? 166 : 0);

			//hexagonal
			case PG::   _6: return 'P' == x ? 168 : 0;
			case PG::  _b6: return 'P' == x ? 174 : 0;
			case PG::  _6m: return 'P' == x ? 175 : 0;
			case PG:: _622: return 'P' == x ? 177 : 0;
			case PG:: _6mm: return 'P' == x ? 183 : 0;
			case PG::_b6m2: return 'P' == x ? 187 : 0;
			case PG::_b62m: return 'P' == x ? 189 : 0;
			case PG::_6mmm: return 'P' == x ? 191 : 0;

			//cubic
			case PG::  _23: return 'P' == x ? 195 : ('F' == x ? 196 : ('I' == x ? 197 : 0 ) ) ;
			case PG:: _mb3: return 'P' == x ? 200 : ('F' == x ? 202 : ('I' == x ? 204 : 0 ) ) ;
			case PG:: _432: return 'P' == x ? 207 : ('F' == x ? 209 : ('I' == x ? 211 : 0 ) ) ;
			case PG::_b43m: return 'P' == x ? 215 : ('F' == x ? 216 : ('I' == x ? 217 : 0 ) ) ;
			case PG::_mb3m: return 'P' == x ? 221 : ('F' == x ? 225 : ('I' == x ? 229 : 0 ) ) ;

			default: return 0;
		}
	}

	//@brief    : get the transformation matrix from the default symmorphic setting to this space group
	//@return   : 3x3 transformation matrix A such that A^T * M * A is the symmetry operation m in the new reference frame (or null if in standard frame)
	//@note     : returns 45@z for 222r type groups, 90@x for 112 type groups, null otherwise
	template <typename Real>
	Real const * PointGroup::symmorphicTrns() const {
		static const Real monoC[9] = {1, 0, 0,  0, 0, -1,  0, 1, 1};//90@x
		static const Real orthR[9] = {//45@z
			Constants<Real>::r1_2, -Constants<Real>::r1_2, 0,
			Constants<Real>::r1_2,  Constants<Real>::r1_2, 0,
			                    0,                      0, 1
        };

		switch(type) {
			case PG::   _1:
			case PG::  _b1:
			case PG:: _121:
			case PG:: _1m1:
			case PG::_12m1:
			case PG:: _222:
			case PG:: _mm2:
			case PG:: _mmm:
			case PG::   _4:
			case PG::  _b4:
			case PG::  _4m:
			case PG:: _422:
			case PG:: _4mm:
			case PG::_b42m:
			case PG::_b4m2:
			case PG::_4mmm:
			case PG::   _3:
			case PG::  _b3:
			case PG:: _321:
			case PG:: _312:
			case PG:: _3m1:
			case PG:: _31m:
			case PG::_b3m1:
			case PG::_b31m:
			case PG::   _6:
			case PG::  _b6:
			case PG::  _6m:
			case PG:: _622:
			case PG:: _6mm:
			case PG::_b6m2:
			case PG::_b62m:
			case PG::_6mmm:
			case PG::  _23:
			case PG:: _mb3:
			case PG:: _432:
			case PG::_b43m:
			case PG::_mb3m: return NULL;

			case PG:: _112:
			case PG:: _11m:
			case PG::_112m: return monoC;

			case PG::_222r:
			case PG::_mm2r:
			case PG::_mmmr: return orthR;
		}
		return NULL;
	}

	////////////////////////////////////////////////////////////////////////
	//                        Symmetry Attributes                         //
	////////////////////////////////////////////////////////////////////////

	//@brief : check if this point group has inversion symmetry
	//@return: true/false if this point group does/doesn't have an inversino center
	bool PointGroup::inversion() const {
		switch(type) {
			case PG::   _1:
			case PG:: _121:
			case PG:: _112:
			case PG:: _1m1:
			case PG:: _11m:
			case PG:: _222:
			case PG::_222r:
			case PG:: _mm2:
			case PG::_mm2r:
			case PG::   _4:
			case PG::  _b4:
			case PG:: _422:
			case PG:: _4mm:
			case PG::_b42m:
			case PG::_b4m2:
			case PG::   _3:
			case PG:: _321:
			case PG:: _312:
			case PG:: _3m1:
			case PG:: _31m:
			case PG::   _6:
			case PG::  _b6:
			case PG:: _622:
			case PG:: _6mm:
			case PG::_b6m2:
			case PG::_b62m:
			case PG::  _23:
			case PG:: _432:
			case PG::_b43m: return false;

			case PG::  _b1:
			case PG::_12m1:
			case PG::_112m:
			case PG:: _mmm:
			case PG::_mmmr:
			case PG::  _4m:
			case PG::_4mmm:
			case PG::  _b3:
			case PG::_b3m1:
			case PG::_b31m:
			case PG::  _6m:
			case PG::_6mmm:
			case PG:: _mb3:
			case PG::_mb3m: return true ;
		}
		return false;
	}

	//@brief : check if this point group has enantiomorphism
	//@return: true/false if this crystal does / doesn't have enantiomorphism
		//@note  : true means no mirror planes or inversion symmetry
	bool PointGroup::enantiomorphism() const {
		switch(type) {
			case PG::   _1:
			case PG:: _121:
			case PG:: _112:
			case PG:: _222:
			case PG::_222r:
			case PG::   _4:
			case PG:: _422:
			case PG::   _3:
			case PG:: _321:
			case PG:: _312:
			case PG::   _6:
			case PG:: _622:
			case PG::  _23:
			case PG:: _432: return true;

			case PG::  _b1:
			case PG:: _1m1:
			case PG:: _11m:
			case PG::_12m1:
			case PG::_112m:
			case PG:: _mm2:
			case PG::_mm2r:
			case PG:: _mmm:
			case PG::_mmmr:
			case PG::  _b4:
			case PG::  _4m:
			case PG:: _4mm:
			case PG::_b42m:
			case PG::_b4m2:
			case PG::_4mmm:
			case PG::  _b3:
			case PG:: _3m1:
			case PG:: _31m:
			case PG::_b3m1:
			case PG::_b31m:
			case PG::  _b6:
			case PG::  _6m:
			case PG:: _6mm:
			case PG::_b62m:
			case PG::_b6m2:
			case PG::_6mmm:
			case PG:: _mb3:
			case PG::_b43m:
			case PG::_mb3m: return false;
		}
		return false;
	}

	//@brief : check if this point group has a mirror plane perpendicular to the z axis
	//@return: true/false if there is/isn't a mirror perpendicular to the z axis
	bool PointGroup::zMirror() const {
		switch(type) {
			case PG::   _1:
			case PG::  _b1:
			case PG:: _121:
			case PG:: _112:
			case PG:: _1m1:
			case PG::_12m1:
			case PG:: _222:
			case PG::_222r:
			case PG:: _mm2:
			case PG::_mm2r:
			case PG::   _4:
			case PG::  _b4:
			case PG:: _422:
			case PG:: _4mm:
			case PG::_b42m:
			case PG::_b4m2:
			case PG::   _3:
			case PG::  _b3:
			case PG:: _321:
			case PG:: _312:
			case PG:: _3m1:
			case PG:: _31m:
			case PG::_b3m1:
			case PG::_b31m:
			case PG::   _6:
			case PG:: _622:
			case PG:: _6mm:
			case PG::  _23:
			case PG:: _432:
			case PG::_b43m: return false;

			case PG:: _11m:
			case PG::_112m:
			case PG:: _mmm:
			case PG::_mmmr:
			case PG::  _4m:
			case PG::_4mmm:
			case PG::  _b6:
			case PG::  _6m:
			case PG::_b6m2:
			case PG::_b62m:
			case PG::_6mmm:
			case PG:: _mb3:
			case PG::_mb3m: return true ;
		}
		return false;
	}

	//@brief : check if this point group has a zRot() planes with normals in the equator (e.g. Nmm where N == zRot() )
	//@return: 0 - there are no mirrors with normals in the equatorial plane
	//       : 1 - there are are mirrors with normals in the equatorial plane with normal  alignment (e.g. 31m and -4m2)
	//       : 2 - there are are mirrors with normals in the equatorial plane with rotated alignment (e.g. 31m and -42m)
	uint_fast8_t PointGroup::mmType() const {
		switch(type) {
			case PG::   _1: return 0;
			case PG::  _b1: return 0;
			case PG:: _121: return 0;
			case PG:: _112: return 0;
			case PG:: _1m1: return 1;
			case PG:: _11m: return 0;
			case PG::_12m1: return 1;
			case PG::_112m: return 0;
			case PG:: _222: return 0;
			case PG::_222r: return 0;
			case PG:: _mm2: return 1;
			case PG::_mm2r: return 2;
			case PG:: _mmm: return 1;
			case PG::_mmmr: return 2;
			case PG::   _4: return 0;
			case PG::  _b4: return 0;
			case PG::  _4m: return 0;
			case PG:: _422: return 0;
			case PG:: _4mm: return 1;
			case PG::_b42m: return 2;
			case PG::_b4m2: return 1;
			case PG::_4mmm: return 1;
			case PG::   _3: return 0;
			case PG::  _b3: return 0;
			case PG:: _321: return 0;
			case PG:: _312: return 0;
			case PG:: _3m1: return 2;
			case PG:: _31m: return 1;
			case PG::_b3m1: return 2;
			case PG::_b31m: return 1;
			case PG::   _6: return 0;
			case PG::  _b6: return 0;
			case PG::  _6m: return 0;
			case PG:: _622: return 0;
			case PG:: _6mm: return 1;
			case PG::_b6m2: return 2;
			case PG::_b62m: return 1;
			case PG::_6mmm: return 1;
			case PG::  _23: return 0;
			case PG:: _mb3: return 1;
			case PG:: _432: return 0;
			case PG::_b43m: return 2;
			case PG::_mb3m: return 1;
		}
		return 0;
	}

	//@brief : get rotational symmetry about z axis
	//@return: degree of rotational symmetry about z
	uint_fast8_t PointGroup::zRot() const {
		switch(rotationGroup().type) {
			case PG::   _1:
			case PG:: _121: return 1;
			case PG:: _112: return 2;
			case PG:: _222:
			case PG::_222r:
			case PG::  _23: return 2;
			case PG::   _3:
			case PG:: _321:
			case PG:: _312: return 3;
			case PG::   _4:
			case PG:: _422:
			case PG:: _432: return 4;
			case PG::   _6:
			case PG:: _622: return 6;
			default       : throw std::logic_error("rotationGroup() '" + rotationGroup().name() + "' not handled for numRotOps()");
		}
	}

	//@brief : get the closed set of rotational symmetry operators
	//@return: number of rotational symmetry operators
	uint_fast8_t PointGroup::numRotOps() const {
		switch(rotationGroup().type) {
			case PG::   _1: return  1;
			case PG:: _121:
			case PG:: _112: return  2;
			case PG:: _222:
			case PG::_222r:
			case PG::   _4: return  4;
			case PG:: _422: return  8;
			case PG::   _3: return  3;
			case PG:: _321:
			case PG:: _312:
			case PG::   _6: return  6;
			case PG:: _622:
			case PG::  _23: return 12;
			case PG:: _432: return 24;
			default       : throw std::logic_error("rotationGroup() '" + rotationGroup().name() + "' not handled for numRotOps()");
		}
	}

	//@brief : get the rotational symmetry operators
	//@return: pointer to rotational symmetry operators (as w,x,y,z quaternions)
	template <typename Real>
	Real const * PointGroup::rotOps() const {
		//cubic type symmetry symmetry operators
		//1, 112, 121, 211, 222, 222r, 4, 422, 23, and 432 area all contigous subsets (with only 3 extra values beyond 432)
		static Real const cub[27*4] = {
			          Real(0.00) ,           Real(1.00) ,           Real(0.00) ,           Real(0.00) ,//180 @ x        [          211                      ] (duplicate for contigous 211     )
			          Real(1.00) ,           Real(0.00) ,           Real(0.00) ,           Real(0.00) ,//identity       [          211     222r             ] (duplicate for contigous 211 222r)
			          Real(0.00) ,           Real(0.00) ,           Real(0.00) ,           Real(1.00) ,//180 @ z        [                  222r             ] (duplicate for contigous     222r)
			          Real(0.00) , Constants<Real>::r1_2, Constants<Real>::r1_2,           Real(0.00) ,//180 @  x,y     [                  222r   422    432]
			          Real(0.00) ,-Constants<Real>::r1_2, Constants<Real>::r1_2,           Real(0.00) ,//180 @ -x,y     [                  222r   422    432]
			Constants<Real>::r1_2,           Real(0.00) ,           Real(0.00) , Constants<Real>::r1_2,// 90 @  z       [                       4 422    432]
			Constants<Real>::r1_2,           Real(0.00) ,           Real(0.00) ,-Constants<Real>::r1_2,// 90 @ -z       [                       4 422    432]
			          Real(0.00) ,           Real(0.00) ,           Real(0.00) ,           Real(1.00) ,//180 @ z        [  112         222      4 422 23 432] {    222r}
			          Real(1.00) ,           Real(0.00) ,           Real(0.00) ,           Real(0.00) ,//identity       [1 112 121     222      4 422 23 432] {211 222r}
			          Real(0.00) ,           Real(0.00) ,           Real(1.00) ,           Real(0.00) ,//180 @ y        [      121     222        422 23 432] 
			          Real(0.00) ,           Real(1.00) ,           Real(0.00) ,           Real(0.00) ,//180 @ x        [              222        422 23 432] {211     }
			          Real(0.50) ,           Real(0.50) ,           Real(0.50) ,           Real(0.50) ,//120 @  x, y, z [                             23 432]
			          Real(0.50) ,-          Real(0.50) ,-          Real(0.50) ,-          Real(0.50) ,//120 @ -x,-y,-z [                             23 432]
			          Real(0.50) ,-          Real(0.50) ,           Real(0.50) ,           Real(0.50) ,//120 @ -x, y, z [                             23 432]
			          Real(0.50) ,           Real(0.50) ,-          Real(0.50) ,-          Real(0.50) ,//120 @  x,-y,-z [                             23 432]
			          Real(0.50) ,           Real(0.50) ,-          Real(0.50) ,           Real(0.50) ,//120 @  x,-y, z [                             23 432]
			          Real(0.50) ,-          Real(0.50) ,           Real(0.50) ,-          Real(0.50) ,//120 @ -x, y,-z [                             23 432]
			          Real(0.50) ,           Real(0.50) ,           Real(0.50) ,-          Real(0.50) ,//120 @  x, y,-z [                             23 432]
			          Real(0.50) ,-          Real(0.50) ,-          Real(0.50) ,           Real(0.50) ,//120 @ -x,-y, z [                             23 432]
			          Real(0.00) ,           Real(0.00) , Constants<Real>::r1_2, Constants<Real>::r1_2,//180 @  y,z     [                                432]
			          Real(0.00) ,           Real(0.00) ,-Constants<Real>::r1_2, Constants<Real>::r1_2,//180 @ -y,z     [                                432]
			          Real(0.00) , Constants<Real>::r1_2,           Real(0.00) , Constants<Real>::r1_2,//180 @  z,x     [                                432]
			          Real(0.00) , Constants<Real>::r1_2,           Real(0.00) ,-Constants<Real>::r1_2,//180 @ -z,x     [                                432]
			Constants<Real>::r1_2, Constants<Real>::r1_2,           Real(0.00) ,           Real(0.00) ,// 90 @  x       [                                432]
			Constants<Real>::r1_2,-Constants<Real>::r1_2,           Real(0.00) ,           Real(0.00) ,// 90 @ -x       [                                432]
			Constants<Real>::r1_2,           Real(0.00) , Constants<Real>::r1_2,           Real(0.00) ,// 90 @  y       [                                432]
			Constants<Real>::r1_2,           Real(0.00) ,-Constants<Real>::r1_2,           Real(0.00) ,// 90 @ -y       [                                432]
		};

		//hexagonal type symmetry operators
		//1, 3, 321, 312, 6, and 622 are all contigous subests (with only 3 extra values beyond 622)
		static Real const hex[(12+3)*4] = {
			Constants<Real>::r3_4,           Real(0.00) ,           Real(0.00) ,           Real(0.50) ,// 60 @  z                  [        6 622    ]
			Constants<Real>::r3_4,           Real(0.00) ,           Real(0.00) ,          -Real(0.50) ,// 60 @ -z                  [        6 622    ]
			          Real(0.00) ,           Real(0.00) ,           Real(0.00) ,           Real(1.00) ,//180 @  z                  [        6 622    ] (also in cub)
			          Real(1.00) ,           Real(0.00) ,           Real(0.00) ,           Real(0.00) ,//identity                  [1 3 321 6 622    ] (also in cub) {312}
			          Real(0.50) ,           Real(0.00) ,           Real(0.00) , Constants<Real>::r3_4,//120 @  z                  [  3 321 6 622    ]               {312}
			          Real(0.50) ,           Real(0.00) ,           Real(0.00) ,-Constants<Real>::r3_4,//120 @ -z                  [  3 321 6 622    ]               {312}
			          Real(0.00) ,           Real(1.00) ,           Real(0.00) ,           Real(0.00) ,//180 @  x                  [    321   622    ] (also in cub)
			          Real(0.00) ,           Real(0.50) , Constants<Real>::r3_4,           Real(0.00) ,//180 @ (x rotated 60 @  z) [    321   622    ]
			          Real(0.00) ,           Real(0.50) ,-Constants<Real>::r3_4,           Real(0.00) ,//180 @ (x rotated 60 @ -z) [    321   622    ]
			          Real(0.00) , Constants<Real>::r3_4,           Real(0.50) ,           Real(0.00) ,//180 @ (x rotated 30 @  z) [          622 312]
			          Real(0.00) , Constants<Real>::r3_4,          -Real(0.50) ,           Real(0.00) ,//180 @ (x rotated 30 @ -z) [          622 312]
			          Real(0.00) ,           Real(0.00) ,           Real(1.00) ,           Real(0.00) ,//180 @  y                  [          622 312] (also in cub)
			          Real(1.00) ,           Real(0.00) ,           Real(0.00) ,           Real(0.00) ,//identity                  [              312] (duplicate for 312 block)
			          Real(0.50) ,           Real(0.00) ,           Real(0.00) , Constants<Real>::r3_4,//120 @  z                  [              312] (duplicate for 312 block)
			          Real(0.50) ,           Real(0.00) ,           Real(0.00) ,-Constants<Real>::r3_4,//120 @ -z                  [              312] (duplicate for 312 block)
		};

		//get correct operators based on symmetry
		switch(rotationGroup().type) {
			case PG::   _1: return cub + 8*4;
			case PG:: _121: return cub + 8*4;
			case PG:: _112: return cub + 7*4;
			case PG:: _222: return cub + 7*4;
			case PG::_222r: return cub + 1*4;
			case PG::   _4: return cub + 5*4;
			case PG:: _422: return cub + 3*4;
			case PG::   _3: return hex + 3*4;
			case PG:: _321: return hex + 3*4;
			case PG:: _312: return hex + 9*4;
			case PG::   _6: return hex + 0*4;
			case PG:: _622: return hex + 0*4;
			case PG::  _23: return cub + 7*4;
			case PG:: _432: return cub + 3*4;
			default       : throw std::logic_error("rotationGroup() '" + rotationGroup().name() + "' not handled for rotOps()");
		}
	}

	//@brief : get the number of mirror planes
	//@return: number of mirror planes
	uint_fast8_t PointGroup::numMirror() const {
		switch(type) {
			case PG::   _1:
			case PG::  _b1:
			case PG:: _121:
			case PG:: _112:
			case PG:: _222:
			case PG::_222r:
			case PG::   _4:
			case PG::  _b4:
			case PG:: _422:
			case PG::   _3:
			case PG::  _b3:
			case PG:: _321:
			case PG:: _312:
			case PG::   _6:
			case PG:: _622:
			case PG::  _23:
			case PG:: _432: return 0;
			case PG:: _1m1:
			case PG:: _11m:
			case PG::_12m1:
			case PG::_112m:
			case PG::  _4m:
			case PG::  _b6:
			case PG::  _6m: return 1;
			case PG:: _mm2:
			case PG::_mm2r:
			case PG::_b42m:
			case PG::_b4m2: return 2;
			case PG:: _mmm:
			case PG::_mmmr:
			case PG:: _3m1:
			case PG:: _31m:
			case PG::_b3m1:
			case PG::_b31m:
			case PG:: _mb3: return 3;
			case PG:: _4mm:
			case PG::_b62m:
			case PG::_b6m2: return 4;
			case PG::_4mmm: return 5;
			case PG:: _6mm:
			case PG::_b43m: return 6;
			case PG::_6mmm: return 7;
			case PG::_mb3m: return 9;
		}
		return 0;
	}

	//@brief : get mirror plane normals
	//@return: pointer to list of mirror plane normals (as x,y,z unit vectors)
	template <typename Real>
	Real const * PointGroup::mirrors() const {
		//mirror planes for cubic type groups
		//all cubic groups are a contigous subset
		static Real const cub[11*3] = {
			-Constants<Real>::r1_2, Constants<Real>::r1_2,                0     ,//-xy [                          mmmr                                         ] (duplicate for mmr)
			 Constants<Real>::r1_2, Constants<Real>::r1_2,                0     ,// xy [                          mmmr                                         ] (duplicate for mmr)
			                0     ,                0     ,                1     ,// z  [   11m      112/m mmm     mmmr 4/m                    4/mmm m3      m3m] {          2mm m2m}
			                0     ,                1     ,                0     ,// y  [1m1    12/m1      mmm mm2               4mm      -4m2 4/mmm m3      m3m] {          2mm    }
			                1     ,                0     ,                0     ,// x  [                  mmm mm2               4mm      -4m2 4/mmm m3      m3m] {m11 2/m11     m2m}
			-Constants<Real>::r1_2, Constants<Real>::r1_2,                0     ,//-xy [                         (mmmr)    mm2r 4mm -42m      4/mmm    -43m m3m]
			 Constants<Real>::r1_2, Constants<Real>::r1_2,                0     ,// xy [                         (mmmr)    mm2r 4mm -42m      4/mmm    -43m m3m]
			                0     ,-Constants<Real>::r1_2, Constants<Real>::r1_2,//-yz [                                                               -43m m3m]
			                0     , Constants<Real>::r1_2, Constants<Real>::r1_2,// yz [                                                               -43m m3m]
			-Constants<Real>::r1_2,                0     , Constants<Real>::r1_2,//-xz [                                                               -43m m3m]
			 Constants<Real>::r1_2,                0     , Constants<Real>::r1_2,// xz [                                                               -43m m3m]
		};

		//mirror planes for hexagonal type groups
		//all 3/6 fold groups are a contigous subset
		static Real const hex[8*3] = {
			                0     ,                0     ,  1,// z  [                  -6 6/m          -62m 6/mmm]
			                0     ,                1     ,  0,// y  [    31m      -31m        6mm      -62m 6/mmm]
			 Constants<Real>::r3_4,           Real(0.5)  ,  0,// 30 [    31m      -31m        6mm      -62m 6/mmm]
			-Constants<Real>::r3_4,           Real(0.5)  ,  0,//-30 [    31m      -31m        6mm      -62m 6/mmm]
			                1     ,                0     ,  0,// x  [3m1     -3m1             6mm -6m2      6/mmm]
			           Real(0.5)  , Constants<Real>::r3_4,  0,// 60 [3m1     -3m1             6mm -6m2      6/mmm]
			-          Real(0.5)  , Constants<Real>::r3_4,  0,//-60 [3m1     -3m1             6mm -6m2      6/mmm]
			                0     ,                0     ,  1,// z  [                             -6m2           ] (duplicate for -6m2)
		};

		switch(type) {
			//no mirrors
			case PG::   _1:
			case PG::  _b1:
			case PG:: _121:
			case PG:: _112:
			case PG:: _222:
			case PG::_222r:
			case PG::   _4:
			case PG::  _b4:
			case PG:: _422:
			case PG::   _3:
			case PG::  _b3:
			case PG:: _321:
			case PG:: _312:
			case PG::   _6:
			case PG:: _622:
			case PG::  _23:
			case PG:: _432: return NULL;

			//cubic types
			case PG:: _1m1: return cub + 3*3;
			case PG:: _11m: return cub + 2*3;
			case PG::_12m1: return cub + 3*3;
			case PG::_112m: return cub + 2*3;
			case PG:: _mm2: return cub + 3*3;
			case PG::_mm2r: return cub + 5*3;
			case PG:: _mmm: return cub + 2*3;
			case PG::_mmmr: return cub + 0*3;
			case PG::  _4m: return cub + 2*3;
			case PG:: _4mm: return cub + 3*3;
			case PG::_b42m: return cub + 5*3;
			case PG::_b4m2: return cub + 3*3;
			case PG::_4mmm: return cub + 2*3;
			case PG:: _mb3: return cub + 2*3;
			case PG::_b43m: return cub + 5*3;
			case PG::_mb3m: return cub + 2*3;

			//hex types
			case PG:: _3m1: return hex + 4*3;
			case PG:: _31m: return hex + 1*3;
			case PG::_b3m1: return hex + 4*3;
			case PG::_b31m: return hex + 1*3;
			case PG::  _b6: return hex + 0*3;
			case PG::  _6m: return hex + 0*3;
			case PG:: _6mm: return hex + 1*3;
			case PG::_b6m2: return hex + 4*3;
			case PG::_b62m: return hex + 0*3;
			case PG::_6mmm: return hex + 0*3;
		}
		return NULL;
	}

	//@brief : get the number of rotational symmetry operators
	//@return: number of rotational symmetry operators
	uint_fast8_t PointGroup::numRotAxis() const {
		switch(rotationGroup().type) {
			case PG::   _1: return  0;
			case PG:: _121:
			case PG:: _112:
			case PG::   _3:
			case PG::   _4:
			case PG::   _6: return  1;
			case PG:: _222:
			case PG::_222r: return  3;
			case PG:: _321:
			case PG:: _312: return  4;
			case PG:: _422: return  5;
			case PG:: _622:
			case PG::  _23: return  7;
			case PG:: _432: return 13;
			default       : throw std::logic_error("rotationGroup() '" + rotationGroup().name() + "' not handled for numRotOps()");
		}
	}
	
	//@brief : get the rotational symmetry axis
	//@return: pointer to rotational symmetry axis (as n,x,y,z where xyz is a unit axis and n is the order (negative for rotoinversion))
	template <typename Real>
	Real const * PointGroup::rotAxis() const {
		//cubic 432 symmetry operators (1, 2, 222, 4, 422, and 23 are all a subset of these operators)
		//1, 121, 222, and 23 are in contiguous blocks as ordered
		static Real const cub[64] = {
			2, Constants<Real>::r1_2, Constants<Real>::r1_2,                     0,//2 fold @  xy
			2,-Constants<Real>::r1_2, Constants<Real>::r1_2,                     0,//2 fold @ -xy
			2, Constants<Real>::r1_2,                     0, Constants<Real>::r1_2,//2 fold @  xz
			2,-Constants<Real>::r1_2,                     0, Constants<Real>::r1_2,//2 fold @ -xz
			2,                     0, Constants<Real>::r1_2, Constants<Real>::r1_2,//2 fold @  yz
			2,                     0,-Constants<Real>::r1_2, Constants<Real>::r1_2,//2 fold @ -yz
			4,                     0,                     0,                     1,//4 fold @ z
			4,                     0,                     1,                     0,//4 fold @ y
			4,                     1,                     0,                     0,//4 fold @ x
			3, Constants<Real>::r1_3, Constants<Real>::r1_3, Constants<Real>::r1_3,//3 fold @ xyz
			3,-Constants<Real>::r1_3, Constants<Real>::r1_3, Constants<Real>::r1_3,//3 fold @ -x yz
			3,-Constants<Real>::r1_3,-Constants<Real>::r1_3, Constants<Real>::r1_3,//3 fold @ -x-yz
			3, Constants<Real>::r1_3,-Constants<Real>::r1_3, Constants<Real>::r1_3,//3 fold @  x-yz
			2,                     1,                     0,                     0,//2 fold @ x
			2,                     0,                     1,                     0,//2 fold @ y
			2,                     0,                     0,                     1,//2 fold @ z
		};

		static Real const cubI[64] = {
			 cub[4* 0+0], cub[4* 0+1], cub[4* 0+2], cub[4* 0+3],
			 cub[4* 1+0], cub[4* 1+1], cub[4* 1+2], cub[4* 1+3],
			 cub[4* 2+0], cub[4* 2+1], cub[4* 2+2], cub[4* 2+3],
			 cub[4* 3+0], cub[4* 3+1], cub[4* 3+2], cub[4* 3+3],
			 cub[4* 4+0], cub[4* 4+1], cub[4* 4+2], cub[4* 4+3],
			 cub[4* 5+0], cub[4* 5+1], cub[4* 5+2], cub[4* 5+3],
			 cub[4* 6+0], cub[4* 6+1], cub[4* 6+2], cub[4* 6+3],
			 cub[4* 7+0], cub[4* 7+1], cub[4* 7+2], cub[4* 7+3],
			 cub[4* 8+0], cub[4* 8+1], cub[4* 8+2], cub[4* 8+3],
			-cub[4* 9+0], cub[4* 9+1], cub[4* 9+2], cub[4* 9+3],
			-cub[4*10+0], cub[4*10+1], cub[4*10+2], cub[4*10+3],
			-cub[4*11+0], cub[4*11+1], cub[4*11+2], cub[4*11+3],
			-cub[4*12+0], cub[4*12+1], cub[4*12+2], cub[4*12+3],
			 cub[4*13+0], cub[4*13+1], cub[4*13+2], cub[4*13+3],
			 cub[4*14+0], cub[4*14+1], cub[4*14+2], cub[4*14+3],
			 cub[4*15+0], cub[4*15+1], cub[4*15+2], cub[4*15+3],
		};

		static Real const hex[48] = {
			6,                      0,                     0,  1,//z
			2,                      0,                     1,  0,//y
			2,  Constants<Real>::r3_4,           Real(0.5)  ,  0,// 60
			2, -Constants<Real>::r3_4,           Real(0.5)  ,  0,//-60
			2,                      1,                     0,  0,//x
			2,            Real(0.5)  , Constants<Real>::r3_4,  0,// 30
			2, -          Real(0.5)  , Constants<Real>::r3_4,  0,//-30
			3,                      0,                     0,  1,//z (duplicate so more point groups are subset of this)
			2,                      0,                     1,  0,//y (duplicate so more point groups are subset of this)
			2,  Constants<Real>::r3_4,           Real(0.5)  ,  0,// 60 (duplicate so more point groups are subset of this)
			2, -Constants<Real>::r3_4,           Real(0.5)  ,  0,//-60 (duplicate so more point groups are subset of this)
		};

		static Real const hexI[48] = {
			-hex[4* 0+0], hex[4* 0+1], hex[4* 0+2], hex[4* 0+3],
			 hex[4* 1+0], hex[4* 1+1], hex[4* 1+2], hex[4* 1+3],
			 hex[4* 2+0], hex[4* 2+1], hex[4* 2+2], hex[4* 2+3],
			 hex[4* 3+0], hex[4* 3+1], hex[4* 3+2], hex[4* 3+3],
			 hex[4* 4+0], hex[4* 4+1], hex[4* 4+2], hex[4* 4+3],
			 hex[4* 5+0], hex[4* 5+1], hex[4* 5+2], hex[4* 5+3],
			 hex[4* 6+0], hex[4* 6+1], hex[4* 6+2], hex[4* 6+3],
			-hex[4* 7+0], hex[4* 7+1], hex[4* 7+2], hex[4* 7+3],
			 hex[4* 8+0], hex[4* 8+1], hex[4* 8+2], hex[4* 8+3],
			 hex[4* 9+0], hex[4* 9+1], hex[4* 9+2], hex[4* 9+3],
			 hex[4*10+0], hex[4*10+1], hex[4*10+2], hex[4*10+3],
		};

		static Real const b62[28] = {
			 hex[4* 4+0], hex[4* 4+1], hex[4* 4+2], hex[4* 4+3],// 2 @  30
			 hex[4* 5+0], hex[4* 5+1], hex[4* 5+2], hex[4* 5+3],// 2 @ -30
			 hex[4* 6+0], hex[4* 6+1], hex[4* 6+2], hex[4* 6+3],// 2 @  x
			-hex[4* 0+0], hex[4* 0+1], hex[4* 0+2], hex[4* 0+3],// 6 @  z
			 hex[4* 1+0], hex[4* 1+1], hex[4* 1+2], hex[4* 1+3],// 2 @  60
			 hex[4* 2+0], hex[4* 2+1], hex[4* 2+2], hex[4* 2+3],// 2 @ -60
			 hex[4* 3+0], hex[4* 3+1], hex[4* 3+2], hex[4* 3+3],// 2 @  y
		};

		static Real const tet[24] = {
			cub[4*13+0], cub[4*13+1], cub[4*13+2], cub[4*13+3],//2 @  x
			cub[4*14+0], cub[4*14+1], cub[4*14+2], cub[4*14+3],//2 @  y
			cub[4* 6+0], cub[4* 6+1], cub[4* 6+2], cub[4* 6+3],//4 @  z
			cub[4* 0+0], cub[4* 0+1], cub[4* 0+2], cub[4* 0+3],//2 @  xy
			cub[4* 1+0], cub[4* 1+1], cub[4* 1+2], cub[4* 1+3],//2 @ -xy
			cub[4*15+0], cub[4*15+1], cub[4*15+2], cub[4*15+3],//2 @  z (extra for 222r)
		};

		static Real const b43[28] = {
			-cub[4* 6+0], cub[4* 6+1], cub[4* 6+2], cub[4* 6+3],//4 @  z
			-cub[4* 7+0], cub[4* 7+1], cub[4* 7+2], cub[4* 7+3],//4 @  y
			-cub[4* 8+0], cub[4* 8+1], cub[4* 8+2], cub[4* 8+3],//4 @  z
			 cub[4* 9+0], cub[4* 9+1], cub[4* 9+2], cub[4* 9+3],//3 fold @ xyz
			 cub[4*10+0], cub[4*10+1], cub[4*10+2], cub[4*10+3],//3 fold @ -x yz
			 cub[4*11+0], cub[4*11+1], cub[4*11+2], cub[4*11+3],//3 fold @ -x-yz
			 cub[4*12+0], cub[4*12+1], cub[4*12+2], cub[4*12+3],//3 fold @  x-yz
		};

		static Real const b42[20] = {
			 cub[4*13+0], cub[4*13+1], cub[4*13+2], cub[4*13+3],//2 @  x
			 cub[4*14+0], cub[4*14+1], cub[4*14+2], cub[4*14+3],//2 @  y
			-cub[4* 6+0], cub[4* 6+1], cub[4* 6+2], cub[4* 6+3],//4 @  z
			 cub[4* 0+0], cub[4* 0+1], cub[4* 0+2], cub[4* 0+3],//2 @  xy
			 cub[4* 1+0], cub[4* 1+1], cub[4* 1+2], cub[4* 1+3],//2 @ -xy
		};

		//get correct operators based on symmetry
		switch(type) {
			case PG::   _1: return NULL;
			case PG:: _1m1: return NULL;
			case PG:: _11m: return NULL;
			case PG::  _b1: return NULL;
			case PG:: _121: return cub +56;//14
			case PG:: _112: return cub +60;//15
			case PG::_12m1: return cub +56;//14
			case PG::_112m: return cub +60;//15
			case PG:: _222: return cub +52;//last 3
			case PG::_222r: return tet +12;//last 3
			case PG:: _mm2: return cub +60;//15
			case PG::_mm2r: return cub +60;//15
			case PG:: _mmm: return cub +52;//last 3
			case PG::_mmmr: return tet +12;//last 3
			case PG::   _4: return cub +24;//2
			case PG::  _b4: return b43    ;//1
			case PG::  _4m: return cub +24;//2
			case PG:: _422: return tet    ;//first 5
			case PG:: _4mm: return tet  +8;//2
			case PG::_b42m: return b42    ;//first 3
			case PG::_b4m2: return b42  +8;//last 3
			case PG::_4mmm: return tet    ;//first 5
			case PG::   _3: return hex +28;//7
			case PG::  _b3: return hexI+28;//7
			case PG:: _321: return hex +16;//4->7
			case PG:: _312: return hex +28;//7->10
			case PG:: _3m1: return hex +28;//7
			case PG:: _31m: return hex +28;//7
			case PG::_b3m1: return hexI+16;//4->7
			case PG::_b31m: return hexI+28;//7->10
			case PG::   _6: return hex    ;//first 1
			case PG::  _6m: return hex    ;//first 1
			case PG::  _b6: return hexI   ;//first 1
			case PG:: _622: return hex    ;//first 7
			case PG:: _6mm: return hex    ;//first 1
			case PG::_b62m: return b62    ;//first 4
			case PG::_b6m2: return b62 +12;//last 4
			case PG::_6mmm: return hex    ;//first 7
			case PG::  _23: return cub +36;//last 7
			case PG:: _mb3: return cubI+36;//last 7
			case PG:: _432: return cub    ;//first 13
			case PG::_b43m: return b43    ;//all
			case PG::_mb3m: return cubI   ;//first 13
		}
		return NULL;
	}

	////////////////////////////////////////////////////////////////////////
	//                        Symmetry Operations                         //
	////////////////////////////////////////////////////////////////////////

	//@brief   : check if a rodrigues vector is in the fundamental zone
	//@param ro: orientation to check (x, y, z, tan(w/2))
	//@return  : true/false if ro is/isn't in the fundamental zone
	template <typename Real>
	bool PointGroup::roInFz(Real const * const ro) const {
		switch(rotationGroup().type) {
			case PG::   _1: return PointGroup::FZ1   (ro);
			case PG:: _121: return PointGroup::FZ121 (ro);
			case PG:: _112: return PointGroup::FZ112 (ro);
			case PG:: _222: return PointGroup::FZ222 (ro);
			case PG::_222r: return PointGroup::FZ222r(ro);
			case PG::   _4: return PointGroup::FZ4   (ro);
			case PG:: _422: return PointGroup::FZ422 (ro);
			case PG::   _3: return PointGroup::FZ3   (ro);
			case PG:: _321: return PointGroup::FZ321 (ro);
			case PG:: _312: return PointGroup::FZ312 (ro);
			case PG::   _6: return PointGroup::FZ6   (ro);
			case PG:: _622: return PointGroup::FZ622 (ro);
			case PG::  _23: return PointGroup::FZ23  (ro);
			case PG:: _432: return PointGroup::FZ432 (ro);
			default       : throw std::logic_error("rotationGroup() '" + rotationGroup().name() + "' not handled for roInFz()");
		}
		return false;
	}

	//@brief   : compute the symmetric equivalent orientation in the fundamental zone
	//@param qu: orientation to compute symmetric equivalent of
	//@param fz: location to write symmetric equivalent
	template <typename Real>
	void PointGroup::fzQu(Real const * const qu, Real * const fz) const {
		//Real ro[4];//use rodrigues vector to check if we're inside the FZ
		const Real q0[4] = {qu[0], qu[1], qu[2], qu[3]};//save input in case qu and fz are the same
		const uint_fast8_t num = numRotOps();//get number of symmetry operators once
		Real const * const ops = rotOps<Real>();//get symmetry operators once
		
		std::copy(qu, qu+4, fz);
		Real quI[4];
		for(uint_fast8_t i = 0; i < num; i++) {//loop over symmetry operators
			quat::mul (ops+4*i, q0, quI);//compute operator * qu (this is for passive rotations, active would be q * O_sym)
			quat::expl(quI, quI);//select north hemisphere quat (0 <= rotation <= pi)
			if(quI[0] > fz[0]) {
				std::copy(quI, quI+4, fz);
				//this is currently commented since the unit test checks against roInFz
				// qu2ro(fz, ro);//convert to rodrigues vector (for FZ check)
				// if(roInFz(ro)) return;//stop if we're in the fundamental zone
			}
		}
	}

	//@brief    : compute the disorientation between 2 orientations
	//@param qu1: first  orientation to compute disorientation of (as w,x,y,z quaternion)
	//@param qu2: second orientation to compute disorientation of (as w,x,y,z quaternion)
	//@param dis: location to write disorientation (as w,x,y,z quaternion)
	//@note     : qu2 is assumed to belong to the same point group
	template <typename Real>
	void PointGroup::disoQu(Real const * const qu1, Real const * const qu2, Real * const dis) const {
		switch(rotationGroup().type) {
			case PG:: _432: return Diso432(qu1, qu2, dis);
			default       : return DisoQu<Real>(qu1, qu2, dis, rotOps<Real>(), numRotOps(), rotOps<Real>(), numRotOps());
		}
	}

	//@brief: compute the symmetric equivalent orientation of qu2 closest to qu1
	//@param qu1: orientation to search near
	//@param qu2: orientation to compute symmetric equivalent of
	//@param qu3: location to write symmetric equivalent of qu2 closest to qu1
	template <typename Real>
	void PointGroup::nearbyQu(Real const * const qu1, Real const * const qu2, Real * const qu3) const {
		Real work[4];
		Real maxDot = Real(-1);
		const uint_fast8_t num = numRotOps();//get number of symmetry operators
		Real const * const ops = rotOps<Real>();//get symmetry operators
		for(uint_fast8_t i = 0; i < num; i++) {//loop over symmetry operators
			quat::mul(ops + 4 * i, qu2, work);//comput3e symmetric equivalent of qu2
			const Real dot = std::fabs(quat::dot(qu1, work));//compute how close we are to qu1
				if(dot > maxDot) {//this is the closest symmetric equivalent of qu2 so far
					maxDot = dot;//update best
					std::copy(work, work + 4, qu3);//save equivalent
				}
		}
		quat::expl(qu3, qu3);//restrict rotation to [0,pi]
	}

	//@brief   : reduce a unit direction to the fundamental sector
	//@param n : unit direction to reduce (magnitude is assumed to be 1)
	//@param fs: location to write symmetric equivalent of n in the fundamental sector (can be the same as n)
	//@param   : true if n was already in the fundamental sector
	template <typename Real>
	bool PointGroup::fsDir(Real const * const n, Real * const fs) const {
		std::copy(n, n+3, fs);
		switch(type) {
			case PG::   _1: return true;
			case PG::  _b1: return fs::b1   (fs);
			case PG:: _121: return fs::_121 (fs);
			case PG:: _112: return fs::_112 (fs);
			case PG:: _1m1: return fs::_1m1 (fs);
			case PG:: _11m: return fs::_11m (fs);
			case PG::_12m1: return fs::_12m1(fs);
			case PG::_112m: return fs::_112m(fs);
			case PG:: _222: return fs::_222 (fs);
			case PG::_222r: return fs::_222r(fs);
			case PG:: _mm2: return fs:: mm2 (fs);
			case PG::_mm2r: return fs:: mm2r(fs);
			case PG:: _mmm: return fs:: mmm (fs);
			case PG::_mmmr: return fs:: mmmr(fs);
			case PG::   _4: return fs::_4   (fs);
			case PG::  _b4: return fs::b4   (fs);
			case PG::  _4m: return fs::_4m  (fs);
			case PG:: _422: return fs::_422 (fs);
			case PG:: _4mm: return fs::_4mm (fs);
			case PG::_b42m: return fs::b42m (fs);
			case PG::_b4m2: return fs::b4m2 (fs);
			case PG::_4mmm: return fs::_4mmm(fs);
			case PG::   _3: return fs::_3   (fs);
			case PG::  _b3: return fs::b3   (fs);
			case PG:: _321: return fs::_321 (fs);
			case PG:: _312: return fs::_312 (fs);
			case PG:: _3m1: return fs::_3m1 (fs);
			case PG:: _31m: return fs::_31m (fs);
			case PG::_b3m1: return fs::b3m1 (fs);
			case PG::_b31m: return fs::b31m (fs);
			case PG::   _6: return fs::_6   (fs);
			case PG::  _b6: return fs::b6   (fs);
			case PG::  _6m: return fs::_6m  (fs);
			case PG:: _622: return fs::_622 (fs);
			case PG:: _6mm: return fs::_6mm (fs);
			case PG::_b6m2: return fs::b6m2 (fs);
			case PG::_b62m: return fs::b62m (fs);
			case PG::_6mmm: return fs::_6mmm(fs);
			case PG::  _23: return fs::_23  (fs);
			case PG:: _mb3: return fs:: mb3 (fs);
			case PG:: _432: return fs::_432 (fs);
			case PG::_b43m: return fs::b43m (fs);
			case PG::_mb3m: return fs::mb3m (fs);
		}
		return true;
	}

	//@param n  : crystal direction to color to color (magnitude is assumed to be 1)
	//@param rgb: location to write color [0,1]
	//@param h2r: hsl2rgb like coloring function to use with h, s, and l in [0,1] and output as [0,1] rgb: void(Real const * const hsl, Real * const rgb)
	template <typename Real>
	void PointGroup::ipfColor(Real const * const n, Real * const rgb, std::function<void(Real const*const, Real*const)> h2r) const {
		//build a list of commonly needed ipf color corners
		static const Real nz  [3]    = {           Real(0.0)  ,           Real(0.0)  ,           Real(-1)   };
		static const Real y12 [3]    = { Constants<Real>::r3_4,           Real(0.5)  ,           Real( 0)   };// 30 degrees
		static const Real y8  [3]    = { Constants<Real>::r1_2, Constants<Real>::r1_2,           Real( 0)   };// 45 degrees
		static const Real n8  [3]    = {-Constants<Real>::r1_2, Constants<Real>::r1_2,           Real( 0)   };//135 degrees
		static const Real y6  [3]    = {           Real(0.5)  , Constants<Real>::r3_4,           Real( 0)   };// 60 degrees
		static const Real y3  [3]    = {-          Real(0.5)  , Constants<Real>::r3_4,           Real( 0)   };//120 degrees
		static const Real  xy [3]    = { Constants<Real>::r1_3, Constants<Real>::r1_3, Constants<Real>::r1_3};// 111
		static const Real nxy [3]    = {-Constants<Real>::r1_3, Constants<Real>::r1_3, Constants<Real>::r1_3};//-111
		static const Real xz  [3]    = { Constants<Real>::r1_2,           Real(0.0)  , Constants<Real>::r1_2};// 101
		static const Real dirs[4][3] = {
			Real(0), Real( 0), Real(1),// z
			Real(0), Real(-1), Real(0),//-y
			Real(1), Real( 0), Real(0),// x
			Real(0), Real( 1), Real(0),// y
		};

		//triangles (when there are 2 or more rotation axis)
		static const detail::SphericalTriangle<Real> b3 (dirs[0], dirs[2], y3     );//120 degrees
		static const detail::SphericalTriangle<Real> c2 (dirs[0], dirs[2], dirs[3]);// 90 degrees
		static const detail::SphericalTriangle<Real> c2b(dirs[0], y8     , n8     );// c2 rotated 45 @ z
		static const detail::SphericalTriangle<Real> c3 (dirs[0], dirs[2], y6     );// 60 degrees
		static const detail::SphericalTriangle<Real> c3b(dirs[0], y12    , dirs[3]);// c3 rotated 30 @ z
		static const detail::SphericalTriangle<Real> c4 (dirs[0], dirs[2], y8     );// 45 degrees
		static const detail::SphericalTriangle<Real> c6 (dirs[0], dirs[2], y12    );// 30 degrees
		static const detail::SphericalTriangle<Real> m3 (dirs[0], xy     , nxy    );//cubic low triangle
		static const detail::SphericalTriangle<Real> m3m(dirs[0], xz     , xy     );//cubic triangle

		//wedges (when there is only 1 rotation axis)
		static const detail::SphericalWedge<  Real> w2 (dirs[2], dirs[3]);//2 fold wedge
		static const detail::SphericalWedge<  Real> w2b(y8     , n8     );//2 fold wedge
		static const detail::SphericalPatch<4,Real> w2c(dirs   , xz     );//alternate 2 fold wedge
		static const detail::SphericalWedge<  Real> w4 (dirs[2], y8     );//4 fold wedge
		static const detail::SphericalWedge<  Real> w3a(dirs[2], y6     );//3 fold wedge
		static const detail::SphericalWedge<  Real> w3b(y12    , dirs[3]);//w3a rotated 30 @ z
		static const detail::SphericalWedge<  Real> w6 (dirs[2], y12    );//6 fold wedge

		//@brief    : convert from unit direction to fractional hsl
		//@param n  : unit direction to convert to color [in place]
		//@param inv: true if inversion symmetry should be applied by doubling the polar angle
		static const std::function<void(Real*const,const bool)> n2hsl = [](Real*const n, const bool inv) {
				//convert from unit direction to fractional spherical coordinates
				Real theta = std::atan2(n[1], n[0]) / Constants<Real>::pi2;//[-0.5,0.5]
				if(std::signbit(theta)) theta += Real(1);//[0,1]
				Real phi = std::acos(n[2]) / Constants<Real>::pi;//[0,1]
				if(inv) phi *= Real(2);//inversion via doubling polar angle

				//convert from fractional sphercal coordinate to fractional hs;
				n[0] = theta        ;//h
				n[1] = Real(1)      ;//s
				n[2] = Real(1) - phi;//l
		};

		//compute ipf color
		fsDir<Real>(n, rgb);//temporarily store fundamental sector direction in rgb
		bool wCen = true, iTht = false, inv = false;
		switch(type) {
			case PG::  _b1: inv = true;
			case PG::   _1:
			case PG:: _11m: n2hsl(rgb, inv); return h2r(rgb, rgb);

			case PG:: _121: wCen = fs::_m11(rgb);
			                return w2c.toColor(rgb, rgb, h2r, wCen);

			case PG:: _112: wCen = fs::_m11(rgb);
			case PG:: _mm2: return w2 .toColor(rgb, rgb, h2r, wCen);

			case PG:: _1m1: std::rotate(rgb, rgb+2, rgb+3); n2hsl(rgb, inv); return h2r(rgb, rgb);

			case PG::_12m1:
			case PG::_112m:
			case PG:: _222: wCen = fs::mmm (rgb);
			case PG:: _mmm: return c2 .toColor(rgb, rgb, h2r, wCen);

			case PG::_222r: wCen = fs::mmmr(rgb);
			case PG::_mmmr: return c2b.toColor(rgb, rgb, h2r, wCen);

			case PG::_mm2r: return w2b.toColor(rgb, rgb, h2r, wCen);

			case PG::   _4: wCen = fs::_4mm(rgb);
			case PG:: _4mm: return w4 .toColor(rgb, rgb, h2r, wCen);

			case PG::  _b4: inv = true;
							if(std::signbit(rgb[2])) {
								rgb[2] = -rgb[2];
								if(std::signbit(rgb[0])) {
									rgb[0] = -rgb[0];
									rgb[1] = -rgb[1];
								}
							} else {
								std::swap(rgb[0], rgb[1]);
								rgb[1] = -rgb[1];
							}
							w2c.toHemi(rgb, rgb[0], rgb[2]);
							rgb[1] = Real(1); 
							rgb[2] = Real(1) - rgb[2] * 2;
							return h2r(rgb, rgb);

			case PG::  _4m:
			case PG:: _422: wCen = fs::_4mm(rgb);
			                return c4 .toColor(rgb, rgb, h2r, wCen);

			case PG::_b42m: wCen = !std::signbit(rgb[0]); fs::_4mm(rgb);
			                return c4 .toColor(rgb, rgb, h2r, wCen);

			case PG::_b4m2: wCen = fs::_4mm(rgb) && wCen;
			case PG::_4mmm: return c4 .toColor(rgb, rgb, h2r, wCen);

			case PG::   _3: wCen = fs::_31m(rgb);
			case PG:: _31m: return w3a.toColor(rgb, rgb, h2r, wCen);

			case PG::  _b3: b3.toHemi(rgb, rgb[0], rgb[2]);
						    rgb[1] = Real(1); 
						    rgb[2] = Real(1) - rgb[2] * 2;
						    return h2r(rgb, rgb);

			case PG:: _321:
			case PG::  _b6: wCen = fs::_31m(rgb);
			case PG::_b62m: return c3 .toColor(rgb, rgb, h2r, wCen);

			case PG:: _312: wCen = fs::_m11(rgb);
			case PG::_b6m2: return c3b.toColor(rgb, rgb, h2r, wCen, true);

			case PG:: _3m1: return w3b.toColor(rgb, rgb, h2r, wCen, true);

			case PG::_b3m1: wCen = fs::_31m(rgb); fs::_6mm(rgb);
			                return c6 .toColor(rgb, rgb, h2r, wCen);

			case PG::_b31m: iTht = true;
			case PG::  _6m:
			case PG:: _622: wCen = fs::_6mm(rgb);
			case PG::_6mmm: return c6 .toColor(rgb, rgb, h2r, wCen, iTht);

			case PG::   _6: wCen = fs::_6mm(rgb);
			case PG:: _6mm: return w6 .toColor(rgb, rgb, h2r, wCen);

			case PG::  _23: wCen = fs::mm2r(rgb);
			case PG::_b43m: return m3 .toColor(rgb, rgb, h2r, wCen);

			case PG:: _mb3:
			case PG:: _432: wCen = fs::_4mm(rgb);
			case PG::_mb3m: return m3m.toColor(rgb, rgb, h2r, wCen);
		}
	}


	////////////////////////////////////////////////////////////////////////
	//            Static Functions for Symmetry Specific Code             //
	////////////////////////////////////////////////////////////////////////

	namespace detail {
		//@brief: check if |a * b| <= c handling a == 0 and b == inf
		//@param a: a (may be 0)
		//@param b: b (may be inf)
		//@param c: c (assumed to be >= 0)
		//@return : |a * b| <= c with 0 * inf == 0 instead of nan
		template <typename Real>
		bool cmpInf(const Real a, const Real b, const Real c) {
			return std::isinf(b) ? 0 == a : a * b <= c;
		}
	}
	
	//@brief   : check if an orientation is in the fundamental zone
	//@param ro: orientation to check as [x, y, z, tan(w/2)]
	//@return  : true/false if the orientation is/isn't in the fundamental zone
	template <typename Real> bool PointGroup::FZ121 (Real const * const ro) {return FZ1(ro) && detail::cmpInf( std::fabs(ro[1]), ro[3], Real(1));}// | ro(y) | <= 1 (correcting for 0 * inf == nan)
	template <typename Real> bool PointGroup::FZ112 (Real const * const ro) {return FZ1(ro) && detail::cmpInf( std::fabs(ro[2]), ro[3], Real(1));}// | ro(z) | <= 1 (correcting for 0 * inf == nan)
	template <typename Real> bool PointGroup::FZ222 (Real const * const ro) {return FZ1(ro) && detail::cmpInf( std::max( std::fabs(ro[0]), std::max( std::fabs(ro[1]), std::fabs(ro[2]) ) ), ro[3], Real(1));}// max(ro) <= 1
	template <typename Real> bool PointGroup::FZ222r(Real const * const ro) {return FZ112(ro) && detail::cmpInf( std::fabs(ro[0]) + std::fabs(ro[1]), ro[3], Constants<Real>::r2);}
	template <typename Real> bool PointGroup::FZ3   (Real const * const ro) {return FZ1(ro) && detail::cmpInf( std::fabs(ro[2]), ro[3], Constants<Real>::r1_3);}//z cutting planes at +/- tan(pi/6)
	template <typename Real> bool PointGroup::FZ4   (Real const * const ro) {return FZ1(ro) && detail::cmpInf( std::fabs(ro[2]), ro[3], Constants<Real>::tp_8);}//z cutting planes at +/- tan(pi / (2 * N))
	template <typename Real> bool PointGroup::FZ422 (Real const * const ro) {return FZ4(ro) && FZ222(ro) & FZ222r(ro);}
	template <typename Real> bool PointGroup::FZ321 (Real const * const ro) {return FZ3(ro) && FZ3xy(ro, false);}
	template <typename Real> bool PointGroup::FZ312 (Real const * const ro) {return FZ3(ro) && FZ3xy(ro, true );}
	template <typename Real> bool PointGroup::FZ6   (Real const * const ro) {return FZ1(ro) && detail::cmpInf(std::fabs(ro[2]), ro[3], Constants<Real>::tp_12);}//z cutting planes at +/- tan(pi / (2 * N))
	template <typename Real> bool PointGroup::FZ622 (Real const * const ro) {return FZ6(ro) && FZ3xy(ro, false) & FZ3xy(ro, true);}

	//@brief   : check if an orientation is in the fundamental zone
	//@param ro: orientation to check as [x, y, z, tan(w/2)]
	//@return  : true/false if the orientation is/isn't in the fundamental zone
	//@note    : 23 group
	template <typename Real>
	bool PointGroup::FZ23 (Real const * const ro) {
		if(std::signbit(ro[3])) return false;//rotation <= pi
		const Real xyz[3] = {std::fabs(ro[0]), std::fabs(ro[1]), std::fabs(ro[2])};//element wise absolute value
		return detail::cmpInf( std::accumulate(xyz, xyz+3, Real(0)), ro[3], Real(1));//check against octohedron with faces at tan(pi/8)
	}

	//@brief   : check if an orientation is in the fundamental zone
	//@param ro: orientation to check as [x, y, z, tan(w/2)]
	//@return  : true/false if the orientation is/isn't in the fundamental zone
	//@note    : 432 group
	template <typename Real>
	bool PointGroup::FZ432(Real const * const ro) {
		if(std::signbit(ro[3])) return false;//rotation <= pi
		const Real xyz[3] = {std::fabs(ro[0]), std::fabs(ro[1]), std::fabs(ro[2])};//element wise absolute value
		const Real vMax = *std::max_element(xyz, xyz+3);//get largest absolute value element
		if(vMax * ro[3] > Constants<Real>::tp_8) return false;//check against cube faces of cube octohedron at tan(pi/8)
		return detail::cmpInf( std::accumulate(xyz, xyz+3, Real(0)), ro[3], Real(1));//check against octohedral faces cube octohedron at tan(pi/6)
	}

	//@brief   : test x and y against 60 degree ro fz cutting plane
	//@param ro: orientation to check as [x, y, z, tan(w/2)]
	//@param yx: true to swap x and y
	//@return  : true if inside fundamental zone
	template <typename Real>
	bool PointGroup::FZ3xy (Real const * const ro, const bool yx) {
		const Real x2  = std::fabs(ro[yx ? 1 : 0] * ro[3]) / 2;
		const Real y34 = std::fabs(ro[yx ? 0 : 1] * ro[3]) * Constants<Real>::r3_4;
		return y34 > x2 ? x2 + y34 < Real(1) : x2 <= Real(0.5);
	}

	//@brief    : compute the disorientation between 2 orientations
	//@param qu1: first  orientation to compute disorientation of
	//@param qu2: second orientation to compute disorientation of
	//@param dis: location to write disorientation
	//@param op1: rotational symmetry operators as quaternions (wxyz) for the first  orientation
	//@param no1: number of rotational symmetry operators for the first  orientation
	//@param op2: rotational symmetry operators as quaternions (wxyz) for the second orientation
	//@param no2: number of rotational symmetry operators for the second orientation
	template <typename Real>
	void PointGroup::DisoQu(Real const * const qu1, Real const * const qu2, Real * const dis, Real const * const op1, const uint_fast8_t no1, Real const * const op2, const uint_fast8_t no2) {
		Real delQ[4], work[4];//qu1 * conjg(qu2) iDeltaQ[4];
		quat::conj(qu2, work);//store conj(qu2) in work
		quat::mul (qu1, work, delQ);//compute qu1 * conjg(qu2)
		std::copy(delQ, delQ+4, dis);//start with the result for no symmetry
		//this loop (O_sym * qu1 * conj(qu2) * conj(O_sym)) is for passive rotations
		//active rotations are qu1 * O_sym * conj(O_sym) * conj(qu2) ==> qu1 * O_sym * conj(qu2)
		for(uint_fast8_t i = 0; i < no1; i++) {//loop over qu1 symmetry
			quat::mul (op1+4*i, delQ, work);//compute (symmetric equivalent of qu1) * conj(qu2)
			for(uint_fast8_t j = 0; j < no2; j++) {//loop over qu2 symmetry
				quat::mul (work, op2+4*j, work);//compute (symmetric equivalent of qu1) * (symmetric equivalent of conj(qu2))
				quat::expl(work, work);//restrict rotation to [0,pi]
				if(work[0] > dis[0]) std::copy(work, work+4, dis);//save the result if it is the lowest rotation angle so far
			}
		}
	}

	//@brief    : compute the disorientation between 2 cubic (432) orientations
	//@param qu1: first  orientation to compute disorientation of
	//@param qu2: second orientation to compute disorientation of
	//@param dis: location to write disorientation
	//@note     : this shortcut is significantly faster than the general algorithm
	template <typename Real> 
	void PointGroup::Diso432(Real const * const qu1, Real const * const qu2, Real * const dis) {
		//compute 2 potential rotations once
		static const Real qu3[4] = {          Real(0.5)  ,-          Real(0.5)  , -Real(0.5), -Real(0.5)};
		static const Real qu4[4] = {Constants<Real>::r1_2,-Constants<Real>::r1_2,  Real(0.0),  Real(0.0)};

		//compute initial misorientation
		quat::conj(qu2, dis);//store conj(qu2) in output temporarily
		quat::mul (qu1, dis, dis);//compute \Delta{qu} = qu1 * qu2^-1
		quat::cAbs(dis, dis);//element wise absolute value of misorientation
		std::sort(dis, dis+4, std::greater<Real>());//sort misorientation from largest to smallest (from 2 fold axis)

		//now check for improvements from 3 and 4 fold axis
		const Real s10 =  dis[1] + dis[0];
		const Real op3 = (dis[3] + dis[2] + s10) * qu3[0];//rotate by 3 fold axis
		const Real op4 =                    s10  * qu4[0];//rotate by 4 fold axis
		if(op3 >= op4 && op3 > dis[0]) {//3 fold is best
			quat::mul (dis, qu3, dis);//apply 3 fold axis
			quat::cAbs(dis, dis);//element wise abs
			std::sort(dis, dis+4, std::greater<Real>());//sort large -> small
		} else if(op4 > op3 && op4 > dis[0]) {//4 fold is best
			quat::mul (dis, qu4, dis);//apply 4 fold axis
			quat::cAbs(dis, dis);//element wise abs
			std::sort(dis, dis+4, std::greater<Real>());//sort large -> small
		}
	}
}

#endif//_symmetry_h_
