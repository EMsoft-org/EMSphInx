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

#ifndef _position_h_
#define _position_h_

#include <string>
#include <vector>

namespace xtal {
	//bitmasks for wyckoff position string formatting
	namespace wyckoff {
		enum Format {
			SftX   = 1,//should space be left for x translations e.g. "x    ,y,z+1/2" instead of "x,y,z+1/2"
			SftY   = 2,//should space be left for y translations
			SftZ   = 4,//should space be left for z translations
			SftXYZ = 7,//SftX | SftY | SftZ
			Hex    = 8,//should space be left for 3 fold rotations e.g. "  -x,y-x,z" instead of "-x,y-x,z"
		};
	}

	//@brief: helper to compactly hold a 4x4 general position matrix
	struct GenPos {

		////////////////////////////////////////////////
		// methods to build a general position matrix //
		////////////////////////////////////////////////

		//@brief: construct a default general position (identity matrix)
		GenPos() : GenPos(0) {}

		//@brief  : construct a GenPos from a set of bits
		//@param b: bits to construct from
		//@note   : you'd better know exactly what you're doing with the internal structure to use this
		GenPos(const uint32_t b) {u.i = b;}

		//@brief: construct an identity matrix
		static GenPos Identity() {return GenPos(0x0000);}

		//@brief: construct an inversion matrix
		static GenPos Inversion() {return GenPos(0x0007);}

		//@brief    : construct a mirror plane with the given normal
		//@param n  : mirror plane normal
		//@param hex: true for 3/6 fold type groups, e.g. 2@100 == -x,y,z for false, -x+y,y,z for true
		//@return   : mirror matrix e.g. x mirror for (1,0,0)
		static GenPos Mirror(int8_t const * const n, const bool hex = false);

		//@brief    : construct a two fold rotation axis about the given direction
		//@param n  : two fold axis
		//@param hex: true for 3/6 fold type groups, e.g. 2@100 == x,-y,-z for false, x-y,-y,-z for true
		//@return   : two fold rtation matrix
		static GenPos Two(int8_t const * const n, const bool hex = false);

		//@brief    : construct a three fold rotation axis about the given direction
		//@param n  : three fold axis (must by 001 or 111 type direction)
		//@param inv: true for roto inversion (-3)
		//@return   : three fold rtation matrix
		static GenPos Three(int8_t const * const n, const bool inv = false);

		//@brief  : construct an N fold rotation about the Z axis
		//@param n: rotational order (must be +/- 1,2,3,4, or 6)
		//@return : rotation matrix
		//@note   : 1, -1 and 2 return idenity, inversion, and mirror respectively
		static GenPos Z(const int8_t n);

		////////////////////////////////////////////////
		//       3x3 + augmenting vector access       //
		////////////////////////////////////////////////

		//@brief : get the upper left 3x3 sub matrix
		//@return: pointer to matrix in row major order
		//@note  : read only
		int8_t const * getMat3() const {return IdxToMat3(get3x3());}

		//@brief  : set the upper left 3x3 sub matrix
		//@param m: matrix to use (must be one of the 64 valid 3x3 position matrices)
		void setMat3(int8_t const * m);

		//@brief  : get the translation in [0,1) * 24
		//@param t: location to write translation {a, b, c}
		int8_t const * getTrans() const {return u.c+1;}
		void getTrans(int8_t * t) const {std::copy(u.c+1, u.c+4, t);}

		//@brief  : set the translation * 24
		//@param t: translation for x,y,z each in [0,1) * 24
		void setTrans(int8_t* t);

		//@brief: set translation to {0,0,0}
		void removeTrans() {std::fill(u.c+1, u.c+4, 0);}

		//@brief  : get the translation
		//@param t: location to write translation {a, b, c}
		template <typename Real> void getTransReal(Real * const t) const;

		//@brief  : set the translation
		//@param t: translation to set {a, b, c}
		//@note   : translations must be one of (0, 1/6, 1/4, 1/3, 1/2, 2/3, 3/4, or 5/6) +/- n*1 for normal
		template <typename Real> void setTransReal(Real const * const t);

		//@brief   : transform a hexagonal general position matrix from a 3 fold coordinate system to a cartesian frame
		//@param om: location to write 3x3 matrix in cartesian frame
		//@note    : throws for cubic matricies
		template <typename Real> void getMat3HexCart(Real * const om);

		////////////////////////////////////////////////
		//              property queries              //
		////////////////////////////////////////////////

		//@brief : get the determinant of the 3x3 matrix
		//@return: determinant
		int_fast8_t det() const;

		//@brief : get the trace of the 3x3 matrix
		//@return: trace
		int_fast8_t tr() const;

		//@brief : get the rotational order of the 3x3 matrix
		//@return: n * det() such that matrix^n == I, must be +/- 1, 2, 3, 4, 6
		//@note  : e.g. 1 for identity, -2 for mirror plane, 3 for 3 fold rotation etc
		int_fast8_t order() const;

		//@brief : get the axis associated with the 3x3 matrix (i.e. rotation axis or mirror plane normal)
		//@return: rotation axis as e.g. {2,1,0}
		//@note  : identity and inversion return {1,0,0}
		int8_t const* axis() const;

		//@brief : determine if this matrix contains a translation component
		//@return: true if there is a translation, false if the translation is {0,0,0}
		//@note  : in conjunction with order this enables all possible element types to be distinguished
		//         switch(order()) {
		//           case  1      : hasTranslation() ? "translation" : "identity"
		//           case  2,3,4,6: hasTranslation() ? "screw"       : "rotation"
		//           case -1      : hasTranslation() ? ""            : "inversion"
		//           case -2      : hasTranslation() ? "glide"       : "mirror"
		//           case -3,-4,-6: hasTranslation() ? ""            : "rotoinversion"
		//         }
		bool hasTranslation() const {return 0 != getTrans()[0] || 0 != getTrans()[1] || 0 != getTrans()[2];}

		//@brief : check if a general position is identity
		//@return: true/false if this does/doesn't represent identity
		bool isIdent() const {return *this == Identity();}

		////////////////////////////////////////////////
		//             matrix operations              //
		////////////////////////////////////////////////

		//@brief    : multiply this general position with another
		//@param rhs: other general position to multiply with
		//@return   : this
		GenPos& operator*=(const GenPos& rhs);

		//@brief    : multiply two position with another
		//@param rhs: other general position to multiply with
		//@return   : GenPos(this * rhs)
		GenPos operator*(const GenPos& rhs) const {return GenPos(*this) *= rhs;}

		//@brief    : get the inverse of this general position
		//@return   : this^-1
		GenPos inverse() const;

		//@brief  : apply a transformation matrix to this general position
		//@param p: transformation matrix
		//@note   : applys basis transform A' = P^-1 * A * P for P == I | p
		GenPos transform(const GenPos& p) const {return p.inverse().operator*=(*this) * p;}

		//@brief  : shift the origin by a given vector
		//@param p: origin translation for x,y,z each in [0,1) * 24
		//@note   : this applys a special case of basis transform A' = P^-1 * A * P for P == I | -p
		GenPos shiftOrigin(int8_t const*const p) const;

		//@brief    : since general positions can be described by a single int comparison is easy
		//@param rhs: other general position to compare against
		//@return   : this < rhs
		bool operator<(const GenPos& rhs) const {return u.i < rhs.u.i;}

		//@brief    : since general positions can be described by a single int comparison is easy
		//@param rhs: other general position to compare against
		//@return   : this == rhs
		bool operator==(const GenPos& rhs) const {return u.i == rhs.u.i;}

		//@brief    : since general positions can be described by a single int comparison is easy
		//@param rhs: other general position to compare against
		//@return   : this == rhs
		bool operator!=(const GenPos& rhs) const {return u.i != rhs.u.i;}

		//@brief    : build closed set of general positions
		//@param gen: set of general positions to close
		//@return   : closed set of general positions
		static std::vector<GenPos> CloseSet(std::vector<GenPos> gen);

		////////////////////////////////////////////////
		//                     IO                     //
		////////////////////////////////////////////////

		//@brief    : convert to string representation
		//@param pre: prefix for each line of matrix
		//@return   : string representation of 4x3 general position matrix
		std::string to_string(std::string pre = "") const;

		//@brief    : convert to string representation
		//@param fmt: bitmask of format options
		//@return   : string representation of 4x3 general position matrix
		std::string to_wyckoff(const int fmt = wyckoff::SftX | wyckoff::SftY | wyckoff::SftZ | wyckoff::Hex) const;

		//@brief    : build general position matrix from a piece of an EMsoft generator string
		//@param str: 4 character EMsoft generator string e.g. "eBFF"
		void fromEMsoft(char const * str);

		//@brief    : convert stored position matrix to EMsoft generator string if possible
		//@param str: location to write 4 character EMsoft generator string e.g. "eBFF"
		void toEMsoft(char * str, const bool ori = false) const;

		////////////////////////////////////////////////
		//  detail (exposed primariily for testing)   //
		////////////////////////////////////////////////

		//@brief  : get a 3x3 position matrix from an index
		//@param i: index of matrix to get, must be [0,64)
		//@return : pointer to 3x3 matrix corresponding to i (row major)
		static int8_t const * IdxToMat3(const int_fast8_t i);

		//@brief  : get the index of a 3x3 position matrix
		//@param m: 3x3 matrix to get index of
		//@return : index of matrix m (such that IdxToMat3(index) has the same matrix as m)
		//@note   : returns 0xFF if m is invalid (not in 64 possible returns of IdxToMat3)
		//@note   : entries of m should be exclusively -1, 0, or 1
		static uint_fast8_t Mat3ToIdx(int8_t const * const m);

		private:
			//a general position is a augmented 3x3 matrix
			//there are only 64 crystallographically valid 3x3 matrix types
			//there are only 12 crystallographically valid translations (all 24ths)
			//that means we need at least log(64 * 12^3, 2) == 17 bits to store any possible matrix
			union {
				uint32_t i   ;
				int8_t   c[4];//index of matrix + 3x translation * 24 - matrix must be [0,64), translations should be [0,24) {really [0,21]}
			} u;

			//@brief : get the 6 bits representing the 3x3 matrix
			//@return: 3x3 matrix index [0,64)
			int_fast8_t get3x3() const {return u.c[0];}

			//@brief  : set the 6 bits representing the 3x3 matrix
			//@param i: 3x3 matrix index [0,64)
			void set3x3(const int_fast8_t i) {u.c[0] = i;}
	};
}

////////////////////////////////////////////////////////////////////////
//                          Implementations                           //
////////////////////////////////////////////////////////////////////////

#include <algorithm>
#include <cmath>
#include <set>
#include <sstream>
#include <iomanip>
#include <stdexcept>

namespace xtal {

	////////////////////////////////////////////////////////////////////////
	//                          Helper Functions                          //
	////////////////////////////////////////////////////////////////////////

	namespace detail {

		//@brief    : get string representation of a 24th fraction
		//@param num: numerator (must be in [0,24])
		//@param uni: should unicode characters be used to make nicer fractions
		//@return   : string representation of num/24
		std::string get24th(const size_t num, const bool uni) {
			static const bool vul = false;//if unicode characters are used should vulgar fractions be used where possible
			switch(num) {
				case  0: return uni ? ( vul ? "0" :  " 0 " ) :  " 0 " ;
				case  1: return uni ? (              "¹⧸₂₄") : "1/24" ;
				case  2: return uni ? (              "¹⧸₁₂") : "1/12" ;
				case  3: return uni ? ( vul ? "⅛" :  "¹⧸₈" ) :  "1/8" ;
				case  4: return uni ? ( vul ? "⅙" :  "¹⧸₆" ) :  "1/6" ;
				case  5: return uni ? (              "⁵⧸₂₄") : "5/24" ;
				case  6: return uni ? ( vul ? "¼" :  "¹⧸₄" ) :  "1/4" ;
				case  7: return uni ? (              "⁷⧸₂₄") : "7/24" ;
				case  8: return uni ? ( vul ? "⅓" :  "¹⧸₃" ) :  "1/3" ;
				case  9: return uni ? ( vul ? "⅜" :  "³⧸₈" ) :  "3/8" ;
				case 10: return uni ? (              "⁵⧸₁₂") : "5/12" ;
				case 11: return uni ? (             "¹¹⧸₂₄") : "11/24";
				case 12: return uni ? ( vul ? "½" :  "¹⧸₂" ) :  "1/2" ;
				case 13: return uni ? (             "¹³⧸₂₄") : "13/24";
				case 14: return uni ? (              "⁷⧸₁₂") : "7/12" ;
				case 15: return uni ? ( vul ? "⅝" :  "⁵⧸₈" ) :  "5/8" ;
				case 16: return uni ? ( vul ? "⅔" :  "²⧸₃" ) :  "2/3" ;
				case 17: return uni ? (             "¹⁷⧸₂₄") : "17/24";
				case 18: return uni ? ( vul ? "¾" :  "³⧸₄" ) :  "3/4" ;
				case 19: return uni ? (             "¹⁹⧸₂₄") : "19/24";
				case 20: return uni ? ( vul ? "⅚" :  "⁵⧸₆" ) :  "5/6" ;
				case 21: return uni ? ( vul ? "⅞" :  "⁷⧸₈" ) :  "7/8" ;
				case 22: return uni ? (             "¹¹⧸₁₂") : "11/12";
				case 23: return uni ? (             "²³⧸₂₄") : "23/24";
				default: throw std::runtime_error("get24th argument must lie in [0,24)");
			}
		}
	}

	////////////////////////////////////////////////////////////////////////
	//                      General Position Members                      //
	////////////////////////////////////////////////////////////////////////

	////////////////////////////////////////////////
	// methods to build a general position matrix //
	////////////////////////////////////////////////

	//@brief    : construct a mirror plane with the given normal
	//@param n  : mirror plane normal
	//@param hex: true for 3/6 fold type groups, e.g. 2@100 == -x,y,z for false, -x+y,y,z for true
	//@return   : mirror matrix e.g. x mirror for (1,0,0)
	GenPos GenPos::Mirror(int8_t const * const n, const bool hex) {
		if     ( ( 1 == n[0] &&  0 == n[1] &&  0 == n[2]) || (-1 == n[0] &&  0 == n[1] &&  0 == n[2]) ) return GenPos(hex ? 0x0036 : 0x0001);//x
		else if( ( 0 == n[0] &&  1 == n[1] &&  0 == n[2]) || ( 0 == n[0] && -1 == n[1] &&  0 == n[2]) ) return GenPos(hex ? 0x0038 : 0x0002);//y
		else if( ( 0 == n[0] &&  0 == n[1] &&  1 == n[2]) || ( 0 == n[0] &&  0 == n[1] && -1 == n[2]) ) return GenPos(0x0004);//z   (hex and cubic)
		else if( ( 0 == n[0] &&  1 == n[1] && -1 == n[2]) || ( 0 == n[0] && -1 == n[1] &&  1 == n[2]) ) return GenPos(0x0018);//y-z
		else if( ( 0 == n[0] &&  1 == n[1] &&  1 == n[2]) || ( 0 == n[0] && -1 == n[1] && -1 == n[2]) ) return GenPos(0x001E);//y+z
		else if( ( 1 == n[0] && -1 == n[1] &&  0 == n[2]) || (-1 == n[0] &&  1 == n[1] &&  0 == n[2]) ) return GenPos(0x0020);//x-y (hex and cubic)
		else if( ( 1 == n[0] &&  1 == n[1] &&  0 == n[2]) || (-1 == n[0] && -1 == n[1] &&  0 == n[2]) ) return GenPos(0x0023);//x+y (hex and cubic)
		else if( (-1 == n[0] &&  0 == n[1] &&  1 == n[2]) || ( 1 == n[0] &&  0 == n[1] && -1 == n[2]) ) return GenPos(0x0028);//x-z
		else if( ( 1 == n[0] &&  0 == n[1] &&  1 == n[2]) || (-1 == n[0] &&  0 == n[1] && -1 == n[2]) ) return GenPos(0x002D);//x+z
		else if( ( 1 == n[0] &&  2 == n[1] &&  0 == n[2]) || (-1 == n[0] && -2 == n[1] &&  0 == n[2]) ) return GenPos(0x0034);//(hex only)
		else if( ( 2 == n[0] &&  1 == n[1] &&  0 == n[2]) || (-2 == n[0] && -1 == n[1] &&  0 == n[2]) ) return GenPos(0x003A);//(hex only)
		else throw std::runtime_error("unsupported mirror");
	}

	//@brief    : construct a two fold rotation axis about the given direction
	//@param n  : two fold axis
	//@param hex: true for 3/6 fold type groups, e.g. 2@100 == x,-y,-z for false, x-y,-y,-z for true
	//@return   : two fold rtation matrix
	GenPos GenPos::Two(int8_t const * const n, const bool hex) {
		if     ( ( 1 == n[0] &&  0 == n[1] &&  0 == n[2]) || (-1 == n[0] &&  0 == n[1] &&  0 == n[2]) ) return GenPos(hex ? 0x0035 : 0x0006);//x
		else if( ( 0 == n[0] &&  1 == n[1] &&  0 == n[2]) || ( 0 == n[0] && -1 == n[1] &&  0 == n[2]) ) return GenPos(hex ? 0x003B : 0x0005);//y
		else if( ( 0 == n[0] &&  0 == n[1] &&  1 == n[2]) || ( 0 == n[0] &&  0 == n[1] && -1 == n[2]) ) return GenPos(0x0003);//z   (hex and cubic)
		else if( ( 0 == n[0] &&  1 == n[1] && -1 == n[2]) || ( 0 == n[0] && -1 == n[1] &&  1 == n[2]) ) return GenPos(0x001F);//y-z
		else if( ( 0 == n[0] &&  1 == n[1] &&  1 == n[2]) || ( 0 == n[0] && -1 == n[1] && -1 == n[2]) ) return GenPos(0x0019);//y+z
		else if( ( 1 == n[0] && -1 == n[1] &&  0 == n[2]) || (-1 == n[0] &&  1 == n[1] &&  0 == n[2]) ) return GenPos(0x0027);//x-y (hex and cubic)
		else if( ( 1 == n[0] &&  1 == n[1] &&  0 == n[2]) || (-1 == n[0] && -1 == n[1] &&  0 == n[2]) ) return GenPos(0x0024);//x+y (hex and cubic)
		else if( (-1 == n[0] &&  0 == n[1] &&  1 == n[2]) || ( 1 == n[0] &&  0 == n[1] && -1 == n[2]) ) return GenPos(0x002F);//x-z
		else if( ( 1 == n[0] &&  0 == n[1] &&  1 == n[2]) || (-1 == n[0] &&  0 == n[1] && -1 == n[2]) ) return GenPos(0x002A);//x+z
		else if( ( 1 == n[0] &&  2 == n[1] &&  0 == n[2]) || (-1 == n[0] && -2 == n[1] &&  0 == n[2]) ) return GenPos(0x0037);//(hex only)
		else if( ( 2 == n[0] &&  1 == n[1] &&  0 == n[2]) || (-2 == n[0] && -1 == n[1] &&  0 == n[2]) ) return GenPos(0x0039);//(hex only)
		else throw std::runtime_error("unsupported two");
	}

	//@brief    : construct a three fold rotation axis about the given direction
	//@param n  : three fold axis (must by 001 or 111 type direction)
	//@param inv: true for roto inversion (-3)
	//@return   : three fold rtation matrix
	GenPos GenPos::Three(int8_t const * const n, const bool inv) {
		if(0 == n[0] && 0 == n[1] && 1 == std::abs(n[2])) return Z(inv ? -3 : 3);//handle z specially
		if     ( 1 == n[0] &&  1 == n[1] &&  1 == n[2] ) return GenPos(inv ? 0x000F : 0x0008);// 1 1 1
		else if(-1 == n[0] &&  1 == n[1] &&  1 == n[2] ) return GenPos(inv ? 0x0012 : 0x0015);//-1 1 1
		else if( 1 == n[0] && -1 == n[1] &&  1 == n[2] ) return GenPos(inv ? 0x0014 : 0x0013);// 1-1 1
		else if( 1 == n[0] &&  1 == n[1] && -1 == n[2] ) return GenPos(inv ? 0x0011 : 0x0016);// 1 1-1
		else if(-1 == n[0] && -1 == n[1] && -1 == n[2] ) return GenPos(inv ? 0x0017 : 0x0010);//-1-1-1
		else if( 1 == n[0] && -1 == n[1] && -1 == n[2] ) return GenPos(inv ? 0x000C : 0x000B);// 1-1-1
		else if(-1 == n[0] &&  1 == n[1] && -1 == n[2] ) return GenPos(inv ? 0x0009 : 0x000E);//-1 1-1
		else if(-1 == n[0] && -1 == n[1] &&  1 == n[2] ) return GenPos(inv ? 0x000A : 0x000D);//-1-1 1
		else throw std::runtime_error("unsupported three");
	}

	//@brief  : construct an N fold rotation about the Z axis
	//@param n: rotational order (must be +/- 1,2,3,4, or 6)
	//@return : rotation matrix
	//@note   : 1, -1 and 2 return idenity, inversion, and mirror respectively
	GenPos GenPos::Z(const int8_t n) {
		switch(n) {
			case  1: return Identity();
			case  2: return GenPos(0x0003);
			case  3: return GenPos(0x003C);
			case  4: return GenPos(0x0021);
			case  6: return GenPos(0x0030);
			case -1: return Inversion();
			case -2: return GenPos(0x0004);
			case -3: return GenPos(0x003F);
			case -4: return GenPos(0x0026);
			case -6: return GenPos(0x0033);
			default: throw std::runtime_error("unsupported Z");
		}
	}

	////////////////////////////////////////////////
	//       3x3 + augmenting vector access       //
	////////////////////////////////////////////////

	//@brief  : set the upper left 3x3 sub matrix
	//@param m: matrix to use (must be one of the 64 valid 3x3 position matrices)
	void GenPos::setMat3(int8_t const * m) {
		uint_fast8_t i = Mat3ToIdx(m);
		if(0xFF == i) throw std::runtime_error("matrix isn't one of 64 valid 3x3 general positions");
		set3x3(i);
	}

	//@brief  : set the translation * 24
	//@param t: translation for x,y,z each in [0,1) * 24
	void GenPos::setTrans(int8_t* t) {
		//copy translations
		u.c[1] = t[0];
		u.c[2] = t[1];
		u.c[3] = t[2];

		//bring to [0,1)
		for(size_t i = 0; i < 3; i++) {
			while(u.c[i+1] < 0) u.c[i+1] += 24;
			u.c[i+1] = u.c[i+1] % 24;
		}
	}

	//@brief  : get the translation
	//@param t: location to write translation {a, b, c}
	template <typename Real> void GenPos::getTransReal(Real * const t) const {
		std::transform(getTrans(), getTrans()+3, t, [](const int v){return Real(v) / 24;});
	}

	//@brief  : set the translation
	//@param t: translation to set {a, b, c}
	//@note   : translations must be one of (0, 1/6, 1/4, 1/3, 1/2, 2/3, 3/4, or 5/6) +/- n*1 for normal
	template <typename Real> void GenPos::setTransReal(Real const * const t) {
		int8_t iTrans[3];//integer translation
		for(size_t i = 0; i < 3; i++) {
			Real v = t[i] * 24;//multiply by 24
			int vi = (int)std::round(v);//convert to nearest fraction
			if(std::fabs(v - vi) > std::numeric_limits<Real>::epsilon() * 10) throw std::runtime_error("translation too far from fraction of 24");
			while(vi <  0) vi += 24;
			while(vi > 23) vi -= 24;
			iTrans[i] = (int8_t)vi;//save integer version - cast is safe since we brought vi to [0,24)
		}
		setTrans(iTrans);
	}

	//@brief   : transform a hexagonal general position matrix from a 3 fold coordinate system to a cartesian frame
	//@param om: location to write 3x3 matrix in cartesian frame
	//@note    : throws for cubic matricies
	template <typename Real> void GenPos::getMat3HexCart(Real * const om) {
		//compute A * m * A^-1 where A = {1, -0.5, 0}, {0, sqrt(3)/2, 0}, {0, 0, 1} ==> A^-1 = {1, 1/sqrt(3), 0}, {0, 2/sqrt(3), 0}, {0, 0, 1}
		//this is a simplification assuming m[2,5,6,7] == 0
		int8_t const * m = getMat3();
		if(!(0 == m[2] && 0 == m[5] && 0 == m[6] && 0 == m[7])) throw std::runtime_error("getMat3HexCart requires hexagonal type matrix");
		static const Real k32 = std::sqrt(Real(3))/2;
		om[0] = Real(m[0] * 2 - m[3]) / 2;
		om[1] =(Real(m[0] * 2 - m[3] + m[1] * 4 - m[4] * 2) / 3) * k32;
		om[3] = Real(m[3]) * k32;
		om[4] = Real(m[3] + m[4] * 2) / 2;
		om[2] = om[5] = om[6] = om[7] = 0;
		om[8] = m[8];
	}

	////////////////////////////////////////////////
	//              property queries              //
	////////////////////////////////////////////////

	//@brief : get the determinant of the 3x3 matrix
	//@return: determinant
	int_fast8_t GenPos::det() const {
		int8_t const* m = getMat3();
		return (m[0] * m[4] * m[8] + m[1] * m[5] * m[6] + m[2] * m[3] * m[7]) -
		       (m[0] * m[5] * m[7] + m[1] * m[3] * m[8] + m[2] * m[4] * m[6]);
	}

	//@brief : get the trace of the 3x3 matrix
	//@return: trace
	int_fast8_t GenPos::tr() const {
		int8_t const* m = getMat3();
		return m[0] + m[4] + m[8];
	}

	//@brief : get the type of the 3x3 matrix
	//@return: type
	//@note  : type is n * det() such that matrix^n == I and must be +/- 1, 2, 3, 4, 6
	//         e.g. 1 for identity, -2 for mirror plane, 3 for 3 fold rotation etc
	int_fast8_t GenPos::order() const {
		//get matrix and trace/determinant
		int8_t const* m = getMat3();
		int_fast8_t d = (m[0] * m[4] * m[8] + m[1] * m[5] * m[6] + m[2] * m[3] * m[7])
		              - (m[0] * m[5] * m[7] + m[1] * m[3] * m[8] + m[2] * m[4] * m[6]);
		int_fast8_t t = m[0] + m[4] + m[8];

		//get type from both (international tables table 11.2.1.1)
		if(1 == d) {
			switch(t) {
				case -1: return 2;
				case  0: return 3;
				case  1: return 4;
				case  2: return 6;
				case  3: return 1;
				default: throw std::logic_error("unexpected trace for det == +1");
			}
		} else if(-1 == d) {
			switch(t) {
				case -3: return -1;
				case -2: return -6;
				case -1: return -4;
				case  0: return -3;
				case  1: return -2;
				default: throw std::logic_error("unexpected trace for det == -1");
			}
		} else {
			throw std::logic_error("non +/-1 determinant");
		}
		return 0;
	}

	//@brief : get the axis associated with the 3x3 matrix (i.e. rotation axis or mirror plane normal)
	//@return: rotation axis as e.g. {2,1,0}
	//@note  : identity and inversion return {1,0,0}
	int8_t const* GenPos::axis() const {
		//get matrix and principal eigenvector
		int8_t const* m = getMat3();//get matrix we need principal eigenvector of
		const int8_t d = det();//determine principal eigenvalue: +/-1

		//build A - I * y (all elements are 0, +/-1, +/-2)
		int_fast8_t a[9] = {
			(int_fast8_t)(m[0] - d), (int_fast8_t)(m[1]    ), (int_fast8_t)(m[2]    ),
			(int_fast8_t)(m[3]    ), (int_fast8_t)(m[4] - d), (int_fast8_t)(m[5]    ),
			(int_fast8_t)(m[6]    ), (int_fast8_t)(m[7]    ), (int_fast8_t)(m[8] - d),
		};

		//there actually extremely limited options for eigen vectors
		//the easiest solution is just to check all of them
		static const int8_t v[15][3] = {
			{ 1, 0, 0}, { 0, 1, 0}, { 0, 0, 1}, { 0, 1, 1}, { 1, 0, 1},
			{ 1, 1, 0}, { 0,-1, 1}, {-1, 0, 1}, { 1,-1, 0}, { 1, 1, 1},
			{ 1,-1, 1}, {-1,-1, 1}, {-1, 1, 1}, { 2, 1, 0}, { 1, 2, 0},
		};

		//loop over possible eigen vectors searching for a match
		for(size_t i = 0; i < 15; i++) {
			if(0 != v[i][0] * a[0] + v[i][1] * a[1] + v[i][2] * a[2]) continue;//does v * the first row == 0
			if(0 != v[i][0] * a[3] + v[i][1] * a[4] + v[i][2] * a[5]) continue;//does v * the second row == 0
			if(0 != v[i][0] * a[6] + v[i][1] * a[7] + v[i][2] * a[8]) continue;//does v * the third row == 0
			return v[i];//if we made it this far (A-I*y)*v == 0
		}
		throw std::logic_error("no suitable eigenvector found");//we're in trouble if we didn't find a match
	}

	////////////////////////////////////////////////
	//             matrix operations              //
	////////////////////////////////////////////////

	//@brief    : multiply this general position with another
	//@param rhs: other general position to multiply with
	//@return   : this
	GenPos& GenPos::operator*=(const GenPos& rhs) {
		//first extract both 3x3 matrices and the translations
		int8_t const * a  =     getMat3 ();
		int8_t const * b  = rhs.getMat3 ();
		int8_t const * ta =     getTrans();
		int8_t const * tb = rhs.getTrans();

		//do the 3x3 multiplication (c = a * b)
		int8_t c[9] = { 
			(int8_t)(a[0] * b[0] + a[1] * b[3] + a[2] * b[6]), (int8_t)(a[0] * b[1] + a[1] * b[4] + a[2] * b[7]), (int8_t)(a[0] * b[2] + a[1] * b[5] + a[2] * b[8]),
			(int8_t)(a[3] * b[0] + a[4] * b[3] + a[5] * b[6]), (int8_t)(a[3] * b[1] + a[4] * b[4] + a[5] * b[7]), (int8_t)(a[3] * b[2] + a[4] * b[5] + a[5] * b[8]),
			(int8_t)(a[6] * b[0] + a[7] * b[3] + a[8] * b[6]), (int8_t)(a[6] * b[1] + a[7] * b[4] + a[8] * b[7]), (int8_t)(a[6] * b[2] + a[7] * b[5] + a[8] * b[8]),
		};

		//compute the new translation
		int8_t t[3] = {
			(int8_t)(ta[0] + tb[0] * a[0] + tb[1] * a[1] + tb[2] * a[2]),
			(int8_t)(ta[1] + tb[0] * a[3] + tb[1] * a[4] + tb[2] * a[5]),
			(int8_t)(ta[2] + tb[0] * a[6] + tb[1] * a[7] + tb[2] * a[8]),
		};

		//bring translation to [0,1)
		for(size_t i = 0; i < 3; i++) {
			while(t[i] <  0) t[i] += 24;
			while(t[i] > 23) t[i] -= 24;
		}

		//save results and return
		setMat3(c);
		setTrans(t);
		return *this;
	}

	//@brief    : get the inverse of this general position
	//@return   : this^-1
	GenPos GenPos::inverse() const {
		//lookup table for inverse of 3x3 component
		static const int8_t lut[64] = {
			0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x10, 0x14, 0x11, 0x15, 0x12, 0x16, 0x13, 0x17,
			0x08, 0x0A, 0x0C, 0x0E, 0x09, 0x0B, 0x0D, 0x0F, 0x18, 0x19, 0x1C, 0x1D, 0x1A, 0x1B, 0x1E, 0x1F,
			0x20, 0x22, 0x21, 0x23, 0x24, 0x26, 0x25, 0x27, 0x28, 0x2C, 0x2A, 0x2E, 0x29, 0x2D, 0x2B, 0x2F,
			0x3E, 0x3F, 0x3C, 0x3D, 0x34, 0x35, 0x36, 0x37, 0x38, 0x39, 0x3A, 0x3B, 0x32, 0x33, 0x30, 0x31, 
		};

		//create matrix from inverse of 3x3 part and get 3x3 part
		GenPos p((uint16_t)lut[get3x3()]);
		int8_t const* mInv = p.getMat3();

		//extra transformation and compute mInv * t
		int8_t const * t = getTrans();
		int8_t tNew[3] = {
			(int8_t) -(mInv[0] * t[0] + mInv[1] * t[1] + mInv[2] * t[2]),
			(int8_t) -(mInv[3] * t[0] + mInv[4] * t[1] + mInv[5] * t[2]),
			(int8_t) -(mInv[6] * t[0] + mInv[7] * t[1] + mInv[8] * t[2]),
		};

		//bring new translation back to [0,1) and save
		for(size_t i = 0; i < 3; i++) {
			while(tNew[i] < 0) tNew[i] += 24;
			tNew[i] = tNew[i] % 24;
		}
		p.setTrans(tNew);
		return p;
	}

	//@brief  : shift the origin by a given vector
	//@param p: origin translation for x,y,z each in [0,1) * 24
	//@note   : this applys a special case of basis transform A' = P^-1 * A * P for P == I | p
	GenPos GenPos::shiftOrigin(int8_t const*const p) const {
		//get the matrix and translation
		int8_t const * W = getMat3 ();
		int8_t const * w = getTrans();

		//now compute new translation as w' = w + (W - I) * p
		int8_t wp[3] = {
			int8_t(w[0] + (W[0] - 1) * p[0] +  W[1]      * p[1] +  W[2]      * p[2]),
			int8_t(w[1] +  W[3]      * p[0] + (W[4] - 1) * p[1] +  W[5]      * p[2]),
			int8_t(w[2] +  W[6]      * p[0] +  W[7]      * p[1] + (W[8] - 1) * p[2]),
		};

		//bring w' back to [0,24) and save
		for(size_t i = 0; i < 3; i++) {
			while(wp[i] <  0) wp[i] += 24;
			while(wp[i] > 23) wp[i] -= 24;
		};
		GenPos pos(*this);
		pos.setTrans(wp);
		return pos;
	}

	//@brief    : build closed set of general positions
	//@param gen: set of general positions to close
	//@return   : closed set of general positions
	std::vector<GenPos> GenPos::CloseSet(std::vector<GenPos> gen) {
		if(gen.empty()) return gen;//handle empty set

		//initialize closed set with identity, gen, and gen^2
		std::set<GenPos> s;
		s.insert(GenPos::Identity());
		for(const GenPos& p : gen) s.insert(p  );
		for(const GenPos& p : gen) s.insert(p*p);

		//now repeatedly multiply until set is closed
		size_t num = 0;
		while(s.size() > num) {//are there new positions since the last loop?
			num = s.size();//how many positions did we start with
			for(const GenPos& a : s) {
				for(const GenPos& b : s) {
					s.insert(a * b);//add all combination of existing general positions
				}
			}
		}
		return std::vector<GenPos>(s.cbegin(), s.cend());
	}

	////////////////////////////////////////////////
	//                     IO                     //
	////////////////////////////////////////////////

	//@brief    : convert to string representation
	//@param pre: prefix for each line of matrix
	//@return   : string representation of 4x3 general position matrix
	std::string GenPos::to_string(std::string pre) const {
		//get 3x3 matrix and translation
		int8_t const * m = getMat3 ();
		int8_t const * t = getTrans();

		//now write out string
		std::stringstream ss;
		const bool uni = true;
		for(size_t j = 0; j < 3; j++) {
			ss << pre;
			for(size_t i = 0; i < 3; i++) ss << std::setw(3) << (int)m[3*j+i];
			ss << " | " << detail::get24th(t[j], uni) << '\n';
		}
		return ss.str();
	}

	//@brief    : convert to string representation
	//@param fmt: bitmask of format options
	//@return   : string representation of 4x3 general position matrix
	std::string GenPos::to_wyckoff(const int fmt) const {
		static const bool uni = true;//should unicode be used to make nicer \bar{x}
		static const std::string fracPad   = "    ";//only 2 spaces for vulgar fractions
		static const std::string    xyz[3] = { uni ? "x" : " x", uni ? "y" : " y", uni ? "z" : " z"};//extra space without unicode to stay aligned with \bar{x} -> -x
		static const std::string negXYZ[3] = {             "-x",             "-y",             "-z"};
		static const std::string posXYZ[3] = {             "+x",             "+y",             "+z"};
		static const std::string barXYZ[3] = { uni ? "x̅" : "-x", uni ? "y̅" : "-y", uni ? "z̅" : "-z"};// overline is \u0305 so \bar{x} is "x\u0305"
		const bool hex = fmt & wyckoff::Hex;
		const bool sft[3] = {
			(fmt & wyckoff::SftX) ? true : false,
			(fmt & wyckoff::SftY) ? true : false,
			(fmt & wyckoff::SftZ) ? true : false,
		};

		//get matrix + translation
		int8_t const * const m = getMat3 ();
		int8_t const * const t = getTrans();

		std::stringstream ss;
		for(size_t j = 0; j < 3; j++) {
			//print out letters and count number printed
			size_t count = 0;
			for(size_t i = 0; i < 3; i++) {
				switch(m[j*3+i]) {
					case -1: ss << ((count++ > 0) ? negXYZ[i] : barXYZ[i]); break; 
					case  1: ss << ((count++ > 0) ? posXYZ[i] :    xyz[i]); break;
					default: break;
				}
			}

			//pad for hexagonal type matrices
			if(hex && j < 2) {//there is only 3 fold about z
				switch(count) {
					case 0: //intentional fall through
					case 3: throw std::logic_error("there should be 1 or 2 non-zero entries per matrix row in wyckoff conversion");
					case 1: ss << "  "; break;//to align with e.g. x-y
					case 2: break;
				}
			}

			//now add fraction if needed
			if(0 == t[j]) {
				if(sft[j]) ss << fracPad;
			} else {
				ss << "+" << detail::get24th(t[j], uni);
			}

			//add separator between rows
			if(j < 2) ss << ',';
		}
		return ss.str();
	}

	//@brief    : convert an EMsoft generator string to a general position matrix
	//@param str: 4 character EMsoft generator string e.g. "eBFF"
	//@return   : 16 bit encoded general position matrix
	void GenPos::fromEMsoft(char const * str) {
		//parse 3x3 generator matrix
		int8_t m[9] = {0,0,0, 0,0,0, 0,0,0};
		switch(str[0]) {
			case '1': //intentional fall through (alternate origin)
			case 'a': m[0] =  1; m[4] =  1;            m[8] =  1; break;// 0x00: identity
			case 'b': m[0] = -1; m[4] = -1;            m[8] =  1; break;// 0x03: {-1, 0, 0,  0,-1, 0,  0, 0, 1} 2@z
			case 'c': m[0] = -1; m[4] =  1;            m[8] = -1; break;// 0x05: {-1, 0, 0,  0, 1, 0,  0, 0,-1} 2@y
			case 'd': m[2] =  1; m[3] =  1;            m[7] =  1; break;// 0x08: { 0, 0, 1,  1, 0, 0,  0, 1, 0} 3@111
			case 'e': m[1] =  1; m[3] =  1;            m[8] = -1; break;// 0x24: { 0, 1, 0,  1, 0, 0,  0, 0,-1} 2@110
			case 'f': m[1] = -1; m[3] = -1;            m[8] = -1; break;// 0x27: { 0,-1, 0, -1, 0, 0,  0, 0,-1} 2@1-10
			case 'g': m[1] = -1; m[3] =  1;            m[8] =  1; break;// 0x21: { 0,-1, 0,  1, 0, 0,  0, 0, 1} 4@z
			case 'h': m[0] = -1; m[4] = -1;            m[8] = -1; break;// 0x07: {-1, 0, 0,  0,-1, 0,  0, 0,-1} inversion
			case 'i': m[0] =  1; m[4] =  1;            m[8] = -1; break;// 0x04: { 1, 0, 0,  0, 1, 0,  0, 0,-1} mz
			case 'j': m[0] =  1; m[4] = -1;            m[8] =  1; break;// 0x02: { 1, 0, 0,  0,-1, 0,  0, 0, 1} my
			case 'k': m[1] = -1; m[3] = -1;            m[8] =  1; break;// 0x23: { 0,-1, 0, -1, 0, 0,  0, 0, 1} m110
			case 'l': m[1] =  1; m[3] =  1;            m[8] =  1; break;// 0x20: { 0, 1, 0,  1, 0, 0,  0, 0, 1} m1-10
			case 'm': m[1] =  1; m[3] = -1;            m[8] = -1; break;// 0x26: { 0, 1, 0, -1, 0, 0,  0, 0,-1} -4@z
			case 'n': m[1] = -1; m[3] =  1; m[4] = -1; m[8] =  1; break;// 0x3C: { 0,-1, 0,  1,-1, 0,  0, 0, 1} 3@z

			//additions for monoclinic a
			case 'o': m[0] =  1; m[4] = -1;            m[8] = -1; break;// 0x06: { 1, 0, 0,  0,-1, 0,  0, 0,-1} 2@x
			case 'p': m[0] = -1; m[4] =  1;            m[8] =  1; break;// 0x01: {-1, 0, 0,  0, 1, 0,  0, 0, 1} mx

			default: throw std::runtime_error("failed to parse 3x3 matrix from EMsoft generator string `" + std::string(str, str + 4) + "'");
		}

		//parse the translation
		int8_t t[3];
		for(size_t i = 0; i < 3; i++) {//loop over x, y, z
			switch(str[1+i]) {
				case 'O': t[i] =  0; break;//no translation
				case 'A': t[i] =  4; break;//  1/6
				case 'B': t[i] =  6; break;//  1/4
				case 'C': t[i] =  8; break;//  1/3
				case 'D': t[i] = 12; break;//  1/2
				case 'E': t[i] = 16; break;//  2/3
				case 'F': t[i] = 18; break;//  3/4
				case 'G': t[i] = 20; break;//  5/6
				case 'X': t[i] = -9; break;// -3/8
				case 'Y': t[i] = -6; break;// -1/4
				case 'Z': t[i] = -3; break;// -1/8
				default: throw std::runtime_error("failed to parse translation from EMsoft generator string `" + std::string(str, str + 4) + "'");
			}
		}

		//save result
		setMat3(m);
		setTrans(t);
	}

	//@brief    : convert a general position matrix to an EMsoft generator string
	//@param gen: 16 bit encoded general position matrix
	//@param str: location to write 4 character EMsoft generator string e.g. "eBFF"
	void GenPos::toEMsoft(char * str, const bool ori) const {
		int8_t const * t = getTrans();

		if(ori) {
			//write 3x3 character
			if(1 != order()) throw std::runtime_error("origin shift has non-identity matrix");
			str[0] = '1';

			//write translation characters
			for(size_t i = 0; i < 3; i++) {
				switch(t[i]) {
					case  0: str[i+1] = 'O'; break;//no translation
					case 15: str[i+1] = 'X'; break;// -3/8
					case 18: str[i+1] = 'Y'; break;// -1/4
					case 21: str[i+1] = 'Z'; break;// -1/8
					default: throw std::runtime_error("couldn't map translation to EMsoft character");
				}
			}

		} else {
			//write 3x3 character
			switch(get3x3()) {
				case 0x00: str[0] = 'a'; break;
				case 0x03: str[0] = 'b'; break;
				case 0x05: str[0] = 'c'; break;
				case 0x08: str[0] = 'd'; break;
				case 0x24: str[0] = 'e'; break;
				case 0x27: str[0] = 'f'; break;
				case 0x21: str[0] = 'g'; break;
				case 0x07: str[0] = 'h'; break;
				case 0x04: str[0] = 'i'; break;
				case 0x02: str[0] = 'j'; break;
				case 0x23: str[0] = 'k'; break;
				case 0x20: str[0] = 'l'; break;
				case 0x26: str[0] = 'm'; break;
				case 0x3C: str[0] = 'n'; break;
				//added for monoclinic a
				case 0x06: str[0] = 'o'; break;
				case 0x01: str[0] = 'p'; break;
				default: throw std::runtime_error("couldn't map matrix to EMsoft character");
			}

			//write translation characters
			for(size_t i = 0; i < 3; i++) {
				switch(t[i]) {
					case  0: str[i+1] = 'O'; break;//no translation
					case  4: str[i+1] = 'A'; break;//  1/6
					case  6: str[i+1] = 'B'; break;//  1/4
					case  8: str[i+1] = 'C'; break;//  1/3
					case 12: str[i+1] = 'D'; break;//  1/2
					case 16: str[i+1] = 'E'; break;//  2/3
					case 18: str[i+1] = 'F'; break;//  3/4
					case 20: str[i+1] = 'G'; break;//  5/6
					default: throw std::runtime_error("couldn't map translation to EMsoft character");
				}
			}
		}
	}

	////////////////////////////////////////////////
	//  detail (exposed primariily for testing)   //
	////////////////////////////////////////////////

	//@brief  : get a 3x3 position matrix from an index
	//@param i: index of matrix to get, must be [0,64)
	//@return : pointer to 3x3 matrix corresponding to i (row major)
	int8_t const * GenPos::IdxToMat3(const int_fast8_t i) {
		static const int8_t lut[64][9] = {
			//m-3m general positions
			{ 1, 0, 0,  0, 1, 0,  0, 0, 1}, // 0x00 |   x  , y  , z :  1
			{-1, 0, 0,  0, 1, 0,  0, 0, 1}, // 0x01 |  -x  , y  , z :  m_x
			{ 1, 0, 0,  0,-1, 0,  0, 0, 1}, // 0x02 |   x  ,-y  , z :  m_y
			{-1, 0, 0,  0,-1, 0,  0, 0, 1}, // 0x03 |  -x  ,-y  , z :  2  _z
			{ 1, 0, 0,  0, 1, 0,  0, 0,-1}, // 0x04 |   x  , y  ,-z :  m_z
			{-1, 0, 0,  0, 1, 0,  0, 0,-1}, // 0x05 |  -x  , y  ,-z :  2  _y
			{ 1, 0, 0,  0,-1, 0,  0, 0,-1}, // 0x06 |   x  ,-y  ,-z :  2  _x
			{-1, 0, 0,  0,-1, 0,  0, 0,-1}, // 0x07 |  -x  ,-y  ,-z : -1

			{ 0, 0, 1,  1, 0, 0,  0, 1, 0}, // 0x08 |   z  , x  , y :  3^+_{ 1 1 1}
			{ 0, 0,-1,  1, 0, 0,  0, 1, 0}, // 0x09 |  -z  , x  , y : -3^+_{-1 1-1}
			{ 0, 0, 1, -1, 0, 0,  0, 1, 0}, // 0x0A |   z  ,-x  , y : -3^+_{-1-1 1}
			{ 0, 0,-1, -1, 0, 0,  0, 1, 0}, // 0x0B |  -z  ,-x  , y :  3^+_{ 1-1-1}
			{ 0, 0, 1,  1, 0, 0,  0,-1, 0}, // 0x0C |   z  , x  ,-y : -3^+_{ 1-1-1}
			{ 0, 0,-1,  1, 0, 0,  0,-1, 0}, // 0x0D |  -z  , x  ,-y :  3^+_{-1-1 1}
			{ 0, 0, 1, -1, 0, 0,  0,-1, 0}, // 0x0E |   z  ,-x  ,-y :  3^+_{-1 1-1}
			{ 0, 0,-1, -1, 0, 0,  0,-1, 0}, // 0x0F |  -z  ,-x  ,-y : -3^+_{ 1 1 1}

			{ 0, 1, 0,  0, 0, 1,  1, 0, 0}, // 0x10 |   y  , z  , x :  3^-_{ 1 1 1} 
			{ 0,-1, 0,  0, 0, 1,  1, 0, 0}, // 0x11 |  -y  , z  , x : -3^-_{-1-1 1}
			{ 0, 1, 0,  0, 0,-1,  1, 0, 0}, // 0x12 |   y  ,-z  , x : -3^-_{ 1-1-1}
			{ 0,-1, 0,  0, 0,-1,  1, 0, 0}, // 0x13 |  -y  ,-z  , x :  3^-_{-1 1-1} 
			{ 0, 1, 0,  0, 0, 1, -1, 0, 0}, // 0x14 |   y  , z  ,-x : -3^-_{-1 1-1}
			{ 0,-1, 0,  0, 0, 1, -1, 0, 0}, // 0x15 |  -y  , z  ,-x :  3^-_{ 1-1-1} 
			{ 0, 1, 0,  0, 0,-1, -1, 0, 0}, // 0x16 |   y  ,-z  ,-x :  3^-_{-1-1 1} 
			{ 0,-1, 0,  0, 0,-1, -1, 0, 0}, // 0x17 |  -y  ,-z  ,-x : -3^-_{ 1 1 1}

			{ 1, 0, 0,  0, 0, 1,  0, 1, 0}, // 0x18 |   x  , z  , y :  m  _{ 0 1-1}
			{-1, 0, 0,  0, 0, 1,  0, 1, 0}, // 0x19 |  -x  , z  , y :  2  _{ 0 1 1}
			{ 1, 0, 0,  0, 0,-1,  0, 1, 0}, // 0x1A |   x  ,-z  , y :  4^+_{ 1 0 0}
			{-1, 0, 0,  0, 0,-1,  0, 1, 0}, // 0x1B |  -x  ,-z  , y : -4^-_{ 1 0 0}
			{ 1, 0, 0,  0, 0, 1,  0,-1, 0}, // 0x1C |   x  , z  ,-y :  4^-_{ 1 0 0}
			{-1, 0, 0,  0, 0, 1,  0,-1, 0}, // 0x1D |  -x  , z  ,-y : -4^+_{ 1 0 0}
			{ 1, 0, 0,  0, 0,-1,  0,-1, 0}, // 0x1E |   x  ,-z  ,-y :  m  _{ 0 1 1}
			{-1, 0, 0,  0, 0,-1,  0,-1, 0}, // 0x1F |  -x  ,-z  ,-y :  2  _{ 0 1-1}

			{ 0, 1, 0,  1, 0, 0,  0, 0, 1}, // 0x20 |   y  , x  , z :  m  _{ 1-1 0} (also in hex)
			{ 0,-1, 0,  1, 0, 0,  0, 0, 1}, // 0x21 |  -y  , x  , z :  4^+_{ 0 0 1} 
			{ 0, 1, 0, -1, 0, 0,  0, 0, 1}, // 0x22 |   y  ,-x  , z :  4^-_{ 0 0 1} 
			{ 0,-1, 0, -1, 0, 0,  0, 0, 1}, // 0x23 |  -y  ,-x  , z :  m  _{ 1 1 0} (also in hex)
			{ 0, 1, 0,  1, 0, 0,  0, 0,-1}, // 0x24 |   y  , x  ,-z :  2  _{ 1 1 0} (also in hex)
			{ 0,-1, 0,  1, 0, 0,  0, 0,-1}, // 0x25 |  -y  , x  ,-z : -4^-_{ 0 0 1} 
			{ 0, 1, 0, -1, 0, 0,  0, 0,-1}, // 0x26 |   y  ,-x  ,-z : -4^+_{ 0 0 1} 
			{ 0,-1, 0, -1, 0, 0,  0, 0,-1}, // 0x27 |  -y  ,-x  ,-z :  2  _{ 1-1 0} (also in hex)

			{ 0, 0, 1,  0, 1, 0,  1, 0, 0}, // 0x28 |   z  , y  , x :  m  _{-1 0 1}
			{ 0, 0,-1,  0, 1, 0,  1, 0, 0}, // 0x29 |  -z  , y  , x :  4^-_{ 0 1 0}
			{ 0, 0, 1,  0,-1, 0,  1, 0, 0}, // 0x2A |   z  ,-y  , x :  2  _{ 1 0 1}
			{ 0, 0,-1,  0,-1, 0,  1, 0, 0}, // 0x2B |  -z  ,-y  , x : -4^+_{ 0 1 0}
			{ 0, 0, 1,  0, 1, 0, -1, 0, 0}, // 0x2C |   z  , y  ,-x :  4^+_{ 0 1 0}
			{ 0, 0,-1,  0, 1, 0, -1, 0, 0}, // 0x2D |  -z  , y  ,-x :  m  _{ 1 0 1}
			{ 0, 0, 1,  0,-1, 0, -1, 0, 0}, // 0x2E |   z  ,-y  ,-x : -4^-_{ 0 1 0}
			{ 0, 0,-1,  0,-1, 0, -1, 0, 0}, // 0x2F |  -z  ,-y  ,-x :  2  _{-1 0 1}

			//hexagonal general positions
			{ 1,-1, 0,  1, 0, 0,  0, 0, 1}, // 0x30 |   x-y, x  , z :  6^+_{ 0 0 1}
			{ 1,-1, 0,  1, 0, 0,  0, 0,-1}, // 0x31 |   x-y, x  ,-z : -3^-_{ 0 0 1}
			{-1, 1, 0, -1, 0, 0,  0, 0, 1}, // 0x32 |  -x+y,-x  , z :  3^-_{ 0 0 1}
			{-1, 1, 0, -1, 0, 0,  0, 0,-1}, // 0x33 |  -x+y,-x  ,-z : -6^+_{ 0 0 1}

			{ 1,-1, 0,  0,-1, 0,  0, 0, 1}, // 0x34 |   x-y,-y  , z :  m  _{ 1 2 0}
			{ 1,-1, 0,  0,-1, 0,  0, 0,-1}, // 0x35 |   x-y,-y  ,-z :  2  _{ 1 0 0}
			{-1, 1, 0,  0, 1, 0,  0, 0, 1}, // 0x36 |  -x+y, y  , z :  m  _{ 1 0 0}
			{-1, 1, 0,  0, 1, 0,  0, 0,-1}, // 0x37 |  -x+y, y  ,-z :  2  _{ 1 2 0}

			{ 1, 0, 0,  1,-1, 0,  0, 0, 1}, // 0x38 |   x  , x-y, z :  m  _{ 0 1 0}
			{ 1, 0, 0,  1,-1, 0,  0, 0,-1}, // 0x39 |   x  , x-y,-z :  2  _{ 2 1 0}
			{-1, 0, 0, -1, 1, 0,  0, 0, 1}, // 0x3A |  -x  ,-x+y, z :  m  _{ 2 1 0}
			{-1, 0, 0, -1, 1, 0,  0, 0,-1}, // 0x3B |  -x  ,-x+y,-z :  2  _{ 0 1 0}

			{ 0,-1, 0,  1,-1, 0,  0, 0, 1}, // 0x3C |  -y  , x-y, z :  3^+_{ 0 0 1}
			{ 0,-1, 0,  1,-1, 0,  0, 0,-1}, // 0x3D |  -y  , x-y,-z : -6^-_{ 0 0 1}
			{ 0, 1, 0, -1, 1, 0,  0, 0, 1}, // 0x3E |   y  ,-x+y, z :  6^-_{ 0 0 1}
			{ 0, 1, 0, -1, 1, 0,  0, 0,-1}, // 0x3F |   y  ,-x+y,-z : -3^+_{ 0 0 1}
		};
		return lut[i];
	}

	//@brief  : get the index of a 3x3 position matrix
	//@param m: 3x3 matrix to get index of
	//@return : index of matrix m (such that IdxToMat3(index) has the same matrix as m)
	//@note   : returns 0xFF if m is invalid (not in 64 possible returns of IdxToMat3)
	uint_fast8_t GenPos::Mat3ToIdx(int8_t const * const m) {
		//m-3m general position matrices are permutations of {+/-1, 0, 0}, {0, +/-1, 0}, {0, 0, +/-1}
		//6/mmm general position matrices not in m-3m all have {0,0,+/-1} in the last row
		//that means everything the absolute value of the upper left 2x2 submatrix uniquely determines the location of 1's we just need to handle +/-

		//start by determining the layout of the 1s
		uint_fast8_t idx = std::abs(m[0]) + std::abs(m[1]) * 2 + std::abs(m[3]) * 4 + std::abs(m[4]) * 8;
		switch(idx) {
			case  9: idx = 0x00; break;// { 1, 0, 0,  0, 1, 0,  0, 0, 1}
			case  4: idx = 0x08; break;// { 0, 0, 1,  1, 0, 0,  0, 1, 0}
			case  2: idx = 0x10; break;// { 0, 1, 0,  0, 0, 1,  1, 0, 0}
			case  1: idx = 0x18; break;// { 1, 0, 0,  0, 0, 1,  0, 1, 0}
			case  6: idx = 0x20; break;// { 0, 1, 0,  1, 0, 0,  0, 0, 1}
			case  8: idx = 0x28; break;// { 0, 0, 1,  0, 1, 0,  1, 0, 0}
			case  7: idx = 0x30; break;// { 1,-1, 0,  1, 0, 0,  0, 0, 1}
			case 11: idx = 0x34; break;// { 1,-1, 0,  0, 1, 0,  0, 0, 1}
			case 13: idx = 0x38; break;// { 1, 0, 0,  1,-1, 0,  0, 0, 1}
			case 14: idx = 0x3C; break;// { 0, 1, 0,  1,-1, 0,  0, 0, 1}
			default: return 0xFF;
		};

		//handle cubic and hexagonal type separately
		if(idx < 0x30) {//cubic type (8 possibilities)
			//determine if the 1 in each row is positive or negative
			//the sum of each row should be +/-1 but may not be (we'll have to check for that later)
			const int r[3] = {
				(1 - (m[0] + m[1] + m[2])) / 2,// 0/1  if (m[0] + m[1] + m[2]) == +/-1 
				(1 - (m[3] + m[4] + m[5]))            ,// 0/2 "                               "
				(1 - (m[6] + m[7] + m[8])) * 2,// 0/4 "                               "
			};
			idx += uint_fast8_t(r[0] + r[1] + r[2]);//shift by [0,8)
		} else {//hexagonal type (4 possibilities)
			const int_fast8_t xy = idx < 0x38 ? m[0]-m[1] : m[3]-m[4];//determine if we should add 0/2 to the index (is the xy row x-y or y-x)
			idx += uint_fast8_t( ( 3 - xy - m[8] ) / 2 );//also shift by 0/1 for +/-z
		}

		//finally check if the matrix actually matches
		return std::equal(m, m+9, IdxToMat3(idx)) ? idx : 0xFF;
	}
}

#endif//_position_h_
