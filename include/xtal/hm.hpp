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

#ifndef _HM_H_
#define _HM_H_

#include "position.hpp"

namespace xtal {
	//Hermann-Maguin notation
	//this class can be used to build short or full Hermann-Maguin name and generators for:
	//    * origin choice 1 (or 2 for multple setting groups)
	//    * arbitrary origin choice (user supplied)
	//    * rhombohedral or hexagonal (obverse) setting for R centered groups (reverse is trivially to implement if needed)
	//    * currently only monoclinic unique axis b cell choice 1 abc (should be able to do all 18)
	//    * currently only orhthorhombic abc (should be able to do all 6)
	//    * currently only P/I centered tetragonal (MAY be able to do multiple cell C/F)
	//    * currently only P centered trigonal/hexagonal (MAY be able to do multiple cell H)
	class HermannMaguin {
		public:
			//@brief   : construct a Hermman-Maguin notation from a space group number
			//@param sg: space group number
			HermannMaguin(const size_t sg = 1) {fromNumber(sg);}

			//@brief    : construct a Hermman-Maguin notation from a space group number
			//@param sg : space group number
			//@param alt: should the alternate international (cell choice 2) be used instead of the default (cell choice 1)
			//@note     : throws if alt == true but there isn't a second cell choice, i.e. not one of -
			//              {48,50,59,68,70,85,86,88,125,126,129,130,133,134,137,138,141,142,201,203,222,224,227,228}
			void fromNumber(const size_t sg, const bool alt = false);

			//@brief    : construct a Hermman-Maguin notation from a string
			//@param str: Hermman-Maguin representation
			void fromString(char const * str);

			//@brief : convert the stored symbol to a string
			//@return: string representation of Hermann-Maguin notation
			std::string to_string() const;

			//@brief    : convert to short symbol
			//@param str: Hermman-Maguin representation (e.g. P m m m instead of P 2/m 2/m 2/m)
			HermannMaguin shortSym() const;

			//@brief     : change axis to new setting (orthorhombic groups only)
			//@param cell: cell choice (must be one of 1, 2, or 3)
			//@param axis: new axis, must be one of: {"b", "-b", "c", "-c", "a", "-a"}
			//@return    : updated symbol with new cell
			//@note      : throws for non-monoclinic groups
			HermannMaguin changeMonoCell(const size_t cell, const std::string axis = "b") const;

			//@brief     : change axis to new setting (orthorhombic groups only)
			//@param axis: new axis, must be one of: {"abc", "bac", "cab", "cba", "bca", "acb"}
			//@return    : cell updated to new axis, e.g. Imma("cab") => Ibmm
			//@note      : throws for non-orthorhombic groups
			HermannMaguin changeOrthoAxis(const std::string axis = "abc") const;

			//@brief    : compare 2 symbols to check if they are the same
			//@param rhs: other symbol to compare against
			//@return   : this == rhs
			bool operator==(const HermannMaguin& rhs) const;

			//@brief    : compare 2 symbols to check if they are the same
			//@param rhs: other symbol to compare against
			//@return   : this != rhs
			bool operator!=(const HermannMaguin& rhs) const {return !operator==(rhs);}

			//@brief    : build generators from the Hermann-Maguin name
			//@param xyz: shift from "origin of the symbol" to apply in 24th, e.g. 3 for 1/8 (NULL to use stored origin)
			//@param hex: should the hexagonal setting be used for rhombohedral groups
			std::vector<GenPos> generators(int8_t const * const xyz = NULL, const bool rHex = true) const;

			//@brief: remove the origin shift from this symbol
			void clearOrigin() {ori[0] = ori[1] = ori[2] = 0;}

		private:

			//enumeration of symmetry types
			enum Centering : int_fast8_t {P = 0x0, C = 0xC, A = 0xA, B = 0xB, I = 0x1, F = 0xF, R = 0x3, H = 0x6};//lattice centering types
			enum Glide : int_fast8_t {m = 0x5, g = 0x7, g1 = 0x8, g2 = 0x9, a = 0xa, b = 0xb, c = 0xc, d = 0xd, e = 0xe, n = 0xf};//glide/mirror types
			enum Family : int_fast8_t {//crystal families
				Triclinic     ,
				MonoclinicB   , MonoclinicBarB, MonoclinicC   , MonoclinicBarC, MonoclinicA   , MonoclinicBarA, //ab̲c, cb̲̅a, abc̲, bac̲̅, a̲bc, a̲̅cb
				OrthohombicABC, OrthohombicBAC, OrthohombicCAB, OrthohombicCBA, OrthohombicBCA, OrthohombicACB, //abc, bac̅, cab, c̅ba, bca, ac̅b (could add 2.22)
				Tetragonal    , Trigonal      , Hexagonal     , Cubic 
			};

			//symmetry about a given axis
			struct AxisSym {
				int_fast8_t r;//rotational order (+/-)1, 2, 3, 4, or 6
				int_fast8_t s;//screw subscript (0, 1, 2, 3, 4, or 5)
				int_fast8_t g;//glide plane type
				bool operator==(const AxisSym& rhs) const {return r == rhs.r && s == rhs.s && g == rhs.g;}
			};

			//@brief : get the number of non-empty symmetry elements in this symbol
			//@return: number of elements - 1, 2, or 3 e.g. 1 for "P 1" and 3 for "F m 3 m"
			size_t numSym() const;

			//@brief: sanity check a symbol, determine the crystal family, and unabbreviate if needed
			void validate();

			//@brief     : change axis to new setting (monoclinic groups only)
			//@param axis: new axis, it and fam must be one of:
			//               Family::MonoclinicB   ,//ab̲c
			//               Family::MonoclinicBarB,//    cb̲̅a
			//               Family::MonoclinicC   ,//abc̲
			//               Family::MonoclinicBarC,//    bac̲̅
			//               Family::MonoclinicA   ,//a̲bc
			//               Family::MonoclinicBarA,//     a̲̅cb
			//@param cell: cell choice (must be one of 1, 2, or 3)
			//@return    : updated symbol with new cell
			HermannMaguin changeMonoCell(const size_t cell, const Family axis) const;

			//@brief     : change axis to new setting (orthorhombic groups only)
			//@param axis: new axis, it and fam must be one of:
			//               Family::OrthohombicABC == abc
			//               Family::OrthohombicBAC ==     bac̅
			//               Family::OrthohombicCAB == cab
			//               Family::OrthohombicCBA ==     c̅ba
			//               Family::OrthohombicBCA == bca
			//               Family::OrthohombicACB ==     ac̅b
			//@return    : updated symbol with new cell
			HermannMaguin changeOrthoAxis(const Family axis) const;

			//an extended Hermann-Maguin space group symbol is a centering, up to 3 symmetry symbols, and an origin shift
			Family      fam   ;//easiest not to compute this repeatedly
			int_fast8_t cen   ;//lattice centering
			AxisSym     sym[3];//symmetry for primary, secondary, and tertiary axis
			int_fast8_t ori[3];//shift from origin of the symbol
	};
}

#include <stdexcept>
#include <cctype>
#include <algorithm>

namespace xtal {

	//@brief    : construct a Hermman-Maguin notation from a space group number
	//@param sg : space group number
	//@param alt: should the alternate international (cell choice 2) be used instead of the default (cell choice 1)
	//@note     : throws if alt == true but there isn't a second cell choice, i.e. not one of -
	//              {48,50,59,68,70,85,86,88,125,126,129,130,133,134,137,138,141,142,201,203,222,224,227,228}
	void HermannMaguin::fromNumber(const size_t sg, const bool alt) {
		//
		// this table encodes the full Hermann-Maguin name of each space group along with 2 translations:
		//   -the translation from the 'origin of the symbol' to the international tables origin (choice 1)
		//   -the translation from international tables origin choice 1 ==> origin choice 2 (or 0, 0, 0 if there is only 1 choice)
		// the full symbol for the standard setting is encoded (i.e. axis abc for monoclinic/orthrhombic and unique axis b, cell choice 1 for monoclinic)
		//
		// the symbols are encoded 'nibble wise' (i.e. a single hex value in [0x0, 0xf])
		// it is designed to be trivial for human interpretation with a few exceptions (parenthesis list logic for non less intuitive selections)
		// the following conventions are used to encode different type of symmetry elements
		//
		//   lattice centering - possibilities are enumerated in table 1.2.1 of the international tables
		//     P = 0x0 (no centering => 0)
		//     C = 0xC
		//     A = 0xA
		//     B = 0xB
		//     I = 0x1 (1 looks similar to I)
		//     F = 0xF
		//     R = 0x3 (rhombohedral groups are 3 fold symmetric)
		//     H = 0x6 (hexagonal is 6 fold symmetric)
		//
		//   other symmetry elements - possibilities are enumerated in table 1.3.1 of the international tables
		//
		//     glide types
		//     the logic for these is the most tenuous - to avoid confusion with rotation symbols 1, 2, 3, 4, and 6 are reserved
		//     after assigning like symbols (e.g. a -> 0xa) only 5 values remain: 5, 7, 8, 9, and f
		//     {g, g1, and g2} are assigned to {7, 8, and 9} so they form a block leaving 5/f for m/n
		//     n was somewhat arbitrarily assigned to f since it is conceptually similar to d and e glide
		//     m was somewhat arbitrarily assigned to 5 since m is very common and 5 isn't used anywhere else
        //       -  =  0x0 (no glide)
        //       m  =  0x5
        //       g  =  0x7
        //       g1 =  0x8
        //       g2 =  0x9
        //       a  =  0xa
        //       b  =  0xb
        //       c  =  0xc
        //       d  =  0xd
        //       e  =  0xe
        //       n  =  0xf
        // 
        //     rotation types require 2 nibbles (1 byte) with the first and second nibble storing the rotational order and modifier flag respectively
        //       proper rotations
        //         1 = 0x10
        //         2 = 0x20
        //         3 = 0x30
        //         4 = 0x40
        //         6 = 0x60
        //       rotoinversions (0xb for \bar)
        //         -1 = 0x1b
        //         -2 = 0x2b
        //         -3 = 0x3b
        //         -4 = 0x4b
        //         -6 = 0x6b
        //       screw axis
        //         2_1                     = 0x21
        //         3_1, 3_2                = 0x31, 0x32
        //         4_1, 4_2, 4_3           = 0x41, 0x42, 0x43
        //         5_1, 5_2, 5_3, 5_4      = 0x51, 0x52, 0x53, 0x54
        //         6_1, 6_2, 6_3, 6_4, 6_5 = 0x61, 0x62, 0x63, 0x64, 0x65
        // 
		//   components of origin translations are stored as 24ths
		//   3/4 (18/24) doesn't fit in 1 nibble so it needs a special value
		//   5/8, 7/8, 2/3, and 5/6 aren't needed for any space group origin shifts
		//     0   = 0x0
		//     1/8 = 0x3 (  3/24 )
		//     1/6 = 0x4 (  4/24 )
		//     1/4 = 0x6 (  6/24 )
		//     1/3 = 0x8 (  8/24 )
		//     3/8 = 0x9 (  9/24 )
		//     1/2 = 0xb ( 12/24 )
		//     3/4 = 0x7 ( minimal conflict with other symbols [7 is only used for 'g'] and 7 == 3+4 )
        // 
        // using the above conventions a complete Hermann-Maguin symbol can be compactly stored in a human readable 64 bit integer:
		//   nibble [ 1, 1]: centering
		//   nibble [ 2, 4]: primary axis symmetry (2 nibbles for rotation followed by 1 for glide)
		//   nibble [ 5, 7]: secondary axis symmetry
		//   nibble [ 8,10]: tertiary axis symmetry
		//   nibble [11,13]: {x,y,z} translation from "origin of the symbol" to international tables origin
		//   nibble [14,16]: {x,y,z} translation from origin choice 2 to choice 1 (000 for no alternate setting)
		//
		//   for example consider space group 227 F 4/d -3 2/m (F d 3 m) [2 origin choices]
		//     the shift from 'origin of the symbol' to origin 1 is {0.375, 0.375, 0.375} --> {9,9,9} /24
		//     the shift from origin 2 to origin 1 is {0.125, 0.125, 0.125} --> {3,3,3} /24
		//      F   --> 0xf
		//      4/d --> 0x40d
		//      -3  --> 0x3b0
		//      2/m --> 0x205
		//     0x(f)(40d)(3b0)(205)(999)(333) ==> 0xf40d3b0205999333
		//
		//   similarly space for group 135 P 42/m 2/b 2/c (P 42/m b c) [1 origin choice]
		//     the shift from 'origin of the symbol' to origin 1 is {0.25, 0.75, 0} --> {6,18,0} /24
		//      P    --> 0x0
		//      42/m --> 0x425
		//      2/b  --> 0x20b
		//      2/c  --> 0x20c
		//     0x(0)(425)(20b)(20c)(670)(000) ==> 0x042520b20c670000
		// 
		static const uint64_t SGLut[230] = {
			0x0100000000000000, 0x01b0000000000000, 0x0100200100000000, 0x0100210100000000, 0xc100200100000000, // [  1,   5]
			0x0100105100000000, 0x010010c100000000, 0xc100105100000000, 0xc10010c100000000, 0x0100205100000000, // [  6,  10]
			0x0100215100060000, 0xc100205100000000, 0x010020c100006000, 0x010021c100066000, 0xc10020c100006000, // [ 11,  15]
			0x0200200200000000, 0x0200200210000000, 0x0210210200660000, 0x0210210210060000, 0xc200200210000000, // [ 16,  20]
			0xc200200200000000, 0xf200200200000000, 0x1200200200000000, 0x1210210210060000, 0x0105105200000000, // [ 21,  25]
			0x010510c210000000, 0x010c10c200000000, 0x010510a200600000, 0x010c10a210600000, 0x010f10c200060000, // [ 26,  30]
			0x010510f210000000, 0x010b10a200660000, 0x010f10a210660000, 0x010f10f200660000, 0xc105105200000000, // [ 31,  35]
			0xc10510c210000000, 0xc10c10c200000000, 0xa105105200000000, 0xa10e105200060000, 0xa10510a200600000, // [ 36,  40]
			0xa10e10a200660000, 0xf105105200000000, 0xf10d10d200990000, 0x1105105200000000, 0x110b10a200660000, // [ 41,  45]
			0x110510a200600000, 0x0205205205000000, 0x020f20f20f666777, 0x020c20c205000000, 0x020b20a20f660770, // [ 46,  50]
			0x021520520a600000, 0x020f21f20a060000, 0x020520f21a006000, 0x021c20c20a600000, 0x021b21a205660000, // [ 51,  55]
			0x021c21c20f660000, 0x020b21c215066000, 0x021f21f205660000, 0x021521520f000770, 0x021b20c21f606000, // [ 56,  60]
			0x021b21c21a666000, 0x021f21521a666000, 0xc20520c215006000, 0xc20520c21e066000, 0xc205205205000000, // [ 61,  65]
			0xc20c20c205000000, 0xc20520520e060000, 0xc20c20c20e666077, 0xf205205205000000, 0xf20d20d20d999333, // [ 66,  70]
			0x1205205205000000, 0x120b20a205660000, 0x121b21c21a666000, 0x121521521a066000, 0x0400000000000000, // [ 71,  75]
			0x0410000000000000, 0x0420000000000000, 0x0430000000000000, 0x1400000000000000, 0x1410000000670000, // [ 76,  80]
			0x04b0000000000000, 0x14b0000000000000, 0x0405000000000000, 0x0425000000000000, 0x040f0000000c0670, // [ 81,  85]
			0x042f0000000c6666, 0x1405000000000000, 0x141a000000673063, 0x0400200200000000, 0x0400210200670000, // [ 86,  90]
			0x0410200200006000, 0x0410210200673000, 0x0420200200000000, 0x0420210200676000, 0x0430200200006000, // [ 91,  95]
			0x0430210200679000, 0x1400200200000000, 0x1410200200679000, 0x0400105105000000, 0x040010b105670000, // [ 96, 100]
			0x042010c105000000, 0x042010f105660000, 0x040010c10c000000, 0x040010f10c670000, 0x042010510c000000, // [101, 105]
			0x042010b10c670000, 0x1400105105000000, 0x140010c1050c0000, 0x141010510d070000, 0x141010c10d060000, // [106, 110]
			0x04b0200105000000, 0x04b020010c006000, 0x04b0210105670000, 0x04b021010c676000, 0x04b0105200000000, // [111, 115]
			0x04b010c200006000, 0x04b010b200660000, 0x04b010f200666000, 0x14b0105200000000, 0x14b010c200006000, // [116, 120]
			0x14b0200105000000, 0x14b020010d079000, 0x0405205205000000, 0x040520c20c000000, 0x040f20b205670660, // [121, 125]
			0x040f20f20c676666, 0x040520b205670000, 0x040520f20c670000, 0x040f2052050c0670, 0x040f20c20c0c0670, // [126, 130]
			0x042520520c000000, 0x042520c205000000, 0x042f20b20c666676, 0x042f20f205666676, 0x042520b20c670000, // [131, 135]
			0x042520f205660000, 0x042f20520c0c6676, 0x042f20c2050c6676, 0x1405205205000000, 0x140520c2050c0000, // [136, 140]
			0x141a20520d073073, 0x141a20c20d063073, 0x0300000000000000, 0x0310000000000000, 0x0320000000000000, // [141, 145]
			0x3300000000000000, 0x03b0000000000000, 0x33b0000000000000, 0x0300100200000000, 0x0300200100000000, // [146, 150]
			0x0310100200004000, 0x0310200100004000, 0x0320100200008000, 0x0320200100008000, 0x3300200000000000, // [151, 155]
			0x0300105100000000, 0x0300100105000000, 0x030010c100000000, 0x030010010c000000, 0x3300105000000000, // [156, 160]
			0x330010c000000000, 0x03b0100205000000, 0x03b010020c000000, 0x03b0205100000000, 0x03b020c100000000, // [161, 165]
			0x33b0205000000000, 0x33b020c000000000, 0x0600000000000000, 0x0610000000000000, 0x0650000000000000, // [166, 170]
			0x0620000000000000, 0x0640000000000000, 0x0630000000000000, 0x06b0000000000000, 0x0605000000000000, // [171, 175]
			0x0635000000006000, 0x0600200200000000, 0x0610200200000000, 0x0650200200000000, 0x0620200200000000, // [176, 180]
			0x0640200200000000, 0x0630200200000000, 0x0600105105000000, 0x060010c10c000000, 0x063010c105000000, // [181, 185]
			0x063010510c000000, 0x06b0105200000000, 0x06b010c200000000, 0x06b0200105000000, 0x06b020010c000000, // [186, 190]
			0x0605105105000000, 0x060510c10c000000, 0x063510c105006000, 0x063510510c006000, 0x0200300000000000, // [191, 195]
			0xf200300000000000, 0x1200300000000000, 0x0210300000000000, 0x1210300000000000, 0x02053b0000000000, // [196, 200]
			0x020f3b0000666666, 0xf2053b0000000000, 0xf20d3b0000999333, 0x12053b0000000000, 0x021a3b0000666000, // [201, 205]
			0x121a3b0000666000, 0x0400300200000000, 0x0420300200000000, 0xf400300200000000, 0xf410300200666000, // [206, 210]
			0x1400300200000000, 0x0430300200000000, 0x0410300200000000, 0x1410300200000000, 0x04b0300105000000, // [211, 215]
			0xf4b0300105000000, 0x14b0300105000000, 0x04b030010f000000, 0xf4b030010c000000, 0x14b030010d000000, // [216, 220]
			0x04053b0205000000, 0x040f3b020f666666, 0x04053b020f000000, 0x040f3b0205666666, 0xf4053b0205000000, // [221, 225]
			0xf4053b020c000000, 0xf40d3b0205999333, 0xf40d3b020c333999, 0x14053b0205000000, 0x140a3b020d666000, // [226, 230]
		};

		//get the encoded name and break into 16 nibbles
		if(sg < 1 || sg > 230) throw std::runtime_error("space group number must be in [1,230]");//sanity check value
		uint64_t enc = SGLut[sg-1];
		int_fast8_t nibbles[16];
		for(size_t i = 0; i < 16; i++) {
			nibbles[15-i] = (int_fast8_t)(enc & 0x000000000000000F);//extract nibble currently in last place
			enc >>= 4;//shift encoded value over 1 nibble
		}

		//correct rotoinversions being split across 2 nibbles
		for(size_t i = 2; i < 10; i+= 3) {//loop over th 2nd nibble of the 3 symmetry elements
			if(0xb == nibbles[i]) {//this is a rotoinversion indicator
				nibbles[i-1] *= -1 ;//negate value from previous nibble (now that it is unpacked)
				nibbles[i  ]  = 0x0;//there is no screw associated with rotoinversions
			}
		}

		//correct 18/24ths being stored as 0x7 in origins (last 6 nibbles)
		for(size_t i = 10; i < 16; i++) {//loop over the translation elements
			if(0x7 == nibbles[i]) {//this is 18/24ths stored as 7
				nibbles[i] = 0x12;//7 -> 18
			}
		}

		//save family based on number and get encoded symbol from table
		if     (sg <=   2) fam = Family::Triclinic     ;
		else if(sg <=  15) fam = Family::MonoclinicB   ;
		else if(sg <=  74) fam = Family::OrthohombicABC;
		else if(sg <= 142) fam = Family::Tetragonal    ;
		else if(sg <= 167) fam = Family::Trigonal      ;
		else if(sg <= 194) fam = Family::Hexagonal     ;
		else               fam = Family::Cubic         ;

		//next save elements
		cen = nibbles[0];//extract centering from first nibble
		for(size_t i = 0; i < 3; i++) {//loop over symmetry elements
			sym[i].r = nibbles[3*i+1];//save rotation
			sym[i].s = nibbles[3*i+2];//save screw
			sym[i].g = nibbles[3*i+3];//save glide
		}

		//finally save origin translation
		std::copy(nibbles + 10, nibbles + 13, ori);//save origin choice 1
		if(alt) {//shift origin if needed
			if(0 == nibbles[13] && 0 == nibbles[14] && 0 == nibbles[15])
				throw std::runtime_error("cannot select alternate origin for space groups with only 1 choice");
			for(size_t i = 0; i < 3; i++) {
				ori[i] -= nibbles[13+i];
				if(ori[i] < 0) ori[i] += 24;
			}
		}
	}

	//@brief    : construct a Hermman-Maguin notation from a string
	//@param str: Hermman-Maguin representation
	void HermannMaguin::fromString(char const * str) {
		// wrap input name as string stream to make processing easier
		std::istringstream is(str);

		// extract first character and convert to upper case
		char c;
		if(!(is >> c)) throw std::runtime_error("failed to extract Bravais lattice centering from `" + is.str() + "'");
		c = std::toupper(c);

		//parse lattice centering from first character of symbol
		switch(c) {
			case 'P': cen = Centering::P; break;
			case 'C': cen = Centering::C; break;
			case 'A': cen = Centering::A; break;
			case 'B': cen = Centering::B; break;
			case 'I': cen = Centering::I; break;
			case 'F': cen = Centering::F; break;
			case 'R': cen = Centering::R; break;
			case 'H': cen = Centering::H; break;
			default : throw std::runtime_error("unknown lattice centering `" + std::string(1, c) + "'");
		}

		//now extract at least 1 and up to 3 additional symbols
		for(size_t i = 0; i < 3; i++) sym[i].r = sym[i].s = sym[i].g = 0;//clear all symmetry in case there aren't 3 symbols in the name
		for(size_t i = 0; i < 3; i++) {//loop over potential symmetry elements
			//get next non-whitespace character
			if(!(is >> c)) break;//if there isn't any characters left to extract we're done

			//the first character in a symbol should be a rotation or a mirror/glide type, handle roto-inversion specially
			if('-' == c) {
				sym[i].r = -1;//save negative sign in rotation axis
				if(!std::isdigit(is.peek())) std::runtime_error("`" + is.str() + "' has `-' that isn't followed by a number");
				is.get(c);//replace c with the next digit
			} else {
				sym[i].r = 1;
			}

			//now c should be either a number or a mirror/glide type
			if(!std::isalnum(c)) throw std::runtime_error("`" + is.str() + " contains an unexpected character");

			//handle numeric component (always first in compound symbols)
			if(std::isdigit(c)) {
				//extract number
				switch(c) {
					case '1'://intentional fall through
					case '2'://intentional fall through
					case '3'://intentional fall through
					case '4'://intentional fall through
					case '6': sym[i].r *= (c - '0'); break;//standard guarantees 0-9 are contigous
					default : throw std::runtime_error("`" + is.str() + " has unexpected digit");
				}

				//now check for a screw axis
				if(sym[i].r > 0) {//roto-inversion axis don't have screws
					if('0' < is.peek() && is.peek() < c) {//the next character is a number in [1, r)
						is.get(c);//replace c with the next digit
						sym[i].s = c - '0';//get screw component
					}
					//else -> not a digit or a higher number (leave for next symbol)
				}

				//lastly check for a compound numeric symbol e.g. 2/m
				if('/' == is.peek()) {
					is.get(c);//skip '/'
					if(!std::isalpha(is.peek())) throw std::runtime_error("`" + is.str() + " has `/' followed by unexpected character");//or no character
					is.get(c);//extract symbol after '/'
				} else {
					//if the next character after the numeric component is anything but '/' it is either the next symbol or a problem
				}
			}

			//there is a glide type symbol to parse (c is first character or first after '/')
			if(std::isalpha(c)) {
				c = std::tolower(c);
				if(sym[i].r < 0) throw std::runtime_error("`" + is.str() + " has a compound rotoinversion symbol (unsupported)");
				switch(c) {
					case 'm': sym[i].g = Glide::m ; break;
					case 'g': {
						switch(is.peek()) {
							case '1': sym[i].g = Glide::g1; break;
							case '2': sym[i].g = Glide::g2; break;
							default : sym[i].g = Glide::g ; break;
						}
					} break;
					case 'a': sym[i].g = Glide::a ; break;
					case 'b': sym[i].g = Glide::b ; break;
					case 'c': sym[i].g = Glide::c ; break;
					case 'd': sym[i].g = Glide::d ; break;
					case 'e': sym[i].g = Glide::e ; break;
					case 'n': sym[i].g = Glide::n ; break;
					default : throw std::runtime_error("`" + is.str() + " contains an unexpected character");
				}
			}
			//else -> the previous check for isalnum means we've already parsed a stand alone rotation
		}

		//now we have the symbol, unabbreviate and sanity check
		validate();
	}

	//@brief : convert the stored symbol to a string
	//@return: string representation of Hermann-Maguin notation
	std::string HermannMaguin::to_string() const {
		const bool abbrev = false;//should abbreiated notation be used? (not implemented yet)
		std::ostringstream os;

		//print lattice type
		switch(cen) {
			case Centering::P: os << 'P'; break;
			case Centering::C: os << 'C'; break;
			case Centering::A: os << 'A'; break;
			case Centering::B: os << 'B'; break;
			case Centering::I: os << 'I'; break;
			case Centering::F: os << 'F'; break;
			case Centering::R: os << 'R'; break;
			case Centering::H: os << 'H'; break;
		}

		//print symmetry elements
		const size_t numEl = numSym();
		for(size_t i = 0; i < numEl; i++) {
			os << ' ';
			const bool hasRot = 1 != sym[i].r;//is there a rotation?
			const bool hasGld = 0 != sym[i].g;//is there a glide?

			if(hasRot || hasGld) {
				if(hasRot) {
					os << (int)sym[i].r;
					if(0 != sym[i].s) {
						if(3 == sym[i].r && 3 == sym[i].s) {
							os << "_{1,2}";
						} else {
							os << (int)sym[i].s;
						}
					}
					if(hasGld) os << '/';
				}

				if(hasGld) {
					switch(sym[i].g) {
						case Glide::m : os << "m" ; break;
						case Glide::g : os << "g" ; break;
						case Glide::g1: os << "g1"; break;
						case Glide::g2: os << "g2"; break;
						case Glide::a : os << "a" ; break;
						case Glide::b : os << "b" ; break;
						case Glide::c : os << "c" ; break;
						case Glide::d : os << "d" ; break;
						case Glide::e : os << "e" ; break;
						case Glide::n : os << "n" ; break;
					}
				}

			} else {
				os << "1";
			}
		}
		return os.str();
	}

	//@brief    : convert to short symbol
	//@param str: Hermman-Maguin representation
	HermannMaguin HermannMaguin::shortSym() const {
		//copy full symbol
		HermannMaguin sSym(*this);

		//handle based on crystal system
		switch(fam) {
			case Family::Triclinic     :
			case Family::MonoclinicB   :
			case Family::MonoclinicBarB:
			case Family::MonoclinicC   :
			case Family::MonoclinicBarC:
			case Family::MonoclinicA   :
			case Family::MonoclinicBarA: break;//these families have no abbrerviation

			case Family::OrthohombicABC:
			case Family::OrthohombicBAC:
			case Family::OrthohombicCAB:
			case Family::OrthohombicCBA:
			case Family::OrthohombicBCA:
			case Family::OrthohombicACB:// 2/m 2/m 2/m => mmm
				if(0 != sym[0].g && 0 != sym[1].g && 0 != sym[2].g) {
					sSym.sym[0].r = 1; sSym.sym[0].s = 0;
					sSym.sym[1].r = 1; sSym.sym[1].s = 0;
					sSym.sym[2].r = 1; sSym.sym[2].s = 0;
				}
			break;

			case Family::Tetragonal   :// 4/m 2/m 2/m => 4/mmm
			case Family::Hexagonal    :// 6/m 2/m 2/m => 6/mmm
				if(0 != sym[0].g && 0 != sym[1].g && 0 != sym[2].g) {
					sSym.sym[1].r = 1; sSym.sym[1].s = 0;
					sSym.sym[2].r = 1; sSym.sym[2].s = 0;
				}
			break;

			case Family::Trigonal      :
				if(-3 == sym[0].r) {// -3, -3 2/m 1, or -3 1 2/m
					if(0 != sym[1].g) {sSym.sym[1].r = 1; sSym.sym[1].s = 0;}// -3 2/m 1 => -3 m 1
					if(0 != sym[2].g) {sSym.sym[2].r = 1; sSym.sym[2].s = 0;}// -3 1 2/m => -3 1 m
				}
			break;

			case Family::Cubic         :// 4/m -3 2/m  => m-3m
				if(0 != sym[0].g){// m-3 or m-3m
					sSym.sym[0].r = 1;// 4/m -3 => m -3 [and 4/m -3 2/m => m -3 2/m]
					sSym.sym[0].s = 0;// 4/m -3 => m -3 [and 4/m -3 2/m => m -3 2/m]
					if(0 != sym[2].g) {sSym.sym[2].r = 1; sSym.sym[2].s = 0;}// m -3 2/m => m -3 m
				}
				if(-3 == sym[1].r) sSym.sym[1].r = 3;
			break;
		}
		return sSym;
	}

	//@brief     : change axis to new setting (orthorhombic groups only)
	//@param cell: cell choice (must be one of 1, 2, or 3)
	//@param axis: new axis, must be one of: {"b", "-b", "c", "-c", "a", "-a"}
	//@return    : updated symbol with new cell
	//@note      : throws for non-monoclinic groups
	HermannMaguin HermannMaguin::changeMonoCell(const size_t cell, const std::string axis) const {
		if     (0 == axis.compare( "b")) return changeMonoCell(cell, Family::MonoclinicB   );
		else if(0 == axis.compare("-b")) return changeMonoCell(cell, Family::MonoclinicBarB);
		else if(0 == axis.compare( "c")) return changeMonoCell(cell, Family::MonoclinicC   );
		else if(0 == axis.compare("-c")) return changeMonoCell(cell, Family::MonoclinicBarC);
		else if(0 == axis.compare( "a")) return changeMonoCell(cell, Family::MonoclinicA   );
		else if(0 == axis.compare("-a")) return changeMonoCell(cell, Family::MonoclinicBarA);
		else throw std::runtime_error("new orhthorhombic cell must be one of {'b', '-b', 'c', '-c', 'a', '-a'}");
	}

	//@brief     : change axis to new setting (orthorhombic groups only)
	//@param axis: new axis, must be one of: {"abc", "bac", "cab", "cba", "bca", "acb"}
	//@return    : cell updated to new axis, e.g. Imma("cab") => Ibmm
	//@note      : throws for non-orthorhombic groups
	HermannMaguin HermannMaguin::changeOrthoAxis(const std::string axis) const {
		if     (0 == axis.compare("abc")) return changeOrthoAxis(Family::OrthohombicABC);
		else if(0 == axis.compare("bac")) return changeOrthoAxis(Family::OrthohombicBAC);
		else if(0 == axis.compare("cab")) return changeOrthoAxis(Family::OrthohombicCAB);
		else if(0 == axis.compare("cba")) return changeOrthoAxis(Family::OrthohombicCBA);
		else if(0 == axis.compare("bca")) return changeOrthoAxis(Family::OrthohombicBCA);
		else if(0 == axis.compare("acb")) return changeOrthoAxis(Family::OrthohombicACB);
		else throw std::runtime_error("new orhthorhombic cell must be one of {'abc', 'bac', 'cab', 'cba', 'bca', 'acb'}");
	}

	//@brief    : compare 2 symbols to check if they are the same
	//@param rhs: other symbol to compare against
	//@return   : this == rhs
	bool HermannMaguin::operator==(const HermannMaguin& rhs) const {
		return cen == rhs.cen &&  fam == rhs.fam && 
		       std::equal(sym, sym+3, rhs.sym)   &&
		       std::equal(ori, ori+3, rhs.ori);
	}

	//@brief    : build generators from the Hermann-Maguin name
	//@param xyz: shift from "origin of the symbol" to apply in 24th, e.g. 3 for 1/8 (NULL to use stored origin)
	//@param hex: should the hexagonal setting be used for rhombohedral groups
	std::vector<GenPos> HermannMaguin::generators(int8_t const * const xyz, const bool rHex) const {
		//get short symbol and number of symmetry elements once
		const HermannMaguin sSym = shortSym();
		const size_t numEl = sSym.numSym();

		//build up reference directions
		int8_t refDirs[3][3] = {0};
		const bool isHexType = Family::Hexagonal == fam || Family::Trigonal == fam;
		switch(fam) {
			case Family::Triclinic     : break;

			//there is still an issue with a few groups for extended orthorhombic, handle it explicitly here
			case Family::OrthohombicBAC:
			case Family::OrthohombicCAB:
			case Family::OrthohombicCBA:
			case Family::OrthohombicBCA:
			case Family::OrthohombicACB:
				//problem groups are those with only a single permutation
				if(sym[0] == sym[1] && sym[1] == sym[2]) return changeOrthoAxis(Family::OrthohombicABC).generators();
			case Family::OrthohombicABC://no problem with default axis

			case Family::MonoclinicB   : // these are all abc
			case Family::MonoclinicC   : // these are all abc
			case Family::MonoclinicA   : // these are all abc
			case Family::MonoclinicBarB://these aren't actually all the same but it doesn't matter
			case Family::MonoclinicBarC://the wrong directions should all be symmetry free
			case Family::MonoclinicBarA:
				refDirs[0][0] = 1;// 100
				refDirs[1][1] = 1;// 010
				refDirs[2][2] = 1;// 001
			break;

			case Family::Tetragonal    : //intensional fall through
			case Family::Hexagonal     : refDirs[0][2] =  1;// 0 01
			                             refDirs[1][0] =  1;// 1 00
			          refDirs[2][0] = 1; refDirs[2][1] = -1;// 1-10
			break;

			case Family::Trigonal      :
				//as hexagonal
				refDirs[0][2] = 1;                    // 0 01
				refDirs[1][0] = 1;                    // 1 00
				refDirs[2][0] = 1; refDirs[2][1] = -1;// 1-10
				if(Centering::R == cen && !rHex) {//rhombohedral setting
					refDirs[0][0] = refDirs[0][1] = 1;// 1 11
					refDirs[1][1] = -1;// 1-10
				}
			break;

			case Family::Cubic         : 
				refDirs[0][2] = 1                                ;// 001
				refDirs[1][0] = refDirs[1][1] = refDirs[1][2] = 1;// 111
				refDirs[2][0] = refDirs[2][1]                 = 1;// 110, there are actually 2 different choices for the last entry
				if(0 != sym[0].g) refDirs[2][1] = -1;// 1-10 for m3m
			break;
		}

		//next determine indicator direction
		size_t indicator = 3;//many groups have no indicator
		switch(fam) {
			//no indicator
			case Family::Triclinic     :
			case Family::MonoclinicB   :
			case Family::MonoclinicBarB:
			case Family::MonoclinicC   :
			case Family::MonoclinicBarC:
			case Family::MonoclinicA   :
			case Family::MonoclinicBarA: 
			case Family::Trigonal      : break;

			//orthorhombic groups use c axis as indicator
			case Family::OrthohombicABC:// 22(2) or mm(2) 
			case Family::OrthohombicBAC: if(0 == sym[2].g) indicator = 2; break;
			case Family::OrthohombicCAB:// (2)22 or (2)mm 
			case Family::OrthohombicCBA: if(0 == sym[0].g) indicator = 0; break;
			case Family::OrthohombicBCA:// 2(2)2 or m(2)m 
			case Family::OrthohombicACB: if(0 == sym[1].g) indicator = 1; break;

			case Family::Tetragonal    :// (4)/mmm (4)22 (4)mm (-4)2m (-4)m2 [no indicator for 4/m 4 and -4]
			case Family::Hexagonal     :// (6)/mmm (6)22 (6)mm (-6)2m (-6)m2 [no indicator for 6/m 6 and -6]
				if(3 == numEl) indicator = 0;
			break;

			case Family::Cubic         :
				if(3 == numEl) {//432, -43m or m-3m
					if(0 == sym[0].g) indicator = 0;//(4)32 (-4)3m
				} else {//23 or m-3
				}
			break;
		}

		//start by accumulating centering translations
		std::set<GenPos> gen;
		GenPos p = GenPos::Identity();
		gen.insert(p);//all space groups have the identity operation [and lattice translations]
		int8_t t[3] = {0, 0, 0};
		switch(cen) {
			case Centering::P: break;//primitive
			case Centering::C: t[0] = t[1]        = 12; break;//C-face centered
			case Centering::A:        t[1] = t[2] = 12; break;//A-face centered
			case Centering::B: t[0]        = t[2] = 12; break;//B-face centered
			case Centering::I: t[0] = t[1] = t[2] = 12; break;//Body centered
			case Centering::F://All-face centered
				t[0] =  0; t[1] = 12; t[2] = 12; p.setTrans(t); gen.insert(p);
				t[0] = 12; t[1] =  0; t[2] = 12;
			break;

			// case Centering::R: t[0] = 8; t[1] = t[2] = 16; break;
			// case Centering::R: break;//hexagonal {8, 16, 16} for 'obverse setting' {16, 8, 16} for 'reverse setting'
			case Centering::R:
				if(rHex) {
					t[0] = 8; t[1] = t[2] = 16;
				}

			 break;//hexagonal {8, 16, 16} for 'obverse setting' {16, 8, 16} for 'reverse setting'
			//D cell section 4.3.5.3
			//S centering chapter 2.1 table 2.1.2.1
			//pg 16 cell relations
			case Centering::H: throw std::runtime_error("H centering not yet supported");//{1/3, 2/3, 0}, transformation matricies in tabel 5.1.3.1
		}
		p.setTrans(t);
		gen.insert(p);
		t[0] = t[1] = t[2] = 0;

		//loop over symmetry elements building generators
		for(size_t i = 0; i < numEl; i++) {

			//extract rotation
			if(i != indicator && 1 != sSym.sym[i].r) {
				if(0 == refDirs[i][0] && 0 == refDirs[i][1] && 1 == refDirs[i][2]) {//z refDirs (can have larger rotations)
					p = GenPos::Z(sSym.sym[i].r);
				} else {//not z refDirs, max rotation is 2 for non-cubic groups
					if       (-1 == sSym.sym[i].r) {
						p = GenPos::Inversion();
					} else if( 2 == sSym.sym[i].r) {
						p = GenPos::Two(refDirs[i], isHexType);
					} else if( 3 == std::abs( sSym.sym[i].r ) ) {
						p = GenPos::Three(refDirs[i], -3 == sSym.sym[i].r);
					} else {
						throw std::runtime_error("invalid refDirs");
					}
				}

				//extract screw
				if(sSym.sym[i].s != 0) {
					size_t idx;
					if       (1 == refDirs[i][0] && 0 == refDirs[i][1] && 0 == refDirs[i][2]) {
						idx = 0;
					} else if(0 == refDirs[i][0] && 1 == refDirs[i][1] && 0 == refDirs[i][2]) {
						idx = 1;
					} else if(0 == refDirs[i][0] && 0 == refDirs[i][1] && 1 == refDirs[i][2]) {
						idx = 2;
					} else {
						throw std::runtime_error("invalid screw");
					}
					t[idx] = 24 / std::abs(sSym.sym[i].r) * sSym.sym[i].s;
					p.setTrans(t);
				}

				//adjust origin with indicator if needed
				//we don't need to worry about orthorhombic settings here since the indicator is just used to skip axis (always 2_0)
				if(1 != sSym.sym[i].r && 3 != indicator) {//there is a rotation and this group has an indicator
					//determine default application rules
					int8_t wl = sym[indicator].s * 24 / sym[indicator].r;//compute translation from indicator
					bool apply = (i+1) % 3 == indicator;//most families apply the indicator to the symbol before the indicator
					size_t idx = 2;//most families apply the indicator to the Z axis

					//handle family specific rules
					if(Family::OrthohombicABC == fam ||
					   Family::OrthohombicBAC == fam ||
					   Family::OrthohombicCAB == fam ||
					   Family::OrthohombicCBA == fam ||
					   Family::OrthohombicBCA == fam ||
					   Family::OrthohombicACB == fam) {//orthorhombic families always apply indicator to the C axis
						idx = indicator;
					} else if(Family::Cubic == fam) {//cubic families have several special rules
						apply = 2 == i;//cubic families always apply the indicator to the last symbol
						if(apply) {//(4)32 or (-4)3m 
							idx = 0;//indicator is applied normally for x translation
							t[1] += wl; if(t[1] > 23) t[1] -= 24;//indicator is also reverse applied to y translation
							t[2] += wl; if(t[2] > 23) t[2] -= 24;//indicator is also reverse applied to z translation
						} else if(2 == numEl) {//23
							throw std::logic_error("indicator for 2 symbol cubic family");
						}
					}	

					//apply indicator if needed
					if(apply) {
						t[idx] -= wl;//apply indicator
						if(t[idx] < 0) t[idx] += 24;//bring translation back to [0,1) 24ths
						p.setTrans(t);
					}
				}

				//handle 23 groups specially (2 is indicator but also needs to be included)
				//special rule: location part is (-m/n, 0, 0)
				if(2 == numEl && Family::Cubic == fam && //23 or m-3 type
				   0 == i && 0 != sym[0].s) {//2_1 axis
					t[0] -= 12;
					if(t[0] < 0) t[0] += 24;
					p.setTrans(t);
				}

				//add rotation/screw
				gen.insert(p);
				t[0] = t[1] = t[2] = 0;
			}

			//extract glide, see international tables 1.3.1
			if(0 != sym[i].g) {
				p = GenPos::Mirror(refDirs[i], isHexType);
				int8_t x = 12;//translation for d/n glide (seed with 1/2)

				switch(sym[i].g) {
					case Glide::m :	break;//no translation

					//axial glide types
					case Glide::a :
						if( (0 == refDirs[i][0] && 1 == refDirs[i][1] && 0 == refDirs[i][2]) || //y
						    (0 == refDirs[i][0] && 0 == refDirs[i][1] && 1 == refDirs[i][2]) ) {//z
							t[0] = 12;// 1/2 a
						} else throw std::runtime_error("invalid direction for a glide");
					break;

					case Glide::b :
						if( (0 == refDirs[i][0] && 0 == refDirs[i][1] && 1 == refDirs[i][2]) || //z
						    (1 == refDirs[i][0] && 0 == refDirs[i][1] && 0 == refDirs[i][2]) ) {//x
							t[1] = 12;// 1/2 b
						} else throw std::runtime_error("invalid direction for b glide");
					break;

					case Glide::c :
						if(Centering::R == cen && !rHex) {//handle rhombohedral setting specially
							if( ( 1 == refDirs[i][0] && 0 == refDirs[i][1] && 0 == refDirs[i][2]) || //  1 00
							    ( 0 == refDirs[i][0] && 1 == refDirs[i][1] && 0 == refDirs[i][2]) || //  0 10
							    (-1 == refDirs[i][0] &&-1 == refDirs[i][1] && 0 == refDirs[i][2]) || // -1-10
							    ( 1 == refDirs[i][0] &&-1 == refDirs[i][1] && 0 == refDirs[i][2]) || //  1-10
							    ( 1 == refDirs[i][0] && 2 == refDirs[i][1] && 0 == refDirs[i][2]) || //  1 20
							    (-2 == refDirs[i][0] &&-1 == refDirs[i][1] && 0 == refDirs[i][2]) ) {// -2-10
								t[0] = t[1] = t[2] = 12;// "1/2 c" == 1/2 (a+b+c) [see international table footnote for details]
							} else throw std::runtime_error("invalid direction for rhombohedral c glide");
						} else {
							if( (1 == refDirs[i][0] && 0 == refDirs[i][1] && 0 == refDirs[i][2]) || //x
							    (0 == refDirs[i][0] && 1 == refDirs[i][1] && 0 == refDirs[i][2]) || //y
							    (1 == refDirs[i][0] &&-1 == refDirs[i][1] && 0 == refDirs[i][2]) || //x-y
							    (1 == refDirs[i][0] && 1 == refDirs[i][1] && 0 == refDirs[i][2]) ) {//xy
								t[2] = 12;// 1/2 c
							} else throw std::runtime_error("invalid direction for c glide");
						}
					break;
					
					//double glide
					case Glide::e :
						if(Centering::P == cen) throw std::runtime_error("e glide requires centered cell");
						if       (0 ==          refDirs[i][0]  && 0 ==          refDirs[i][1]  && 1 ==          refDirs[i][2] ) {//z  -> (  a)/2 and b/2 
							t[0] = 12; p.setTrans(t); gen.insert(p);
							t[0] =  0; t[1] = 12;
						} else if(1 ==          refDirs[i][0]  && 0 ==          refDirs[i][1]  && 0 ==          refDirs[i][2] ) {//x  -> (  b)/2 and c/2 
							t[1] = 12; p.setTrans(t); gen.insert(p);
							t[1] =  0; t[2] = 12;
						} else if(0 ==          refDirs[i][0]  && 1 ==          refDirs[i][1]  && 0 ==          refDirs[i][2] ) {//y  -> (  c)/2 and a/2 
							t[2] = 12; p.setTrans(t); gen.insert(p);
							t[2] =  0; t[0] = 12;
						} else if(1 ==          refDirs[i][0]  && 1 == std::abs(refDirs[i][1]) && 0 ==          refDirs[i][2] ) {//xy -> (a+b)/2 and c/2 
							t[0] = t[1] = 12; p.setTrans(t); gen.insert(p);
							t[0] = t[1] =  0; t[2] = 12;
						} else if(1 ==          refDirs[i][0]  && 1 ==          refDirs[i][1]  && 1 == std::abs(refDirs[i][2])) {//yz -> (b+c)/2 and a/2 
							t[1] = t[2] = 12; p.setTrans(t); gen.insert(p);
							t[1] = t[2] =  0; t[0] = 12;
						} else if(1 == std::abs(refDirs[i][0]) && 0 ==          refDirs[i][1]  && 1 ==          refDirs[i][2] ) {//zx -> (c+a)/2 and b/2 
							t[2] = t[0] = 12; p.setTrans(t); gen.insert(p);
							t[2] = t[0] =  0; t[1] = 12;
						} else {
							throw std::runtime_error("invalid direction for e glide");
						}
					break;

					//diamond/diagonal glide
					case Glide::d :
						if(!((Family::OrthohombicABC == fam && Centering::F == cen) ||
						     (Family::OrthohombicBAC == fam && Centering::F == cen) ||
						     (Family::OrthohombicCAB == fam && Centering::F == cen) ||
						     (Family::OrthohombicCBA == fam && Centering::F == cen) ||
						     (Family::OrthohombicBCA == fam && Centering::F == cen) ||
						     (Family::OrthohombicACB == fam && Centering::F == cen) ||
						     (Family::Tetragonal     == fam && Centering::I == cen) ||
						     (Family::Cubic          == fam && Centering::P != cen)))
						    throw std::runtime_error("invalid lattice for diamond glide (must be orthorhombic F, tetragonal I, or cubic F/I");
						x = 6;// 1/4
						//intentional fall through (diamond and diagonal are the same except 1/4 vs 1/2 lattice vector)
					case Glide::n :
						if       ( 0 == refDirs[i][0] &&  0 == refDirs[i][1] &&  1 == refDirs[i][2]) {//  z -> x * (  a+-b)
							t[0] = t[1] = x;
						} else if( 1 == refDirs[i][0] &&  0 == refDirs[i][1] &&  0 == refDirs[i][2]) {//  x -> x * (  b+-c)
							t[1] = t[2] = x;
						} else if( 0 == refDirs[i][0] &&  1 == refDirs[i][1] &&  0 == refDirs[i][2]) {//  y -> x * (  c+-a)
							t[2] = t[0] = x;
						} else if( 1 == refDirs[i][0] && -1 == refDirs[i][1] &&  0 == refDirs[i][2]) {//-yx -> x * (a+b+-c)
							t[0] = t[1] = t[2] = x;
						} else if( 0 == refDirs[i][0] &&  1 == refDirs[i][1] && -1 == refDirs[i][2]) {//-zy -> x * (b+c+-a)
							t[1] = t[2] = t[0] = x;
						} else if(-1 == refDirs[i][0] &&  0 == refDirs[i][1] &&  1 == refDirs[i][2]) {//-xz -> x * (c+a+-b)
							t[2] = t[0] = t[1] = x;
						} else if( 1 == refDirs[i][0] &&  1 == refDirs[i][1] &&  0 == refDirs[i][2]) {// xy -> x * (b-a+-c) [-4 3 d type groups need +c with current origns]
							t[1] = t[2] = x; t[0] = 24 - x;
						} else if( 0 == refDirs[i][0] &&  1 == refDirs[i][1] &&  1 == refDirs[i][2]) {// yz -> x * (c-b+-a) [not used by normal space group symbols]
							t[2] = t[0] = x; t[1] = 24 - x;
						} else if( 1 == refDirs[i][0] &&  0 == refDirs[i][1] &&  1 == refDirs[i][2]) {// zx -> x * (a-c+-b) [not used by normal space group symbols]
							t[0] = t[1] = x; t[2] = 24 - x;
						} else {
							throw std::runtime_error("invalid direction for n/d glide");
						}
					break;
				}
				p.setTrans(t);
				gen.insert(p);
				t[0] = t[1] = t[2] = 0;
			}
		}

		//update origin and return
		std::vector<GenPos> vGen(gen.cbegin(), gen.cend());
		for(GenPos& p : vGen) p = p.shiftOrigin(NULL == xyz ? ori : xyz);
		return vGen;
	}

	//@brief : get the number of non-empty symmetry elements in this symbol
	//@return: number of elements - 1, 2, or 3 e.g. 1 for "P 1" and 3 for "F m 3 m"
	size_t HermannMaguin::numSym() const {
		const bool isEmpty[3] = {
			0 == sym[0].r && 0 == sym[0].s && 0 == sym[0].g,
			0 == sym[1].r && 0 == sym[1].s && 0 == sym[1].g,
			0 == sym[2].r && 0 == sym[2].s && 0 == sym[2].g,
		};
		return isEmpty[0] ? 0 : ( isEmpty[1] ? 1 : ( isEmpty[2] ? 2 : 3 ) );
	}

	//@brief: sanity check a symbol, determine the crystal family, and unabbreviate if needed
	void HermannMaguin::validate() {
		//check for forbidden individual symbols
		const size_t numEl = numSym();
		for(size_t i = 0; i < numEl; i++) {
			if(0 != sym[i].g && sym[i].r < 0) 
				throw std::runtime_error("cannot compound rotoinverion in Hermann-Maguin symbol");
			if(0 != sym[i].r == 0 || 0 != sym[i].s == 0 || 0 != sym[i].g == 0) {//is the element non-empty
				switch(sym[i].r) {
					//normal rotations have nothing special
					case  1:
					case  2:
					case  3:
					case  4:
					case  6:
					case -1:
					case -3:
					case -4:
					case -6: break;

					//handle m as -2
					case -2:
						if(0 == sym[i].g || Glide::m == sym[i].g) {
							sym[i].g = Glide::m;
							sym[i].r = 1;
						} else {
							throw std::runtime_error("-2/(not m)");
						}
					break;

					default: throw std::runtime_error("non-crystallographic rotation");
				}
			}
		}

		//now determine family and sanity check
		if(0 == numEl) throw std::runtime_error("empty Hermann-Maguin symbol");
		else if(1 == numEl) {//single symbols are easy
			switch(std::abs(int(sym[0].r))) {
				case 2 ://intentional fall through
				case 1 ://fall through for 1/m
					if(1 == std::abs(sym[0].r) && 0 == sym[0].g) {//handle 1/m
						fam = Family::Triclinic;
					} else {
						//unabbreviate with b unique axis
						fam = Family::MonoclinicB;
						sym[1] = sym[0];
						sym[2].r = 1;
						sym[0] = sym[2];
					}
					break;

				//other groups are handled trivially (cell check if after family determination)
				case 3 : fam = Family::Trigonal  ; break;
				case 4 : fam = Family::Tetragonal; break;
				case 6 : fam = Family::Hexagonal ; break;
				default: throw std::logic_error("non-crystallographic rotation");
			}
		} else {//2 or 3 symbols
			if(3 == std::abs(sym[1].r)) {//handle cubic specially
				fam = Family::Cubic;
				if(0 != sym[1].s || 0 != sym[1].g) throw std::runtime_error("cubic Hermann-Maguin symbol _3_ not have screw or glide with 3");
				fam = Family::Cubic;
				const bool m0 = 0 != sym[0].g;//point group is m__

				if(2 == numEl) {// 23 or m-3
					//uncompress m and sanity check rotation
					if(1 == sym[0].r &&  m0) sym[0].r = 2;// m3 => 2/m3
					else if(2 != sym[0].r) throw std::runtime_error("2 element Hermann-Maguin symbol _3 must be point group 23 or m-3");

					// 3 => -3 and handle space groups 205 and 206 specially
					if(m0) {
						sym[1].r = -3;// make sure we have m-3 not m3
						if(Glide::a == (Glide)sym[0].g) sym[0].s = 1;//for these 2 groups only 2_1/a is abbreviated a
					}
				} else {// 432, -43m, or m-3m
					//uncompress m's and sanity check rotations
					const bool m2 = 0 != sym[2].g;//point group is m__
					if(m0 && m2) {
						if(1 == sym[0].r) sym[0].r = 4;// m 3 m => 4/m 3 m
						if(1 == sym[2].r) sym[2].r = 2;// m 3 m => m 3 2/m
						else if(2 != sym[2].r) throw std::runtime_error("3 element Hermann-Maguin symbol m3m must be m 3 2/m");
						sym[1].r = -3;//m 3 m => m -3 m
					}

					//make sure we fall into one of 3 possible point groups
					const bool is4 = 4 == sym[0].r;
					const bool isb4 = -4 == sym[0].r;
					if       (is4  &&  m0 &&  m2) {// 4/m -3 2/m (2 ensured from previous test)
					} else if(isb4 && !m0 &&  m2 && 1 == sym[2].r) {//-43m
					} else if(is4  && !m0 && !m2 && 2 == sym[2].r) {//432
					} else throw std::runtime_error("3 element Hermann-Maguin symbol _3_ must be 432, -43m or m3m");
				}
			} else {
				//2 or 3 symbols with |2nd rotation| != 3
				for(size_t i = 1; i < numEl; i++) {
					switch(sym[i].r) {
						case 1 :
						case 2 : break;
						default: throw std::runtime_error("secondary and tertiary rotational order must be 1 or 2 for non cubic groups");
					}
				}

				//this should be monoclinic, orthorhombic, tetragonal, trigonal, or hexagonal
				switch(std::abs(int(sym[0].r))) {
					case 1 :
					case 2 : {//monoclinic or orthorhombic
						if(2 == numEl) throw std::runtime_error("monoclinic/orthorhombic symbols must have (1 or 3)/3 elements respectively");
						size_t num1 = ( (1 == sym[0].r && 0 == sym[0].g ) ? 1 : 0 )
						            + ( (1 == sym[1].r && 0 == sym[1].g ) ? 1 : 0 )
						            + ( (1 == sym[2].r && 0 == sym[2].g ) ? 1 : 0 );
			            if(3 == num1 || 1 == num1) throw std::runtime_error("unexpected Hermann-Maguin symbol (1 or 3 '1's");
			            if(2 == num1) {//2, m, or 2/m
			            	if     (1 == sym[0].r && 0 == sym[0].g )
				            	fam = Family::MonoclinicA;
			            	else if(1 == sym[1].r && 0 == sym[1].g )
				            	fam = Family::MonoclinicB;
			            	else// (1 == sym[2].r && 0 == sym[2].g )
				            	fam = Family::MonoclinicC;
			            } else {//0 == num1
							fam = Family::OrthohombicABC;//this could be relaxed in an attempt to determine the setting
							size_t numG = (0 != sym[0].g ? 1 : 0)
							            + (0 != sym[1].g ? 1 : 0)
							            + (0 != sym[2].g ? 1 : 0);
							size_t num2 = (2 == sym[0].r ? 1 : 0)
							            + (2 == sym[1].r ? 1 : 0)
							            + (2 == sym[2].r ? 1 : 0);
				            if(3 == num2 && 0 == numG) {
				            	//222
				            } else if(1 == num2 && 2 == numG) {
				            	//mm2
				            } else if(3 == num2 && 3 == numG) {
				            	//mmm (already full symbol)
				            } else if(0 == num2 && 3 == numG) {
				            	//mmm (abbreviated)
				            	sym[0].r = 2;
				            	sym[1].r = 2;
				            	sym[2].r = 2;
				            } else {
				            	throw std::runtime_error("unexpected mix of 2 fold and mirrors in orthorhombic symbol");
				            }
			            }
					} break;

					case 3 : fam = Family::Trigonal;
						// sanity check length
						if( !( (Centering::P == cen && 3 == numEl) ||
						       (Centering::R == cen && 2 == numEl) ) ) throw std::runtime_error("P/R hexagonal symbols must have 1 or 3/2 elements");
						
						if(sym[0].r < 0) {//-3m1 or -31m
							if(0 != sym[1].g) {
								if(1 == sym[1].r) sym[1].r = 2;
							} else if(0 != sym[2].g) {
								if(1 == sym[2].r) sym[2].r = 2;
							}
						} else {
							//321, 312, 3m1, or 31m
						}
					break;//3

					case 4 : fam = Family::Tetragonal;
						if(2 == numEl) throw std::runtime_error("tetragonal symbols must have 1 or 3 elements");
						if(sym[0].r < 0) {
							//-42m or -4m2
						} else {//422, 4mm, or 4/mmm
							if(0 != sym[0].g) {//4/mmm
								if(1 == sym[1].r) sym[1].r = 2;
								if(1 == sym[2].r) sym[2].r = 2;
							} else {
								//422 or 4mm
							}
						}
					break;

					case 6 : fam = Family::Hexagonal;
						if(2 == numEl) throw std::runtime_error("hexagonal symbols must have 1 or 3 elements");
						if(sym[0].r < 0) {//-6/mmm, -6m2, or -62m
							if(0 == sym[0].g) {
								//-6m2, or -62m
							} else {//-6/mmm
								if(1 == sym[1].r) sym[1].r = 2;
								if(1 == sym[2].r) sym[2].r = 2;
							}
						} else {
							//622 or 6mm
						}
					break;

					default: throw std::logic_error("non-crystallographic rotation");
				}
			}
		}

		//now that we have the family also check the centering
		if(P == cen) return;//all lattices can be primitive
		switch(fam) {
			case Family::Triclinic     : throw std::runtime_error("triclinic groups must be P centered");
			case Family::MonoclinicB   :
			case Family::MonoclinicBarB:
			case Family::MonoclinicC   :
			case Family::MonoclinicBarC:
			case Family::MonoclinicA   :
			case Family::MonoclinicBarA://this could be made more strict by splitting out cases
				switch(cen) {
					case Centering::A:
					case Centering::B:
					case Centering::C: break;
					default: throw std::runtime_error("monoclinic groups must be P, A, B, or C centered");
				}
			break;

			case Family::OrthohombicABC:
			case Family::OrthohombicBAC:
			case Family::OrthohombicCAB:
			case Family::OrthohombicCBA:
			case Family::OrthohombicBCA:
			case Family::OrthohombicACB://this could be made more strict by splitting out cases
				switch(cen) {
					case Centering::A:
					case Centering::B:
					case Centering::C:
					case Centering::I:
					case Centering::F: break;
					default: throw std::runtime_error("monoclinic groups must be P, A, B, C, I, or F centered");
				}
			break;

			case Family::Tetragonal    : if(Centering::I != cen) std::runtime_error("tetragonal groups must be P or I centered (C/F not yet implemented)"); break;
			
			case Family::Trigonal      : if(Centering::R != cen) std::runtime_error("trigonal groups must be P or R centered (H not yet implemented)"); break;
			
			case Family::Hexagonal     : throw std::runtime_error("hexagonal groups must be P centered (H not yet implemented)");
			
			case Family::Cubic         :
				switch(cen) {
					case Centering::I: break;//TODO: check consistency w/unique axis
					case Centering::F: break;//TODO: check consistency w/unique axis
					default: throw std::runtime_error("cubic groups must be P, I, or F centered");
				}
			break;
		}
	}

	//@brief     : change axis to new setting (monoclinic groups only)
	//@param axis: new axis, it and fam must be one of:
	//               Family::MonoclinicB   ,//ab̲c
	//               Family::MonoclinicBarB,//    cb̲̅a
	//               Family::MonoclinicC   ,//abc̲
	//               Family::MonoclinicBarC,//    bac̲̅
	//               Family::MonoclinicA   ,//a̲bc
	//               Family::MonoclinicBarA,//     a̲̅cb
	//@param cell: cell choice (must be one of 1, 2, or 3)
	//@return    : updated symbol with new cell
	HermannMaguin HermannMaguin::changeMonoCell(const size_t cell, const Family axis) const {
		if(!(1 == cell || 2 == cell || 3 == cell)) throw std::runtime_error("changeMonoCell new cell must be (1, 2 or 3)");
		const size_t newCell = cell - 1;//convert to 0 indexing

		//the centering permutation is the same for all non-primitive monoclinic groups
		static const int_fast8_t CenABCI[6][3] = {
			{Centering::C, Centering::A, Centering::I},//unique axis  b cell choice 1, 2, 3
			{Centering::A, Centering::C, Centering::I},//unique axis -b cell choice 1, 2, 3
			{Centering::A, Centering::B, Centering::I},//unique axis  c cell choice 1, 2, 3
			{Centering::B, Centering::A, Centering::I},//unique axis -c cell choice 1, 2, 3
			{Centering::B, Centering::C, Centering::I},//unique axis  a cell choice 1, 2, 3
			{Centering::C, Centering::B, Centering::I},//unique axis -a cell choice 1, 2, 3
		};

		//the glide permutation is almost the same for all groups
		//this is the permutation for groups 7, 13, 14 and the preferred permutation for groups 9 and 15
		//the preferred glide for groups 8 and 12 is m but this is the alternate glide for cell choices (3,1,2)
		//the alternate glide for groups 9 and 15 is this for cell choices (2,3,1)
		static const int_fast8_t GldABCN[6][3] = {
			{Glide::c, Glide::n, Glide::a},//unique axis  b cell choice 1, 2, 3
			{Glide::a, Glide::n, Glide::c},//unique axis -b cell choice 1, 2, 3
			{Glide::a, Glide::n, Glide::b},//unique axis  c cell choice 1, 2, 3
			{Glide::b, Glide::n, Glide::a},//unique axis -c cell choice 1, 2, 3
			{Glide::b, Glide::n, Glide::c},//unique axis  a cell choice 1, 2, 3
			{Glide::c, Glide::n, Glide::b},//unique axis -a cell choice 1, 2, 3
		};

		//convert new/old family to index
		static const Family FamABC[6] = {
			MonoclinicB   ,// b unique
			MonoclinicBarB,//-b unique
			MonoclinicC   ,// c unique
			MonoclinicBarC,//-c unique
			MonoclinicA   ,// a unique
			MonoclinicBarA //-a unique
		};
		const size_t uOld = std::distance(FamABC, std::find(FamABC, FamABC+6, fam ));
		const size_t uNew = std::distance(FamABC, std::find(FamABC, FamABC+6, axis));
		if(uOld > 5) throw std::runtime_error("cannot changeMonoCell of non-monoclinic group");
		if(uNew > 5) throw std::runtime_error("cannot changeMonoCell to non-monoclinic group");

		//copy symbol and get unique axis
		HermannMaguin hm(*this);//copy symbol
		hm.fam = axis;//update unique axis
		static const size_t Idx2[6] = {1, 1, 2, 2, 0, 0};
		AxisSym& sym2 = hm.sym[Idx2[uOld]];

		//next update centering/glide handling special cases (groups 8, 9, 12, and 15)
		size_t oldCell = newCell;
		if(Centering::P != cen) {//this is space group 5, 8, 9, 12, or 15
			//determine cell choice from centering and then update
			oldCell = std::distance(CenABCI[uOld], std::find(CenABCI[uOld], CenABCI[uOld]+3, cen));
			if(3 == oldCell) throw std::runtime_error("failed to determine monoclinic cell choice from centering");
			hm.cen = CenABCI[uNew][newCell];//update centering

			//now update glide
			if(0 == sym2.g || Glide::m == sym2.g) {
				//this is space group 5 (0) or default symbol for 8 (m), we only needed to update centering
			} else {
				//group 8 (alternate symbol), 9, 12, or 15
				const size_t idx = (3 + oldCell - std::distance(GldABCN[uOld], std::find(GldABCN[uOld], GldABCN[uOld]+3, sym2.g)) ) % 3;//get permutation 
				if(3 == idx) throw std::runtime_error("failed to determine monoclinic cell choice from centering + glide");//sanity check glide type
				sym2.g = GldABCN[uNew][(newCell+idx)%3];//update glide symbol

				//further sanity check symbol for groups 12 / 15
				if(2 == sym2.r) {
					if(0 == sym2.s && 0 == idx) {
						//space group 15, normal symbol
					} else if(1 == sym2.s) {
						//space group 12/15 alternate symbol
					} else {
						throw std::runtime_error("unexpected monoclinic rotation + glide combination");//sanity check glide type
					}
				}
			}
		} else if(0 != sym2.g) {//only glide needs to be (potentially) updated
			switch(sym2.g) {//check if the glide type needs to be updated
				case Glide::m: break;
				case Glide::a:
				case Glide::b:
				case Glide::c:
				case Glide::n: {
					//update glide type
					oldCell = std::distance(GldABCN[uOld], std::find(GldABCN[uOld], GldABCN[uOld]+3, sym2.g));//0 indexed => 0, 1, or 2
					if(3 == oldCell) throw std::runtime_error("failed to determine monoclinic cell choice from glide");
					sym2.g = GldABCN[uNew][newCell];
				} break;
				default: throw std::runtime_error("unexpected monoclinic glide type");
			}
		} else {
			//this symbol only has 1 cell choice
		}

		//now that we have the old and new cell type build the transformation matrix between them (see ITA 5.1.3.1)
		if(oldCell != newCell) {
			const bool cyc = ((oldCell + 1) % 3) == newCell;// true/false for 1 -> 2 -> 3 -> 1 / 1 -> 3 -> 2 -> 1
			const bool pos = uOld % 2 == 0;//true/false for +unique axis / mirrored (negative) unique axis
			size_t idx = Idx2[uOld] * 2 + (cyc == pos ? 0 : 1);
			switch(idx) {//appy the transformation
				case 0: hm.ori[1] = -(ori[1] + ori[2]); hm.ori[2] = ori[1]; break;// +/-a unique :  x  ,-y-z, y
				case 1: hm.ori[2] = -(ori[1] + ori[2]); hm.ori[1] = ori[2]; break;// +/-a unique :  x  , z  , -y-z
				case 2: hm.ori[2] = -(ori[2] + ori[0]); hm.ori[0] = ori[2]; break;// +/-b unique :  z  , y  ,-z-x
				case 3: hm.ori[0] = -(ori[2] + ori[0]); hm.ori[2] = ori[0]; break;// +/-b unique : -z-x, y  , x
				case 4: hm.ori[0] = -(ori[0] + ori[1]); hm.ori[1] = ori[0]; break;// +/-c unique : -x-y, x  , z
				case 5: hm.ori[1] = -(ori[0] + ori[1]); hm.ori[0] = ori[1]; break;// +/-c unique :  y  ,-x-y, z
			}
			const size_t idxSum = ( ( (idx+1) / 2 ) + 1 ) % 3;
			while(hm.ori[idxSum] <  0) hm.ori[idxSum] += 24;//bring the summed element back to [0,24)
		}

		//now that the symbol elements and origin have been updated, shuffle as needed to change unique axis
		if(uOld != uNew) {//some sort of change
			const bool mir = (uOld % 2) != (uNew % 2);//unique axis is mirrored
			const bool shuf = Idx2[uOld] != Idx2[uNew];//unique axis is moved

			//shuffle symbols as needed
			if(shuf) {//shuffle axis
				std::rotate(hm.sym, hm.sym + (Idx2[uOld] + 3 - Idx2[uNew]) % 3, hm.sym + 3);
				std::rotate(hm.ori, hm.ori + (Idx2[uOld] + 3 - Idx2[uNew]) % 3, hm.ori + 3);
			}

			//negate unique axis if needed
			if(mir) {//this only needs to be applied to the origin
				const size_t idx = Idx2[uNew];//new unique axis
				hm.ori[idx] = (24 - hm.ori[idx]) % 24;//negate c axis
				std::swap(hm.ori[(idx+1) % 3], hm.ori[(idx+2) % 3]);//swap non c axis
			}
		}
		return hm;
	}

	//@brief     : change axis to new setting (orthorhombic groups only)
	//@param axis: new axis, it and fam must be one of:
	//               Family::OrthohombicABC == abc
	//               Family::OrthohombicBAC ==     bac̅
	//               Family::OrthohombicCAB == cab
	//               Family::OrthohombicCBA ==     c̅ba
	//               Family::OrthohombicBCA == bca
	//               Family::OrthohombicACB ==     ac̅b
	//@return    : updated symbol with new cell
	HermannMaguin HermannMaguin::changeOrthoAxis(const Family axis) const {
		//first get old and new type
		int curC, newC;
		switch(fam) {
			case Family::OrthohombicABC: curC = 3; break;
			case Family::OrthohombicBAC: curC =-3; break;
			case Family::OrthohombicCAB: curC = 1; break;
			case Family::OrthohombicCBA: curC =-1; break;
			case Family::OrthohombicBCA: curC = 2; break;
			case Family::OrthohombicACB: curC =-2; break;
			default: throw std::runtime_error("cannot changeOrthoAxis of non-orthorhombic group");
		}
		switch(axis) {
			case Family::OrthohombicABC: newC = 3;   break;
			case Family::OrthohombicBAC: newC =-3;   break;
			case Family::OrthohombicCAB: newC = 1;   break;
			case Family::OrthohombicCBA: newC =-1;   break;
			case Family::OrthohombicBCA: newC = 2;   break;
			case Family::OrthohombicACB: newC =-2;   break;
			default: throw std::runtime_error("cannot changeOrthoAxis to non-orthorhombic group");
		}

		//next build copy of current symbol and update 
		HermannMaguin hm(*this);
		hm.fam = axis;
		const int aCurC = std::abs(curC);
		const int aNewC = std::abs(newC);
		int delta = ( aNewC + (3-aCurC) ) % 3;
		const bool mir = curC * newC < 0;

		//update centering if needed
		switch(cen) {
			case Centering::A:
			case Centering::B:
			case Centering::C: {
				int x = cen - Centering::A;//A,B,C centerin to 0,1,2
				int dx = ( delta + ( mir ? 7 - aCurC - x : x ) ) % 3;//compute new centering type as 0,1,2
				hm.cen = Centering::A + dx;;//convert back to A,B,C
			}
			default: break;
		}

		//bring glide types back to abc
		for(size_t i = 0; i < 3; i++) {//loop over symmetry elements adjusting
			switch(sym[i].g) {
				case Glide::a:
				case Glide::b:
				case Glide::c: {
					int x = sym[i].g - Glide::a;//a,b,c glide to 0,1,2
					int dx = ( delta + ( mir ? 7 - aCurC - x : x ) ) % 3;//compute new glide type as 0,1,2
					hm.sym[i].g = Glide::a + dx;//convert back to a,b,c
				}
				default: break;
			}
		}

		//update axis
		if(delta != 0) std::rotate(hm.sym, hm.sym + 3 - delta, hm.sym + 3);//c axis moves -> rotate
		if(mir) std::swap(hm.sym[aNewC % 3], hm.sym[(aNewC + 1) % 3]);//sign mismatch -> swap non c axis

		//finally update origin shift
		if(delta != 0) std::rotate(hm.ori, hm.ori + 3 - delta, hm.ori + 3);//c axis moves -> rotate
		if(mir) {
			hm.ori[aNewC-1] = (24 - hm.ori[aNewC-1]) % 24;//sign mismatch -> negate c axis
			std::swap(hm.ori[aNewC % 3], hm.ori[(aNewC + 1) % 3]);//sign mismatch -> swap non c axis
		}
		return hm;
	}
}

#endif//_HM_H_
