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

#ifndef _emgen_h_
#define _emgen_h_

#include <string>
#include <vector>

#include "xtal/position.hpp"

namespace emsoft {

	//@brief    : compress an EMsoft generator string into a 64 bit integer
	//@param gen: 40 character EMsoft generator string to compress
	//@return   : uint64_t representation of generator string
	uint64_t encode(std::string gen);

	//@brief    : decompress an EMsoft generator string from a 64 bit integer
	//@param enc: compressed generator string
	//@return   : EMsoft generator string
	std::string decode(uint64_t enc);

	//@brief    : build generator matricies from compressed EMsoft generator
	//@param enc: compressed generator string
	//@param alt: should the alternate origin be used
	//@return   : 4x4 generator matricies (empty if alt = true but there is only 1 origin)
	std::vector<xtal::GenPos> gen_from_enc(uint64_t enc, const bool alt = false);

	//@brief    : build generator matricies for a space group
	//@param sg : space group number [1,230]
	//@param alt: should an alternate setting be selected (rhombohedral instead of hex or origin choice 2 instead of 1 as appropriate)
	//@return   : 4x4 generator matricies (empty if alt = true but there is only 1 setting)
	std::vector<xtal::GenPos> gen_from_num(size_t sg, const bool alt = false);

	//@brief     : build generator matricies for an extended monoclinic space group
	//@param sg  : space group number [3,15]
	//@param cell: cell choice [1,3]
	//@param axs : unique axis {"b", "-b", "c", "-c", "a", or "-a"}
	//@return   : 4x4 generator matricies (empty if alt = true but there is only 1 setting)
	std::vector<xtal::GenPos> mono_gen_from_num(size_t sg, size_t cell, std::string axs = "b");

	//@brief     : build generator matricies for an extended orthorhombic space group
	//@param sg  : space group number [16,74]
	//@param axs : axis configuration {"abc", "bac", "cab", "cba", "bca", or "acb"}
	//@return   : 4x4 generator matricies (empty if alt = true but there is only 1 setting)
	//@warning  : unverified for space groups [51,74] and currently origin 1 only
	std::vector<xtal::GenPos> ortho_gen_from_num(size_t sg, std::string axs);

	//table of (short) space group names
	static const std::vector<std::string> SGNames = {
		"P  1      " ,"P -1      ",
		"P 2       " ,"P 21      " ,"C 2       " ,"P m       ",
		"P c       " ,"C m       " ,"C c       " ,"P 2/m     ",
		"P 21/m    " ,"C 2/m     " ,"P 2/c     " ,"P 21/c    ",
		"C 2/c     ",
		"P 2 2 2   " ,"P 2 2 21  " ,"P 21 21 2 " ,"P 21 21 21",
		"C 2 2 21  " ,"C 2 2 2   " ,"F 2 2 2   " ,"I 2 2 2   ",
		"I 21 21 21" ,"P m m 2   " ,"P m c 21  " ,"P c c 2   ",
		"P m a 2   " ,"P c a 21  " ,"P n c 2   " ,"P m n 21  ",
		"P b a 2   " ,"P n a 21  " ,"P n n 2   " ,"C m m 2   ",
		"C m c 21  " ,"C c c 2   " ,"A m m 2   " ,"A e m 2   ", // A b m 2 ==> A e m 2
		"A m a 2   " ,"A e a 2   " ,"F m m 2   " ,"F d d 2   ", // A b a 2 ==> A e a 2
		"I m m 2   " ,"I b a 2   " ,"I m a 2   " ,"P m m m   ",
		"P n n n   " ,"P c c m   " ,"P b a n   " ,"P m m a   ",
		"P n n a   " ,"P m n a   " ,"P c c a   " ,"P b a m   ",
		"P c c n   " ,"P b c m   " ,"P n n m   " ,"P m m n   ",
		"P b c n   " ,"P b c a   " ,"P n m a   " ,"C m c m   ",
		"C m c e   " ,"C m m m   " ,"C c c m   " ,"C m m e   ",// C m (c/m) a ==> C m (c/m) e
		"C c c e   " ,"F m m m   " ,"F d d d   " ,"I m m m   ",// C c c a ==> C c c e
		"I b a m   " ,"I b c a   " ,"I m m a   ",
		"P 4       " ,"P 41      " ,"P 42      " ,"P 43      ",
		"I 4       " ,"I 41      " ,"P -4      " ,"I -4      ",
		"P 4/m     " ,"P 42/m    " ,"P 4/n     " ,"P 42/n    ",
		"I 4/m     " ,"I 41/a    " ,"P 4 2 2   " ,"P 4 21 2  ",
		"P 41 2 2  " ,"P 41 21 2 " ,"P 42 2 2  " ,"P 42 21 2 ",
		"P 43 2 2  " ,"P 43 21 2 " ,"I 4 2 2   " ,"I 41 2 2  ",
		"P 4 m m   " ,"P 4 b m   " ,"P 42 c m  " ,"P 42 n m  ",
		"P 4 c c   " ,"P 4 n c   " ,"P 42 m c  " ,"P 42 b c  ",
		"I 4 m m   " ,"I 4 c m   " ,"I 41 m d  " ,"I 41 c d  ",
		"P -4 2 m  " ,"P -4 2 c  " ,"P -4 21 m " ,"P -4 21 c ",
		"P -4 m 2  " ,"P -4 c 2  " ,"P -4 b 2  " ,"P -4 n 2  ",
		"I -4 m 2  " ,"I -4 c 2  " ,"I -4 2 m  " ,"I -4 2 d  ",
		"P 4/m m m " ,"P 4/m c c " ,"P 4/n b m " ,"P 4/n n c ",
		"P 4/m b m " ,"P 4/m n c " ,"P 4/n m m " ,"P 4/n c c ",
		"P 42/m m c" ,"P 42/m c m" ,"P 42/n b c" ,"P 42/n n m",
		"P 42/m b c" ,"P 42/m n m" ,"P 42/n m c" ,"P 42/n c m",
		"I 4/m m m " ,"I 4/m c m " ,"I 41/a m d" ,"I 41/a c d",
		"P 3       " ,"P 31      " ,"P 32      " ,"R 3       ",
		"P -3      " ,"R -3      " ,"P 3 1 2   " ,"P 3 2 1   ",
		"P 31 1 2  " ,"P 31 2 1  " ,"P 32 1 2  " ,"P 32 2 1  ",
		"R 3 2     " ,"P 3 m 1   " ,"P 3 1 m   " ,"P 3 c 1   ",
		"P 3 1 c   " ,"R 3 m     " ,"R 3 c     " ,"P -3 1 m  ",
		"P -3 1 c  " ,"P -3 m 1  " ,"P -3 c 1  " ,"R -3 m    ",
		"R -3 c    ",
		"P 6       " ,"P 61      " ,"P 65      " ,"P 62      ",
		"P 64      " ,"P 63      " ,"P -6      " ,"P 6/m     ",
		"P 63/m    " ,"P 6 2 2   " ,"P 61 2 2  " ,"P 65 2 2  ",
		"P 62 2 2  " ,"P 64 2 2  " ,"P 63 2 2  " ,"P 6 m m   ",
		"P 6 c c   " ,"P 63 c m  " ,"P 63 m c  " ,"P -6 m 2  ",
		"P -6 c 2  " ,"P -6 2 m  " ,"P -6 2 c  " ,"P 6/m m m ",
		"P 6/m c c " ,"P 63/m c m" ,"P 63/m m c",
		"P 2 3     " ,"F 2 3     " ,"I 2 3     " ,"P 21 3    ",
		"I 21 3    " ,"P m 3     " ,"P n 3     " ,"F m 3     ",
		"F d 3     " ,"I m 3     " ,"P a 3     " ,"I a 3     ",
		"P 4 3 2   " ,"P 42 3 2  " ,"F 4 3 2   " ,"F 41 3 2  ",
		"I 4 3 2   " ,"P 43 3 2  " ,"P 41 3 2  " ,"I 41 3 2  ",
		"P -4 3 m  " ,"F -4 3 m  " ,"I -4 3 m  " ,"P -4 3 n  ",
		"F -4 3 c  " ,"I -4 3 d  " ,"P m 3 m   " ,"P n 3 n   ",
		"P m 3 n   " ,"P n 3 m   " ,"F m 3 m   " ,"F m 3 c   ",
		"F d 3 m   " ,"F d 3 c   " ,"I m 3 m   " ,"I a 3 d   ",
	};

}

#include <algorithm>
#include <cctype>

namespace emsoft {

	//lookup table of all elements of generator strings
	//even though there are 15 possible 3x3 matricies characters and 11 translation characters there are only 86 unique 4 character codes
	//this is a table of all of them in alphabetical order  
	static const uint8_t NumBlocks = 86+23+2;
	static const char GenLut[NumBlocks][5] = {
		"    ", "1BBB", "1BBO", "1OBB", "1OBZ", "1OYZ", "1XXX", "1YBO", 
		"1YBY", "1YYO", "1YYY", "1ZZZ", "aDDD", "aDDO", "aDOD", "aECC", 
		"aODD", "bDDD", "bDDO", "bDOD", "bDOO", "bODD", "bODO", "bOOD", 
		"bOOO", "cDDB", "cDDD", "cDDF", "cDDO", "cDOB", "cDOD", "cDOF", 
		"cODD", "cODO", "cOOD", "cOOO", "dOOO", "eBFF", "eDDD", "eFBB", 
		"eFBF", "eOOC", "eOOD", "eOOE", "eOOO", "fDDD", "fOOC", "fOOD", 
		"fOOE", "fOOO", "gDDB", "gDDD", "gDDF", "gDDO", "gODB", "gOOB", 
		"gOOD", "gOOF", "gOOO", "hBBB", "hDDD", "hDDO", "hFFF", "hODB", 
		"hODD", "iOOD", "iOOO", "jBBB", "jDDD", "jDDO", "jDOD", "jDOO", 
		"jODD", "jODO", "jOOD", "jOOO", "kOOD", "kOOO", "lBBB", "lDDD", 
		"lOOD", "lOOO", "mOOO", "nOOC", "nOOE", "nOOO",

		//additions for extended monoclinic settings
		                                                "cDOO", "iDDD",
		"iDDO", "iDOD", "iDOO", "iODD", "iODO", "oDDD", "oDDO", "oDOD",
		"oDOO", "oODD", "oODO", "oOOD", "oOOO", "pDDD", "pDDO", "pDOD",
		"pDOO", "pODD", "pODO", "pOOD", "pOOO", 

		//additions for extended orthorhombic settings
		                                        "iBBB", "pBBB",
	};

	//lookup table of encoded generator strings for all space groups (+7 extra for rhombohedral settings)
	//this is a sequence of 8 bytes with each one indexing into the GenLut table (or 0x00 for no element ==> identity)
	//each byte is actually a 7 bit integer with the leading bit reserved
	//currently the first bit on the first byte is the inversion symmetry flag so 0x0 => 0x7 for no inversion, 0x08 => 0xf yes inversion
	//this is ~4x smaller than the full string representation
	static const uint64_t SGLut[237] = {
		0x0000000000000000,0x8000000000000000,0x2300000000000000,0x2100000000000000,0x0d23000000000000,0x4b00000000000000,
		0x4a00000000000000,0x0d4b000000000000,0x0d4a000000000000,0xa300000000000000,0xa100000000000000,0x8d23000000000000,
		0xa200000000000000,0xa000000000000000,0x8d22000000000000,0x1823000000000000,0x1722000000000000,0x181c000000000000,
		0x1320000000000000,0x0d17220000000000,0x0d18230000000000,0x100e182300000000,0x0c18230000000000,0x0c13200000000000,
		0x184b000000000000,0x174a000000000000,0x184a000000000000,0x1847000000000000,0x1747000000000000,0x1848000000000000,
		0x1346000000000000,0x1845000000000000,0x1745000000000000,0x1844000000000000,0x0d184b0000000000,0x0d174a0000000000,
		0x0d184a0000000000,0x10184b0000000000,0x1018490000000000,0x1018470000000000,0x1018450000000000,0x100e184b00000000,
		0x100e184300000000,0x0c184b0000000000,0x0c18450000000000,0x0c18470000000000,0x9823000000000000,0x18233c0100000000,
		0x9822000000000000,0x18233d0200000000,0x9423000000000000,0x941a000000000000,0x931e000000000000,0x9422000000000000,
		0x981c000000000000,0x9220000000000000,0x9720000000000000,0x981a000000000000,0x181c3d0200000000,0x9122000000000000,
		0x9320000000000000,0x9321000000000000,0x8d17220000000000,0x8d15200000000000,0x8d18230000000000,0x8d18220000000000,
		0x8d16210000000000,0x0d12234003000000,0x900e182300000000,0x100e18233b0b0000,0x8c18230000000000,0x8c181c0000000000,
		0x8c13200000000000,0x8c16210000000000,0x183a000000000000,0x1737000000000000,0x1838000000000000,0x1739000000000000,
		0x0c183a0000000000,0x0c11360000000000,0x1852000000000000,0x0c18520000000000,0x983a000000000000,0x9838000000000000,
		0x18353d0700000000,0x18333c0a00000000,0x8c183a0000000000,0x0c11363f05000000,0x183a230000000000,0x18351c0000000000,
		0x1737230000000000,0x1732190000000000,0x1838230000000000,0x18331a0000000000,0x1739230000000000,0x17341b0000000000,
		0x0c183a2300000000,0x0c11361f00000000,0x183a4b0000000000,0x183a450000000000,0x18384a0000000000,0x1833440000000000,
		0x183a4a0000000000,0x183a440000000000,0x18384b0000000000,0x1838450000000000,0x0c183a4b00000000,0x0c183a4a00000000,
		0x0c11364b00000000,0x0c11364a00000000,0x1852230000000000,0x1852220000000000,0x18521c0000000000,0x18521a0000000000,
		0x18524b0000000000,0x18524a0000000000,0x1852450000000000,0x1852440000000000,0x0c18524b00000000,0x0c18524a00000000,
		0x0c18522300000000,0x0c18521f00000000,0x983a230000000000,0x983a220000000000,0x183a233d09000000,0x183a233c0a000000,
		0x983a1c0000000000,0x983a1a0000000000,0x18351c3d07000000,0x18351a3d07000000,0x9838230000000000,0x9838220000000000,
		0x1833223c08000000,0x1833233c08000000,0x98381c0000000000,0x98331a0000000000,0x18331a3c08000000,0x18331c3c08000000,
		0x8c183a2300000000,0x8c183a2200000000,0x0c11361f3f040000,0x0c11361d3f040000,0x5500000000000000,0x5300000000000000,
		0x5400000000000000,0x0f55000000000000,0xd500000000000000,0x8f55000000000000,0x5531000000000000,0x552c000000000000,
		0x5330000000000000,0x532c000000000000,0x542e000000000000,0x542c000000000000,0x0f552c0000000000,0x554d000000000000,
		0x5551000000000000,0x554c000000000000,0x5550000000000000,0x0f554d0000000000,0x0f554c0000000000,0xd531000000000000,
		0xd52f000000000000,0xd52c000000000000,0xd52a000000000000,0x8f552c0000000000,0x8f552a0000000000,0x5518000000000000,
		0x5317000000000000,0x5417000000000000,0x5418000000000000,0x5318000000000000,0x5517000000000000,0x5542000000000000,
		0xd518000000000000,0xd517000000000000,0x55182c0000000000,0x5317290000000000,0x54172b0000000000,0x54182b0000000000,
		0x5318290000000000,0x55172c0000000000,0x55184d0000000000,0x55184c0000000000,0x55174c0000000000,0x55174d0000000000,
		0x55424d0000000000,0x55414c0000000000,0x55422c0000000000,0x55412c0000000000,0xd5182c0000000000,0xd5182a0000000000,
		0xd5172a0000000000,0xd5172c0000000000,0x1823240000000000,0x100e182324000000,0x0c18232400000000,0x1320240000000000,
		0x0c13202400000000,0x9823240000000000,0x1823243c0a000000,0x900e182324000000,0x100e1823243b0b00,0x8c18232400000000,
		0x9320240000000000,0x8c13202400000000,0x1823242c00000000,0x1823242600000000,0x100e1823242c0000,0x100e151c24280000,
		0x0c1823242c000000,0x1320242500000000,0x1320242700000000,0x0c13202427000000,0x1823245100000000,0x100e182324510000,
		0x0c18232451000000,0x1823244f00000000,0x100e1823244f0000,0x0c1320244e000000,0x9823242c00000000,0x1823242c3c0a0000,
		0x9823242600000000,0x182324263c0a0000,0x900e1823242c0000,0x900e182324260000,0x100e151c24283b0b,0x100e151c24283e06,
		0x8c1823242c000000,0x8c13202427000000,0x2400000000000000,0xa400000000000000,0x2431000000000000,0x2451000000000000,
		0x244f000000000000,0xa431000000000000,0xa42d000000000000
	};

	//extended settings for monoclinic space groups
	//these have been verified against the bilbao crystallography server: http://www.cryst.ehu.es/
	static const uint32_t MonoLut[13][6][3] = {
		{  //space group 3
			// cell 1  ,   cell 2  ,   cell 3
			{0x23000000, 0x23000000, 0x23000000},{0x23000000, 0x23000000, 0x23000000},//unique axis +/-b
			{0x18000000, 0x18000000, 0x18000000},{0x18000000, 0x18000000, 0x18000000},//unique axis +/-c
			{0x64000000, 0x64000000, 0x64000000},{0x64000000, 0x64000000, 0x64000000},//unique axis +/-a
		},{//space group 4
			{0x21000000, 0x21000000, 0x21000000},{0x21000000, 0x21000000, 0x21000000},
			{0x17000000, 0x17000000, 0x17000000},{0x17000000, 0x17000000, 0x17000000},
			{0x60000000, 0x60000000, 0x60000000},{0x60000000, 0x60000000, 0x60000000},
		},{//space group 5
			{0x230D0000, 0x23100000, 0x230C0000},{0x23100000, 0x230D0000, 0x230C0000},
			{0x18100000, 0x180E0000, 0x180C0000},{0x180E0000, 0x18100000, 0x180C0000},
			{0x640E0000, 0x640D0000, 0x640C0000},{0x640D0000, 0x640E0000, 0x640C0000},
		},{//space group 6
			{0x4B000000, 0x4B000000, 0x4B000000},{0x4B000000, 0x4B000000, 0x4B000000},
			{0x42000000, 0x42000000, 0x42000000},{0x42000000, 0x42000000, 0x42000000},
			{0x6C000000, 0x6C000000, 0x6C000000},{0x6C000000, 0x6C000000, 0x6C000000},
		},{//space group 7
			{0x4A000000, 0x46000000, 0x47000000},{0x47000000, 0x46000000, 0x4A000000},
			{0x5A000000, 0x58000000, 0x5C000000},{0x5C000000, 0x58000000, 0x5A000000},
			{0x6A000000, 0x69000000, 0x6B000000},{0x6B000000, 0x69000000, 0x6A000000},
		},{//space group 8
			{0x4B0D0000, 0x4B100000, 0x4B0C0000},{0x4B100000, 0x4B0D0000, 0x4B0C0000},
			{0x42100000, 0x420E0000, 0x420C0000},{0x420E0000, 0x42100000, 0x420C0000},
			{0x6C0E0000, 0x6C0D0000, 0x6C0C0000},{0x6C0D0000, 0x6C0E0000, 0x6C0C0000},
		},{//space group 9
			{0x4A0D0000, 0x10460000, 0x470C0000},{0x10470000, 0x460D0000, 0x4A0C0000},
			{0x105A0000, 0x0E580000, 0x5C0C0000},{0x5C0E0000, 0x10580000, 0x5A0C0000},
			{0x6A0E0000, 0x690D0000, 0x6B0C0000},{0x6B0D0000, 0x690E0000, 0x6A0C0000},
		},{//space group 10
			{0xCB230000, 0xCB230000, 0xCB230000},{0xCB230000, 0xCB230000, 0xCB230000},
			{0x98420000, 0x98420000, 0x98420000},{0x98420000, 0x98420000, 0x98420000},
			{0xEC640000, 0xEC640000, 0xEC640000},{0xEC640000, 0xEC640000, 0xEC640000},
		},{//space group 11
			{0xC9210000, 0xC9210000, 0xC9210000},{0xC9210000, 0xC9210000, 0xC9210000},
			{0xC1170000, 0xC1170000, 0xC1170000},{0xC1170000, 0xC1170000, 0xC1170000},
			{0xE8600000, 0xE8600000, 0xE8600000},{0xE8600000, 0xE8600000, 0xE8600000},
		},{//space group 12
			{0xCB230D00, 0xCB231000, 0xCB230C00},{0xCB231000, 0xCB230D00, 0xCB230C00},
			{0x98421000, 0x98420E00, 0x98420C00},{0x98420E00, 0x98421000, 0x98420C00},
			{0xEC640E00, 0xEC640D00, 0xEC640C00},{0xEC640D00, 0xEC640E00, 0xEC640C00},
		},{//space group 13
			{0xA24A0000, 0x9E460000, 0xD6470000},{0xD6470000, 0x9E460000, 0xA24A0000},
			{0x945A0000, 0x92580000, 0x965C0000},{0x965C0000, 0x92580000, 0x945A0000},
			{0xE26A0000, 0xE1690000, 0xE36B0000},{0xE36B0000, 0xE1690000, 0xE26A0000},
		},{//space group 14
			{0xC8200000, 0x9A440000, 0x9C450000},{0x9C450000, 0x9A440000, 0xC8200000},
			{0x93590000, 0x91570000, 0x955B0000},{0x955B0000, 0x91570000, 0x93590000},
			{0xE65E0000, 0xE55D0000, 0xE75F0000},{0xE75F0000, 0xE55D0000, 0xE65E0000},
		},{//space group 15
			{0xA24A0D00, 0x9E104600, 0xD6470C00},{0xD6104700, 0x9E460D00, 0xA24A0C00},
			{0x94105A00, 0x920E5800, 0x965C0C00},{0x965C0E00, 0x92105800, 0x945A0C00},
			{0xE26A0E00, 0xE1690D00, 0xE36B0C00},{0xE36B0D00, 0xE1690E00, 0xE26A0C00},
		}
	};

	//extended settings for orthorhombic space groups (currently only origin choice 1)
	//these have been PARTIALLY verified against the bilbao crystallography server: http://www.cryst.ehu.es/
	static const uint64_t OrthoLut[59][6] = {
		// abc         bac        cab          cba         bca         acb
	//222 type groups
		{0x2364000000000000,0x2364000000000000,0x2364000000000000,0x2364000000000000,0x2364000000000000,0x2364000000000000},//16
		{0x6422000000000000,0x6422000000000000,0x2314000000000000,0x2314000000000000,0x1862000000000000,0x1862000000000000},//17
		{0x5e1c000000000000,0x5e1c000000000000,0x2015000000000000,0x2015000000000000,0x5f13000000000000,0x5f13000000000000},//18
		{0x5e20000000000000,0x5e20000000000000,0x5e20000000000000,0x5e20000000000000,0x5e20000000000000,0x5e20000000000000},//19
		{0x640d220000000000,0x640d220000000000,0x2314100000000000,0x2314100000000000,0x18620e0000000000,0x18620e0000000000},//20
		{0x23640d0000000000,0x23640d0000000000,0x23640d0000000000,0x23640d0000000000,0x23640d0000000000,0x23640d0000000000},//21
		{0x23640e1000000000,0x23640e1000000000,0x23640e1000000000,0x23640e1000000000,0x23640e1000000000,0x23640e1000000000},//22
		{0x23640c0000000000,0x23640c0000000000,0x23640c0000000000,0x23640c0000000000,0x23640c0000000000,0x23640c0000000000},//23
		{0x5e200c0000000000,0x5e200c0000000000,0x5e200c0000000000,0x5e200c0000000000,0x5e200c0000000000,0x5e200c0000000000},//24
	//mm2 type groups
		{0x6c4b000000000000,0x6c4b000000000000,0x4b42000000000000,0x4b42000000000000,0x6c42000000000000,0x6c42000000000000},//25
		{0x6c4a000000000000,0x4b6b000000000000,0x4b5a000000000000,0x4247000000000000,0x426a000000000000,0x6c5c000000000000},//26
		{0x6b4a000000000000,0x6b4a000000000000,0x475a000000000000,0x475a000000000000,0x6a5c000000000000,0x6a5c000000000000},//27
		{0x6847000000000000,0x496a000000000000,0x495c000000000000,0x414a000000000000,0x416b000000000000,0x685a000000000000},//28
		{0x4767000000000000,0x6a48000000000000,0x455c000000000000,0x594a000000000000,0x5b6b000000000000,0x5a66000000000000},//29
		{0x4869000000000000,0x6746000000000000,0x5946000000000000,0x4558000000000000,0x6658000000000000,0x5b69000000000000},//30
		{0x6c46000000000000,0x4b69000000000000,0x4b58000000000000,0x4246000000000000,0x4269000000000000,0x6c58000000000000},//31
		{0x4566000000000000,0x4566000000000000,0x5b48000000000000,0x5b48000000000000,0x5967000000000000,0x5967000000000000},//32
		{0x4565000000000000,0x6644000000000000,0x5b44000000000000,0x5748000000000000,0x5767000000000000,0x5965000000000000},//33
		{0x4465000000000000,0x4465000000000000,0x5744000000000000,0x5744000000000000,0x5765000000000000,0x5765000000000000},//34
		{0x6c4b0d0000000000,0x6c4b0d0000000000,0x4b42100000000000,0x4b42100000000000,0x6c420e0000000000,0x6c420e0000000000},//35
		{0x6c0d4a0000000000,0x4b0d6b0000000000,0x4b5a100000000000,0x4247100000000000,0x426a0e0000000000,0x6c5c0e0000000000},//36
		{0x0d6b4a0000000000,0x0d6b4a0000000000,0x475a100000000000,0x475a100000000000,0x6a5c0e0000000000,0x6a5c0e0000000000},//37
		{0x6c4b100000000000,0x6c4b0e0000000000,0x4b420e0000000000,0x4b420d0000000000,0x6c420d0000000000,0x6c42100000000000},//38
		{0x496a6b1000000000,0x68474a0e00000000,0x41474a0e00000000,0x495a5c0d00000000,0x685a5c0d00000000,0x416a6b1000000000},//39
		{0x6847100000000000,0x496a0e0000000000,0x495c0e0000000000,0x410d4a0000000000,0x410d6b0000000000,0x685a100000000000},//40
		{0x4566671000000000,0x4566480e00000000,0x455b480e00000000,0x595b0d4800000000,0x595b0d6700000000,0x5966671000000000},//41
		{0x6c4b0e1000000000,0x6c4b0e1000000000,0x4b420e1000000000,0x4b420e1000000000,0x6c420e1000000000,0x6c420e1000000000},//42
		{0x436e0e1000000000,0x436e0e1000000000,0x6d430e1000000000,0x6d430e1000000000,0x6d6e0e1000000000,0x6d6e0e1000000000},//43
		{0x6c4b0c0000000000,0x6c4b0c0000000000,0x4b420c0000000000,0x4b420c0000000000,0x6c420c0000000000,0x6c420c0000000000},//44
		{0x45660c0000000000,0x45660c0000000000,0x5b480c0000000000,0x5b480c0000000000,0x59670c0000000000,0x59670c0000000000},//45
		{0x68470c0000000000,0x496a0c0000000000,0x495c0c0000000000,0x414a0c0000000000,0x416b0c0000000000,0x685a0c0000000000},//46
	//mmm type groups
		{0xec4b000000000000,0xec4b000000000000,0xec4b000000000000,0xec4b000000000000,0xec4b000000000000,0xec4b000000000000},//47
		{0x5744650000000000,0x5744650000000000,0x5744650000000000,0x5744650000000000,0x5744650000000000,0x5744650000000000},//48
		{0xc26b000000000000,0xc26b000000000000,0xec47000000000000,0xec47000000000000,0xcb6a000000000000,0xcb6a000000000000},//49
		{0x4566580000000000,0x4566580000000000,0x5b48690000000000,0x5b48690000000000,0x5967460000000000,0x5967460000000000},//50
	//only verified to here
		{0xe84b000000000000,0xec49000000000000,0xc942000000000000,0xcb41000000000000,0xec41000000000000,0xe842000000000000},//51
		{0xda44000000000000,0xdc46000000000000,0xea57000000000000,0xd86b000000000000,0xd84a000000000000,0xc757000000000000},//52
		{0xec59000000000000,0xcb5b000000000000,0xcb66000000000000,0xc267000000000000,0xc248000000000000,0xec45000000000000},//53
		{0xda67000000000000,0xdc6b000000000000,0xc55a000000000000,0xc759000000000000,0xea5b000000000000,0xc766000000000000},//54
		{0xc245000000000000,0xc245000000000000,0xec5b000000000000,0xec5b000000000000,0xcb59000000000000,0xcb59000000000000},//55
		{0xd867000000000000,0xd867000000000000,0xc559000000000000,0xc559000000000000,0xe65b000000000000,0xe65b000000000000},//56
		{0xc16a000000000000,0xc147000000000000,0xe859000000000000,0xe845000000000000,0xc95a000000000000,0xc95b000000000000},//57
		{0xc244000000000000,0xc244000000000000,0xec57000000000000,0xec57000000000000,0xcb57000000000000,0xcb57000000000000},//58
		{0x6c4b580000000000,0x6c4b580000000000,0x4b42690000000000,0x4b42690000000000,0x6c42460000000000,0x6c42460000000000},//59
		{0xe657000000000000,0xc557000000000000,0xda48000000000000,0xc75b000000000000,0xd96a000000000000,0xdc67000000000000},//60
		{0xd966000000000000,0xc55b000000000000,0xd966000000000000,0xc55b000000000000,0xd966000000000000,0xc55b000000000000},//61
		{0xc959000000000000,0xe85b000000000000,0xc166000000000000,0xc957000000000000,0xe857000000000000,0xc145000000000000},//62
		{0xec4a0d0000000000,0xcb6b0d0000000000,0xe8105a0000000000,0xe810470000000000,0xc90e6a0000000000,0xec0e5c0000000000},//63
		{0xec59480d00000000,0xcb59670d00000000,0xcb59106700000000,0xc245106700000000,0xc2450e4800000000,0xec450e4800000000},//64
		{0xec0d420000000000,0xec0d420000000000,0xec0d420000000000,0xec0d420000000000,0xec0d420000000000,0xec0d420000000000},//65
		{0xc20d4a0000000000,0xc20d4a0000000000,0xec47100000000000,0xec47100000000000,0xcb6a0e0000000000,0xcb6a0e0000000000},//66
		{0xec495a0d00000000,0xe84b5a0d00000000,0xcb416a1000000000,0xc9426a1000000000,0xe842470e00000000,0xec41470e00000000},//67
		{0x595b0d6748000000,0x595b0d6748000000,0x4559666710000000,0x4559666710000000,0x45665b480e000000,0x45665b480e000000},//68
		{0xec4b100e00000000,0xec4b100e00000000,0xec4b100e00000000,0xec4b100e00000000,0xec4b100e00000000,0xec4b100e00000000},//69
		{0x6d436e0e10000000,0x6d436e0e10000000,0x6d436e0e10000000,0x6d436e0e10000000,0x6d436e0e10000000,0x6d436e0e10000000},//70
		{0xec0c420000000000,0xec0c420000000000,0xec0c420000000000,0xec0c420000000000,0xec0c420000000000,0xec0c420000000000},//71
		{0xc20c660000000000,0xc20c660000000000,0xec0c480000000000,0xec0c480000000000,0xcb0c670000000000,0xcb0c670000000000},//72
		{0xd90c480000000000,0xc50c670000000000,0xd90c480000000000,0xc50c670000000000,0xd90c480000000000,0xc50c670000000000},//73
		{0xec49590000000000,0xe84b5b0000000000,0xcb41660000000000,0xc942670000000000,0xe842480000000000,0xec41450000000000},//74
	};

	//@brief    : compress an EMsoft generator string into a 64 bit integer
	//@param gen: 40 character EMsoft generator string to compress
	//@return   : uint64_t representation of generator string
	uint64_t encode(std::string gen) {
		//sanity check input
		if(40 != gen.size()) throw std::runtime_error("invalid generator string (must be 40 characters long)");
		if(!('0' == gen[0] || '1' == gen[0])) throw std::runtime_error("invalid generator string (first character must be '0' or '1')");
		if(!std::isdigit(gen[1])) throw std::runtime_error("invalid generator string (second character must be a number)");

		//now loop over 4 character substrings accumulating index
		const size_t num = gen[1] - '0';//C++ guarentees 0-9 are contigous
		uint64_t encoded = 0;
		for(size_t i = 0; i < num; i++) {
			std::string sub(gen.cbegin() + 2 + i * 4, gen.cbegin() + 6 + i * 4);//break out 4 character substring
			const uint64_t idx = std::distance(GenLut, std::find_if(GenLut, GenLut + NumBlocks, [&sub](const char* str) {return sub == str;}));//get substring index
			if(idx == NumBlocks) throw std::runtime_error("block `" + sub + "' not found in library");//make sure it was actually found
			encoded |= idx << (56 - 8 * i);//accumulate indices
		}

		//save inversion flag and return
		if('1' == gen[0]) encoded |= 0x8000000000000000;
		return encoded;
	}

	//@brief    : decompress an EMsoft generator string from a 64 bit integer
	//@param enc: compressed generator string
	//@return   : EMsoft generator string
	std::string decode(uint64_t enc) {
		//build empty string and set inversion flag
		std::string gen(40, ' ');
		gen[0] = ((0x8000000000000000 & enc) == 0x8000000000000000) ? '1' : '0';

		//build up generator string
		for(int i = 7; i >= 0; i--) {
			uint_fast8_t idx = 0x7F & enc;//extract index
			std::copy(GenLut[idx], GenLut[idx]+4, (char*)gen.data() + 2 + 4 * i);//copy 4 character sequence into output
			if(' ' == gen[1] && 0 != idx) {//this is the first non-zero character
				if('1' == GenLut[idx][0]) {//this is an alternate origin
					gen[1] = '0' + i;//just save generator count
				} else {//this is the first empty block
					gen[2+4*i+4] = '0';//0 terminate
					gen[1] = '1' + i;//save generator count
				}
			}
			enc >>= 8;//shift for next byte
		}
		if(' ' == gen[1]) gen[1] = gen[2] = '0';//handle space group 1
		return gen;
	}

	//@brief    : build generator matricies from compressed EMsoft generator
	//@param enc: compressed generator string
	//@param alt: should the alternate origin be used
	//@return   : 4x4 generator matricies (empty if alt = true but there is only 1 origin)
	std::vector<xtal::GenPos> gen_from_enc(uint64_t enc, const bool alt) {
		const bool hasInv = (0x8000000000000000 & enc) == 0x8000000000000000;
		uint_fast8_t bytes[8];
		for(size_t i = 0; i < 8; i++) {
			bytes[7-i] = 0x7F & enc;//get last byte
			enc >>= 8;//shift over 1 byte
		}

		//seed generators with identity (+inversion if needed)
		std::vector<xtal::GenPos> gen(1, xtal::GenPos::Identity());//vector of generators
		if(hasInv) gen.push_back(xtal::GenPos::Inversion());

		//convert bytes to strings and count valid number
		int8_t ori[3] = {0,0,0};//origin shift
		for(int i = 0; i < 8; i++) {//up to 8 generators
			xtal::GenPos p = xtal::GenPos::Identity();
			if(0x00 == bytes[i]) break;//we've used all the valid generators
			if('1' == GenLut[bytes[i]][0]) {//there is an alternate origin to parse
				p.fromEMsoft(GenLut[bytes[i]]);
				p.getTrans(ori);
			} else {//there is a regular element to parse
				p.fromEMsoft(GenLut[bytes[i]]);//parse
				gen.push_back(p);
			}
		}

		//update origin if needed and return
		if(alt) {
			if(0 == ori[0] && 0 == ori[1] && 0 == ori[2]) gen.clear();//there is no alternate origin
			for(xtal::GenPos& p : gen) p = p.shiftOrigin(ori);
		}
		return gen;
	}

	//@brief    : build generator matricies for a space group
	//@param sg : space group number [1,230]
	//@param alt: should an alternate setting be selected (rhombohedral instead of hex or origin choice 2 instead of 1 as appropriate)
	//@return   : 4x4 generator matricies (empty if alt = true but there is only 1 setting)
	std::vector<xtal::GenPos> gen_from_num(size_t sg, const bool alt) {
		if(sg < 1 || sg > 230) throw std::runtime_error("space group number must be [1,230]");

		//update index for rhombohedral settings if needed
		bool useOri2 = false;
		if(alt) {
			switch(sg) {
				case 146: sg = 231; break;
				case 148: sg = 232; break;
				case 155: sg = 233; break;
				case 160: sg = 234; break;
				case 161: sg = 235; break;
				case 166: sg = 236; break;
				case 167: sg = 237; break;
				default: useOri2 = alt;
			}
		}
		return gen_from_enc(SGLut[sg-1], useOri2);
	}

	//@brief     : build generator matricies for an extended monoclinic space group
	//@param sg  : space group number [3,15]
	//@param cell: cell choice [1,3]
	//@param axs : unique axis {"b", "-b", "c", "-c", "a", or "-a"}
	//@return   : 4x4 generator matricies (empty if alt = true but there is only 1 setting)
	std::vector<xtal::GenPos> mono_gen_from_num(size_t sg, size_t cell, std::string axs) {
		if(sg   < 3 || sg   > 15) throw std::runtime_error("monoclinic space group number must be [3,15]");
		if(cell < 1 || cell > 3 ) throw std::runtime_error("monoclinic cell must be [1,3]");
		uint64_t enc = 0;
		if     (axs ==  "b") enc = MonoLut[sg-3][0][cell-1];
		else if(axs == "-b") enc = MonoLut[sg-3][1][cell-1];
		else if(axs ==  "c") enc = MonoLut[sg-3][2][cell-1];
		else if(axs == "-c") enc = MonoLut[sg-3][3][cell-1];
		else if(axs ==  "a") enc = MonoLut[sg-3][4][cell-1];
		else if(axs == "-a") enc = MonoLut[sg-3][5][cell-1];
		else throw std::runtime_error("monoclinic axis must be {b, -b, c, -c, a, -a}");
		return gen_from_enc(enc << 32);
	}

	//@brief     : build generator matricies for an extended orthorhombic space group
	//@param sg  : space group number [16,74]
	//@param axs : axis configuration {"abc", "bac", "cab", "cba", "bca", or "acb"}
	//@return   : 4x4 generator matricies (empty if alt = true but there is only 1 setting)
	//@warning  : unverified for space groups [51,74] and currently origin 1 only
	std::vector<xtal::GenPos> ortho_gen_from_num(size_t sg, std::string axs) {
		if(sg < 16 || sg > 74) throw std::runtime_error("orthorhombic space group number must be [3,15]");
		uint64_t enc = 0;
		if     (axs == "abc") enc = OrthoLut[sg-16][0];
		else if(axs == "bac") enc = OrthoLut[sg-16][1];
		else if(axs == "cab") enc = OrthoLut[sg-16][2];
		else if(axs == "cba") enc = OrthoLut[sg-16][3];
		else if(axs == "bca") enc = OrthoLut[sg-16][4];
		else if(axs == "acb") enc = OrthoLut[sg-16][5];
		else throw std::runtime_error("orthorhombic axis must be {abc, bac, cab, cba, bca, acb}");
		return gen_from_enc(enc);
	}
}


#endif//_emgen_h_