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

#include <cstdint>

//@brief  : function to undo encoding of generator string into 64bit int
//@param v: generator string encoded as 64 bit integer
//@param s: location to write generator string (35 characters = 34 + NULL)
//@note   : look up table generated with this program
void decodeLut(const uint64_t v, char * s) {
	//look up table for 4 character generator code substrings
	static const char lut[88][5] = {
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
		"lOOD", "lOOO", "mOOO", "nOOC", "nOOE", "nOOO"
	};

	//break apart bytes
	uint8_t bytes[8] = {
		(uint8_t) ((v >> 56) & 0x7F),
		(uint8_t) ((v >> 48) & 0x7F),
		(uint8_t) ((v >> 40) & 0x7F),
		(uint8_t) ((v >> 32) & 0x7F),
		(uint8_t) ((v >> 24) & 0x7F),
		(uint8_t) ((v >> 16) & 0x7F),
		(uint8_t) ((v >>  8) & 0x7F),
		(uint8_t) ((v      ) & 0x7F),
	};

	//convert bytes to strings and count valid number
	int count = 0;
	bool hasAlt = false;
	for(int i = 0; i < 8; i++) {
		char* p = s + 2 + i * 4;//get pointer to sequence start
		for(int j = 0; j < 4; j++) p[j] = lut[bytes[i]][j];//copy string
		if(0x00 != bytes[i]) ++count;//accumulate number of generators
		if('1' == lut[bytes[i]][0]) {//check for an alternate origin
			if(hasAlt) {//this shouldn't be possible
				for(int i = 0; i < 34; i++) s[i] = 0;//clear string
				return;
			}
			hasAlt = true;//flag alternate origin
			--count;//doun't include in count
		}
	}

	//set inversion flag
	s[0] = 0x8000000000000000 & v ? '1' : '0';

	//set count
	s[1] = '0' + count;

	//append with 0 if needed
	if(!hasAlt) s[2 + 4 * count] = '0';
}

//@brief: look up table for generators for each space group such that decodeLut(SGLut[i+1], s) writes the generator string for space group i in s
//@note : generated with this program
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

#include <string>
#include <vector>

//the original space group generator strings
const std::vector<std::string> generators = {
	"000                               ", "100                               ", "01cOOO0                           ",
	"01cODO0                           ", "02aDDOcOOO0                       ", "01jOOO0                           ",
	"01jOOD0                           ", "02aDDOjOOO0                       ", "02aDDOjOOD0                       ",
	"11cOOO0                           ", "11cODO0                           ", "12aDDOcOOO0                       ",
	"11cOOD0                           ", "11cODD0                           ", "12aDDOcOOD0                       ",
	"02bOOOcOOO0                       ", "02bOODcOOD0                       ", "02bOOOcDDO0                       ",
	"02bDODcODD0                       ", "03aDDObOODcOOD0                   ", "03aDDObOOOcOOO0                   ",
	"04aODDaDODbOOOcOOO0               ", "03aDDDbOOOcOOO0                   ", "03aDDDbDODcODD0                   ",
	"02bOOOjOOO0                       ", "02bOODjOOD0                       ", "02bOOOjOOD0                       ",
	"02bOOOjDOO0                       ", "02bOODjDOO0                       ", "02bOOOjODD0                       ",
	"02bDODjDOD0                       ", "02bOOOjDDO0                       ", "02bOODjDDO0                       ",
	"02bOOOjDDD0                       ", "03aDDObOOOjOOO0                   ", "03aDDObOODjOOD0                   ",
	"03aDDObOOOjOOD0                   ", "03aODDbOOOjOOO0                   ", "03aODDbOOOjODO0                   ",
	"03aODDbOOOjDOO0                   ", "03aODDbOOOjDDO0                   ", "04aODDaDODbOOOjOOO0               ",
	"04aODDaDODbOOOjBBB0               ", "03aDDDbOOOjOOO0                   ", "03aDDDbOOOjDDO0                   ",
	"03aDDDbOOOjDOO0                   ", "12bOOOcOOO0                       ", "03bOOOcOOOhDDD1BBB                ",
	"12bOOOcOOD0                       ", "03bOOOcOOOhDDO1BBO                ", "12bDOOcOOO0                       ",
	"12bDOOcDDD0                       ", "12bDODcDOD0                       ", "12bDOOcOOD0                       ",
	"12bOOOcDDO0                       ", "12bDDOcODD0                       ", "12bOODcODD0                       ",
	"12bOOOcDDD0                       ", "03bOOOcDDOhDDO1BBO                ", "12bDDDcOOD0                       ",
	"12bDODcODD0                       ", "12bDODcODO0                       ", "13aDDObOODcOOD0                   ",
	"13aDDObODDcODD0                   ", "13aDDObOOOcOOO0                   ", "13aDDObOOOcOOD0                   ",
	"13aDDObODOcODO0                   ", "04aDDObDDOcOOOhODD1OBB            ", "14aODDaDODbOOOcOOO0               ",
	"05aODDaDODbOOOcOOOhBBB1ZZZ        ", "13aDDDbOOOcOOO0                   ", "13aDDDbOOOcDDO0                   ",
	"13aDDDbDODcODD0                   ", "13aDDDbODOcODO0                   ", "02bOOOgOOO0                       ",
	"02bOODgOOB0                       ", "02bOOOgOOD0                       ", "02bOODgOOF0                       ",
	"03aDDDbOOOgOOO0                   ", "03aDDDbDDDgODB0                   ", "02bOOOmOOO0                       ",
	"03aDDDbOOOmOOO0                   ", "12bOOOgOOO0                       ", "12bOOOgOOD0                       ",
	"03bOOOgDDOhDDO1YBO                ", "03bOOOgDDDhDDD1YYY                ", "13aDDDbOOOgOOO0                   ",
	"04aDDDbDDDgODBhODB1OYZ            ", "03bOOOgOOOcOOO0                   ", "03bOOOgDDOcDDO0                   ",
	"03bOODgOOBcOOO0                   ", "03bOODgDDBcDDB0                   ", "03bOOOgOODcOOO0                   ",
	"03bOOOgDDDcDDD0                   ", "03bOODgOOFcOOO0                   ", "03bOODgDDFcDDF0                   ",
	"04aDDDbOOOgOOOcOOO0               ", "04aDDDbDDDgODBcDOF0               ", "03bOOOgOOOjOOO0                   ",
	"03bOOOgOOOjDDO0                   ", "03bOOOgOODjOOD0                   ", "03bOOOgDDDjDDD0                   ",
	"03bOOOgOOOjOOD0                   ", "03bOOOgOOOjDDD0                   ", "03bOOOgOODjOOO0                   ",
	"03bOOOgOODjDDO0                   ", "04aDDDbOOOgOOOjOOO0               ", "04aDDDbOOOgOOOjOOD0               ",
	"04aDDDbDDDgODBjOOO0               ", "04aDDDbDDDgODBjOOD0               ", "03bOOOmOOOcOOO0                   ",
	"03bOOOmOOOcOOD0                   ", "03bOOOmOOOcDDO0                   ", "03bOOOmOOOcDDD0                   ",
	"03bOOOmOOOjOOO0                   ", "03bOOOmOOOjOOD0                   ", "03bOOOmOOOjDDO0                   ",
	"03bOOOmOOOjDDD0                   ", "04aDDDbOOOmOOOjOOO0               ", "04aDDDbOOOmOOOjOOD0               ",
	"04aDDDbOOOmOOOcOOO0               ", "04aDDDbOOOmOOOcDOF0               ", "13bOOOgOOOcOOO0                   ",
	"13bOOOgOOOcOOD0                   ", "04bOOOgOOOcOOOhDDO1YYO            ", "04bOOOgOOOcOOOhDDD1YYY            ",
	"13bOOOgOOOcDDO0                   ", "13bOOOgOOOcDDD0                   ", "04bOOOgDDOcDDOhDDO1YBO            ",
	"04bOOOgDDOcDDDhDDO1YBO            ", "13bOOOgOODcOOO0                   ", "13bOOOgOODcOOD0                   ",
	"04bOOOgDDDcOODhDDD1YBY            ", "04bOOOgDDDcOOOhDDD1YBY            ", "13bOOOgOODcDDO0                   ",
	"13bOOOgDDDcDDD0                   ", "04bOOOgDDDcDDDhDDD1YBY            ", "04bOOOgDDDcDDOhDDD1YBY            ",
	"14aDDDbOOOgOOOcOOO0               ", "14aDDDbOOOgOOOcOOD0               ", "05aDDDbDDDgODBcDOFhODB1OBZ        ",
	"05aDDDbDDDgODBcDOBhODB1OBZ        ", "01nOOO0                           ", "01nOOC0                           ",
	"01nOOE0                           ", "02aECCnOOO0                       ", "11nOOO0                           ",
	"12aECCnOOO0                       ", "02nOOOfOOO0                       ", "02nOOOeOOO0                       ",
	"02nOOCfOOE0                       ", "02nOOCeOOO0                       ", "02nOOEfOOC0                       ",
	"02nOOEeOOO0                       ", "03aECCnOOOeOOO0                   ", "02nOOOkOOO0                       ",
	"02nOOOlOOO0                       ", "02nOOOkOOD0                       ", "02nOOOlOOD0                       ",
	"03aECCnOOOkOOO0                   ", "03aECCnOOOkOOD0                   ", "12nOOOfOOO0                       ",
	"12nOOOfOOD0                       ", "12nOOOeOOO0                       ", "12nOOOeOOD0                       ",
	"13aECCnOOOeOOO0                   ", "13aECCnOOOeOOD0                   ", "02nOOObOOO0                       ",
	"02nOOCbOOD0                       ", "02nOOEbOOD0                       ", "02nOOEbOOO0                       ",
	"02nOOCbOOO0                       ", "02nOOObOOD0                       ", "02nOOOiOOO0                       ",
	"12nOOObOOO0                       ", "12nOOObOOD0                       ", "03nOOObOOOeOOO0                   ",
	"03nOOCbOODeOOC0                   ", "03nOOEbOODeOOE0                   ", "03nOOEbOOOeOOE0                   ",
	"03nOOCbOOOeOOC0                   ", "03nOOObOODeOOO0                   ", "03nOOObOOOkOOO0                   ",
	"03nOOObOOOkOOD0                   ", "03nOOObOODkOOD0                   ", "03nOOObOODkOOO0                   ",
	"03nOOOiOOOkOOO0                   ", "03nOOOiOODkOOD0                   ", "03nOOOiOOOeOOO0                   ",
	"03nOOOiOODeOOO0                   ", "13nOOObOOOeOOO0                   ", "13nOOObOOOeOOD0                   ",
	"13nOOObOODeOOD0                   ", "13nOOObOODeOOO0                   ", "03bOOOcOOOdOOO0                   ",
	"05aODDaDODbOOOcOOOdOOO0           ", "04aDDDbOOOcOOOdOOO0               ", "03bDODcODDdOOO0                   ",
	"04aDDDbDODcODDdOOO0               ", "13bOOOcOOOdOOO0                   ", "04bOOOcOOOdOOOhDDD1YYY            ",
	"15aODDaDODbOOOcOOOdOOO0           ", "06aODDaDODbOOOcOOOdOOOhBBB1ZZZ    ", "14aDDDbOOOcOOOdOOO0               ",
	"13bDODcODDdOOO0                   ", "14aDDDbDODcODDdOOO0               ", "04bOOOcOOOdOOOeOOO0               ",
	"04bOOOcOOOdOOOeDDD0               ", "06aODDaDODbOOOcOOOdOOOeOOO0       ", "06aODDaDODbODDcDDOdOOOeFBF0       ",
	"05aDDDbOOOcOOOdOOOeOOO0           ", "04bDODcODDdOOOeBFF0               ", "04bDODcODDdOOOeFBB0               ",
	"05aDDDbDODcODDdOOOeFBB0           ", "04bOOOcOOOdOOOlOOO0               ", "06aODDaDODbOOOcOOOdOOOlOOO0       ",
	"05aDDDbOOOcOOOdOOOlOOO0           ", "04bOOOcOOOdOOOlDDD0               ", "06aODDaDODbOOOcOOOdOOOlDDD0       ",
	"05aDDDbDODcODDdOOOlBBB0           ", "14bOOOcOOOdOOOeOOO0               ", "05bOOOcOOOdOOOeOOOhDDD1YYY        ",
	"14bOOOcOOOdOOOeDDD0               ", "05bOOOcOOOdOOOeDDDhDDD1YYY        ", "16aODDaDODbOOOcOOOdOOOeOOO0       ",
	"16aODDaDODbOOOcOOOdOOOeDDD0       ", "07aODDaDODbODDcDDOdOOOeFBFhBBB1ZZZ", "07aODDaDODbODDcDDOdOOOeFBFhFFF1XXX",
	"15aDDDbOOOcOOOdOOOeOOO0           ", "15aDDDbDODcODDdOOOeFBB0           ", "01dOOO0                           ",
	"11dOOO0                           ", "02dOOOfOOO0                       ", "02dOOOlOOO0                       ",
	"02dOOOlDDD0                       ", "12dOOOfOOO0                       ", "12dOOOfDDD0                       " 
};

#include <set>
#include <iostream>
#include <iomanip>

int main() {
	//space to hold the possible 4 character generator substrings e.g. "hDDD"
	std::set<std::string> lib;
	lib.insert("    ");//empty space

	//counter for maximum number of 4 character sub strings across all generators
	size_t maxNum = 0;

	//loop over generators
	for(const std::string& gen : generators) {
		//get number of generator matricies
		size_t num = 0;
		switch(gen[1]) {
			case '0': num = 0; break;
			case '1': num = 1; break;
			case '2': num = 2; break;
			case '3': num = 3; break;
			case '4': num = 4; break;
			case '5': num = 5; break;
			case '6': num = 6; break;
			case '7': num = 7; break;
			case '8': num = 8; break;
			case '9': num = 9; break;
			default: throw std::runtime_error("expected number in 2nd character of generator string");
		}

		//check if there is a special origin shift generator
		if('1' == gen[2 + 4 * num]) ++num;
		if(num > maxNum) maxNum = num;

		//now accumulate all 4 letter substrings into library
		for(size_t i = 0; i < num; i++) lib.insert(std::string(gen.cbegin() + 2 + i * 4, gen.cbegin() + 6 + i * 4));
	}

	//make sure we can still nicely fit everything into a single 64 bit integer (in case some weird QC are added later)
	//at the time of writing there are a maximum of 8 matricies per space group + 86 possibilities so we can nicely convert to 8*8 bits
	//we actually need to fit into 7 bits so there is space for the inversion flag in the first bit
	//that are still 7 unused bits if needed (first bit of each byte after the first)
	if(lib.size() > 128 || maxNum > 8) throw std::runtime_error("cannot convert to 64bit ints");

	//print the library lookup table
	std::cout << "static const char lut[" << lib.size() << "][5] = {";
	size_t count = 0;
	for(const std::string& s : lib) {
		if(0 == count++ % 8) std::cout << "\n\t";
		std::cout << '"' << s << "\", ";
	}
	std::cout << "\n};\n";

	//now loop over space groups again building the space group table
	std::vector<uint64_t> lut;
	for(const std::string& gen : generators) {
		//get number of generator matricies
		size_t num = 0;
		switch(gen[1]) {
			case '0': num = 0; break;
			case '1': num = 1; break;
			case '2': num = 2; break;
			case '3': num = 3; break;
			case '4': num = 4; break;
			case '5': num = 5; break;
			case '6': num = 6; break;
			case '7': num = 7; break;
			case '8': num = 8; break;
			case '9': num = 9; break;
			//we shouldn't need a default here since it would have thrown above
		}

		//check if there is a special origin shift generator
		if('1' == gen[2 + 4 * num]) ++num;

		//now loop over all 4 letter substrings
		uint64_t encoded = 0;
		for(size_t i = 0; i < num; i++) {
			std::string sub(gen.cbegin() + 2 + i * 4, gen.cbegin() + 6 + i * 4);//break out 4 character substring
			encoded |= std::distance(lib.begin(), lib.find(sub)) << (56 - 8 * i);
		}

		//save inversion flag
		if('1' == gen[0]) encoded |= (uint64_t(1) << 63);

		//make sure we can resconstruct the string
		std::string rec("..................................");
		decodeLut(encoded, &rec[0]);
		if(rec != gen) throw std::runtime_error("couldn't round trip decode " + gen);

		//save encoded value
		lut.push_back(encoded);
	}

	//print the space group lookup table
	std::cout << std::hex << std::setfill('0');
	std::cout << "static const uint64_t SGLut[" << lut.size() << "] = {";
	for(size_t i = 0; i < lut.size(); i++) {
		if(0 == i % 6) std::cout << "\n\t";
		std::cout << "0x" << std::setw(16) << lut[i] << ',';
	}
	std::cout << "\n};\n";
}
