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

////////////////////////////////////////////////////////////////////////
//        test program for functions in include/sht/wigner.hpp        //
////////////////////////////////////////////////////////////////////////


#include <iostream>

namespace xtal {
	//@brief   : test if the Hermann-Maguin name class works
	//@param os: location to write status
	//@return  : true/false if tests passed/failed
	bool testHM(std::ostream& os);

	//@brief   : test if the Hermann-Maguin name class parses names correctly
	//@param os: location to write status
	//@return  : true/false if tests passed/failed
	bool testParse(std::ostream& os);

	//@brief   : test if the Hermann-Maguin name builds default generators correctly
	//@param os: location to write status
	//@return  : true/false if tests passed/failed
	bool testGen(std::ostream& os);

	//@brief   : test if the Hermann-Maguin name build alternate generators correctly (rhomobohedral or origin choice 2)
	//@param os: location to write status
	//@return  : true/false if tests passed/failed
	bool testAlt(std::ostream& os);

	//@brief   : test if the Hermann-Maguin name build extended monoclinic generators correctly
	//@param os: location to write status
	//@return  : true/false if tests passed/failed
	bool testMono(std::ostream& os);

	//@brief   : test if the Hermann-Maguin name build extended orthorhombic generators correctly
	//@param os: location to write status
	//@return  : true/false if tests passed/failed
	bool testOrtho(std::ostream& os);
	
}

int main() {
	return xtal::testHM(std::cout) ? EXIT_SUCCESS : EXIT_FAILURE;
}

#include <vector>
#include <string>

#include "emsoft_gen.hpp"
#include "xtal/hm.hpp"

namespace xtal {

	//@brief   : test if the Hermann-Maguin name class works
	//@param os: location to write status
	//@return  : true/false if tests passed/failed
	bool testHM(std::ostream& os) {
		return testParse(os) &&
		       testGen  (os) &&
		       testAlt  (os) &&
		       testMono (os) &&
		       testOrtho(os);
	}

	//@brief   : test if the Hermann-Maguin name class parses names correctly
	//@param os: location to write status
	//@return  : true/false if tests passed/failed
	bool testParse(std::ostream& os) {
		//make sure we can build symbols from short names
		for(size_t i = 0; i < 230; i++) {
			HermannMaguin hmTable(i+1);//get name from number
			HermannMaguin hmName;
			hmTable.clearOrigin();
			hmName .fromString(emsoft::SGNames[i].c_str());
			if(!(hmName.shortSym() == hmTable.shortSym())) {
				os << "didn't properly parse: " << emsoft::SGNames[i] << ":\n";
				os << "table :\t" << hmTable.shortSym().to_string() << '\t' << hmTable.to_string() << '\n';
				os << "parsed:\t" << hmName .shortSym().to_string() << '\t' << hmName .to_string() << '\n';
				return false;
			}
		}
		return true;
	}

	//@brief   : test if the Hermann-Maguin name builds default generators correctly
	//@param os: location to write status
	//@return  : true/false if tests passed/failed
	bool testGen(std::ostream& os) {
		os << "==============================\n";
		os << "checking default settings against EMsoft generator strings\n";
		os << "------------------------------\n";
		for(size_t i = 0; i < 230; i++) {
			HermannMaguin hm(i+1);//get name from space group number
			if(i < 99) os << (i < 9 ? "  " : " ");
			os << i+1 << ": " << hm.to_string() << '\t' << hm.shortSym().to_string() << '\n';//print full and short symbol
			std::vector<GenPos> gen = hm.generators();//build generators from name
			std::vector<GenPos> mats = GenPos::CloseSet(gen);//close set

			std::vector<GenPos> mats2 = GenPos::CloseSet(emsoft::gen_from_num(i+1));//get closed set of EMsoft generators
			if(mats != mats2) {//make sure the generator sets match
				mats = gen;
				os << "old:\n";
				for(auto x : mats2) os << x.to_string() << '\n';
				os << "============\n";
				os << "new:\n";
				for(auto x : mats) os << x.to_string() << '\n';
				return false;
			}
		}
		os << '\n';
		return true;
	}

	//@brief   : test if the Hermann-Maguin name build alternate generators correctly (rhomobohedral or origin choice 2)
	//@param os: location to write status
	//@return  : true/false if tests passed/failed
	bool testAlt(std::ostream& os) {
		//test rhombohedral setting
		os << "==============================\n";
		os << "testing rhombohedral settings against EMsoft generator strings\n";
		os << "------------------------------\n";
		for(size_t i = 0; i < 230; i++) {//loop over all groups
			HermannMaguin hm(i+1);//get name from space group number
			if('R' == hm.to_string().front() ) {//check if this is a rhombohedral group
				os << i+1 << ": " << hm.to_string() << '\t' << hm.shortSym().to_string() << '\n';//print full and short symbol
				std::vector<GenPos> gen = hm.generators(NULL, false);//build generators from name
				std::vector<GenPos> mats = GenPos::CloseSet(gen);//close set
				std::vector<GenPos> mats2 = GenPos::CloseSet(emsoft::gen_from_num(i+1, true));//get closed set of EMsoft generators

				//make sure the generator sets match
				if(mats != mats2) {
					os << "old:\n";
					for(auto x : mats2) os << x.to_string() << '\n';
					os << "============\n";
					os << "new:\n";
					for(auto x : mats) os << x.to_string() << '\n';
					return false;
				}
			}
		}
		os << '\n';

		//test origin choice 2
		os << "==============================\n";
		os << "testing origin choice 2 against EMsoft generator strings\n";
		os << "------------------------------\n";
		for(size_t i = 0; i < 230; i++) {//loop over all groups
			try {//try to build with alternate origin
				HermannMaguin hm;
				hm.fromNumber(i+1, true);//get name from space group number (throws if no alternate origin)
				if(i < 99) os << (i+1 < 9 ? "  " : " ");
				os << i+1 << ": " << hm.to_string() << '\t' << hm.shortSym().to_string() << '\n';//print full and short symbol
				std::vector<GenPos> gen = hm.generators();//build generators from name
				std::vector<GenPos> mats = GenPos::CloseSet(gen);//close set
				std::vector<GenPos> mats2 = GenPos::CloseSet(emsoft::gen_from_num(i+1, true));//Generator(i).closeSet(true);//get closed set of EMsoft generators

				if(mats != mats2) {
					os << "old:\n";
					for(auto x : mats2) os << x.to_string() << '\n';
					// for(auto x : mats2) os << x.to_wyckoff() << '\n';
					os << "============\n";
					os << "new:\n";
					for(auto x : mats) os << x.to_string() << '\n';
					return false;
				}
			} catch(...) {
				//only 1 setting
			}
		}
		os << '\n';
		return true;
	}

	//@brief   : test if the Hermann-Maguin name build extended monoclinic generators correctly
	//@param os: location to write status
	//@return  : true/false if tests passed/failed
	bool testMono(std::ostream& os) {
		//test for monoclinic settings	
		os << "==============================\n";
		os << "testing extended monoclinic symbols\n";
		os << "------------------------------\n";
		const std::string abc[6] = {"b","-b","c","-c","a","-a"};
		os << " # :     ";
		for(size_t i = 0; i < 6; i++) os << "\t  " << abc[i] << "      ";
		os << '\n';
		for(size_t i = 3; i <= 15; i++) {
			HermannMaguin hm(i);//get space group
			//build table of 6 alternate settings in standard order
			HermannMaguin extMono[3][6];
			for(size_t j = 0; j < 3; j++) {
				for(size_t k = 0; k < 6; k++) {
					extMono[j][k] = hm.changeMonoCell(j+1, abc[k]);
					std::vector<GenPos> mats2 = GenPos::CloseSet(emsoft::mono_gen_from_num(i, j+1, abc[k]));
					std::vector<GenPos> gen = extMono[j][k].generators();//build generators from name
					std::vector<GenPos> mats = GenPos::CloseSet(gen);//close set
					if(mats != mats2) {
						os << "extended monoclinic inconsistent with table\n";
						return false;
					}
				}
			}

			//determine if this group has multiple cells
			bool hasCells = false;
			for(size_t k = 0; k < 6; k++) {
				if(extMono[0][k] != extMono[1][k] || extMono[1][k] != extMono[2][k]) {
					hasCells = true;
					break;
				}
			}
			if(i < 99) os << (i < 9 ? "  " : " ");

			//print alternate settings
			if(hasCells) {
				os << i << ":\n";
				for(size_t j = 0; j < 3; j++) {
					os << "     cell " << j+1;
					for(size_t k = 0; k < 6; k++) os << " \t" << extMono[j][k].shortSym().to_string();
					os << '\n';
				}
				os << '\n';
			} else {
				os << i << ":    ";
				for(size_t k = 0; k < 6; k++) os << " \t" << extMono[0][k].shortSym().to_string();
				os << '\n';
			}


			//now make sure round trip conversions are self consistent
			for(size_t jj = 0; jj < 3; jj++) {
				for(size_t kk = 0; kk < 6; kk++) {
					for(size_t j = 0; j < 3; j++) {
						for(size_t k = 0; k < 6; k++) {
							HermannMaguin ext = extMono[jj][kk].changeMonoCell(j+1, abc[k]);
							if(!(ext == extMono[j][k])) throw std::runtime_error("extended monoclinic conversions inconsistent");

						}
					}
				}
			}

			//TODO: check that we can produce extended symbols via transformation matrix


			//TODO: check alternate symbols (e.g. C1a1 instead of C1m1 for 8)
		}
		// os << '\n';
		return true;
	}

	//@brief   : test if the Hermann-Maguin name build extended orthorhombic generators correctly
	//@param os: location to write status
	//@return  : true/false if tests passed/failed
	bool testOrtho(std::ostream& os) {
		//test alternate orthorhombic settings
		os << "==============================\n";
		os << "testing extended orthorhombic symbols\n";
		os << "------------------------------\n";
		const std::string abc[6] = {"abc","bac","cab","cba","bca","acb"};
		os << " # :";
		for(size_t i = 0; i < 6; i++) os << "\t  " << abc[i] << "   ";
		os << '\n';
		for(size_t i = 16; i <= 74; i++) {//loop over orhtorhombic groups printing extended symbols
			HermannMaguin hm(i);//get space group

			//build table of 6 alternate settings in standard order
			HermannMaguin extOrtho[6];
			for(size_t j = 0; j < 6; j++) {
				extOrtho[j] = hm.changeOrthoAxis(abc[j]);
				std::vector<GenPos> mats2 = GenPos::CloseSet(emsoft::ortho_gen_from_num(i, abc[j]));
				std::vector<GenPos> gen = extOrtho[j].generators();//build generators from name
				std::vector<GenPos> mats = GenPos::CloseSet(gen);//close set
				if(mats != mats2) {
					os << "extended orthorhombic inconsistent with table\n";
					return false;
				}
			}

			//print alternate symbols
			if(i < 99) os << (i < 9 ? "  " : " ");
			os << i << ":";
			for(size_t j = 0; j < 6; j++) os << " \t" << extOrtho[j].shortSym().to_string();
			os << '\n';

			//now make sure round trip conversions are self consistent
			for(size_t j = 0; j < 6; j++) {
				for(size_t k = 0; k < 6; k++) {
					HermannMaguin extK = extOrtho[j].changeOrthoAxis(abc[k]);//get round trip conversion
					if(!(extK == extOrtho[k])) throw std::runtime_error("extended orthorhombic conversions inconsistent");
				}
			}

			//TODO: check that we can produce extended symbols via transformation matrix

			//TODO: check alternate symbols
		}
		os << '\n';
		return true;
	}

}
