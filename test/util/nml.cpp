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

#include <iostream>

namespace nml {
	//@brief : check if timer works
	//@return: true / false if tests pass / fail
	bool testNML(std::ostream& os);
}

int main() {
	try {
		return nml::testNML(std::cout) ? EXIT_SUCCESS : EXIT_FAILURE;
	} catch(std::exception& e) {
		std::cout << "caught: " << e.what();
		return EXIT_FAILURE;
	}
}

#include <thread>
#include <vector>
#include <numeric>
#include <cmath>

#include "util/nml.hpp"

namespace nml {

	//@brief : check if timer works
	//@return: true / false if tests pass / fail
	bool testNML(std::ostream& os) {
		//make a good namelist with one of each type
		std::string testList;
		testList += "! make sure we test a comment\n";

		// add some booleans
		testList += " vTrue  = .true. ,\n";
		testList += " vFalse = .false.,\n";

		//add some integers
		testList += " vInt    =  12345,\n";
		testList += " vIntPos = +12345,\n";
		testList += " vIntNeg = -12345,\n";

		//add some doubles
		testList += " vDoub    =  1.2345,\n";
		testList += " vDoubPos = +1.2345,\n";
		testList += " vDoubNeg = -1.2345,\n";

		//add some special doubles
		const bool testNanInf = false;//nan and inf don't work on all platforms
		if(testNanInf) {
			testList += " vDoubNan = nan    ,\n";
			testList += " vDoubInf = inf    ,\n";
		}
		testList += " vDoubSci =  1.23e4,\n";

		testList += "! mix a comment in the middle\n";

		//add some strings
		testList += " vStr   = 'str',\n";
		testList += " vStrSgl = 'str \\'with single quotes\\'',\n";
		testList += " vStrDbl = 'str \"with single quotes\"',\n";

		//add a list of booleans
		testList += " vBools   = .true., .false., .false.,\n";
		testList += " vInts    = 1, 2, 3 , 4,\n";
		testList += " vDoubles = 1, 2, 3., 4,\n";//a mix of doubles/integers should be cast up to doubles
		testList += " vStrs    = 'abc', '123', 'XYZ', '!@#',\n";

		////////////////////////////////////////////////////////////////////////
		//                                                                    //
		//             make sure a good namelist works correctly              //
		//                                                                    //
		////////////////////////////////////////////////////////////////////////

		//wrap namelist and extract values
		std::istringstream is(testList);
		NameList nm;
		nm.read(is);

		os << "checking scalar namelist parsing\n";

		//check bool parsing
		if(!nm.getBool("vTrue") || nm.getBool("vFalse")) {
			os << "failed to extract bool from namelist\n";
			return false;
		}

		//check in parsing
		if(12345 != nm.getInt("vInt") || 12345 != nm.getInt("vIntPos") || -12345 != nm.getInt("vIntNeg")) {
			os << "failed to extract int from namelist\n";
			return false;
		}

		//check normal double parsing
		if(1.2345 != nm.getDouble("vDoub") || 1.2345 != nm.getDouble("vDoubPos") ||-1.2345 != nm.getDouble("vDoubNeg") ) {
			os << "failed to extract double from namelist\n";
			return false;
		}

		//check special double parsing
		const bool nanParsed = testNanInf ? std::isnan(nm.getDouble("vDoubNan")) : true;
		const bool infParsed = testNanInf ? std::isinf(nm.getDouble("vDoubInf")) : true;
		if(!nanParsed || !infParsed || 1.23e4 != nm.getDouble("vDoubSci")) {
			os << "failed to extract special double from namelist\n";
			return false;
		}

		//check string parsing
		if("str"                        != nm.getString("vStr"   ) ||
		   "str 'with single quotes'"   != nm.getString("vStrSgl") ||
		   "str \"with single quotes\"" != nm.getString("vStrDbl") ) {
			os << "failed to extract string from namelist\n";
			return false;
		}

		//first check that a scalar parsed as a vector is correct
		os << "checking scalar vector namelist parsing\n";
		std::vector<bool> bVec = nm.getBools("vTrue");
		std::vector<int> iVec = nm.getInts("vInt");
		std::vector<double> dVec = nm.getDoubles("vDoub");
		std::vector<std::string> sVec = nm.getStrings("vStr");
		if(bVec.size() != 1 || iVec.size() != 1 || dVec.size() != 1 || sVec.size() != 1) {
			os << "failed to extract scalar vectors from namelist\n";
			return false;
		}

		//next check that vectors are parsed correctly
		os << "checking scalar vector namelist parsing\n";
		bVec = nm.getBools("vBools");
		iVec = nm.getInts("vInts");
		dVec = nm.getDoubles("vDoubles");
		sVec = nm.getStrings("vStrs");
		if(3 != bVec.size() || 4 != iVec.size() || 4 != dVec.size() || 4 != sVec.size()) {
			os << "failed to extract vector lengthss from namelist\n";
			return false;
		}

		//check bools parsing
		if(!bVec[0] || bVec[1] || bVec[2]) {
			os << "failed to extract bool vec from namelist\n";
			return false;
		}

		//check ints parsing
		if(1 != iVec[0] || 2 != iVec[1] || 3 != iVec[2] || 4 != iVec[3]) {
			os << "failed to extract int vec from namelist\n";
			return false;
		}

		//check doubles parsing
		if(1.0 != dVec[0] || 2.0 != dVec[1] || 3.0 != dVec[2] || 4.0 != dVec[3]) {
			os << "failed to extract double vec from namelist\n";
			return false;
		}

		//check strings parsing
		if("abc" != sVec[0] || "123" != sVec[1] || "XYZ" != sVec[2] || "!@#" != sVec[3]) {
			os << "failed to extract string vec from namelist\n";
			return false;
		}

		////////////////////////////////////////////////////////////////////////
		//                                                                    //
		//             make sure partial parsing detection works              //
		//                                                                    //
		////////////////////////////////////////////////////////////////////////

		os << "checking scalar partial parsing detection\n";

		//check full parsing flags on last file (false positive)
		if(!nm.fullyParsed() || !nm.unusedTokens().empty()) {
			os << "fully parsed check failed for fully parsed namelist\n";
			return false;
		}

		//build a new file and partially parse
		testList = "placeholder\n tokenOne = 1,\n tokenTwo = 2";
		is.str(testList);//update underlying string
		is.clear();//clear bad bit
		nm.read(is);
		nm.getInt("tokenOne");

		//make sure partial detection works (false negative)
		if(nm.fullyParsed() || "tokentwo" != nm.unusedTokens()) {
			os << "fully parsed check failed for partially parsed namelist\n";
			return false;
		}

		////////////////////////////////////////////////////////////////////////
		//                                                                    //
		//                   test some potential edge cases                   //
		//                                                                    //
		////////////////////////////////////////////////////////////////////////

		os << "checking potential edge cases\n";

		//make sure the first line isn't silently ignored
		try {
			testList = " key = 1\n";
			is.str(testList);//update underlying string
			is.clear();//clear bad bit
			nm.read(is);
			os << "silently ignored first line with a value\n";
			return false;
		} catch (...) {}

		//missing comma
		try {
			testList = "placeholder\n key = 1\n key2 = 2";
			is.str(testList);//update underlying string
			is.clear();//clear bad bit
			nm.read(is);
			os << "failed to detect missing comma at end of line\n";
			return false;
		} catch (...) {}

		//missing leading space
		try {
			testList = "placeholder\nkey = 1\n";
			is.str(testList);//update underlying string
			is.clear();//clear bad bit
			nm.read(is);
			os << "failed to detect missing leading space\n";
			return false;
		} catch (...) {}

		//duplicate key
		try {
			testList = "placeholder\n key = 1,\n key = 2";
			is.str(testList);//update underlying string
			is.clear();//clear bad bit
			nm.read(is);
			os << "failed to detect duplicate key\n";
			return false;
		} catch (...) {}

		//missing '='
		try {
			testList = "placeholder\n key = 1,\n key2 2";
			is.str(testList);//update underlying string
			is.clear();//clear bad bit
			nm.read(is);
			os << "failed to detect missing '='\n";
			return false;
		} catch (...) {}

		//missing string delimiter
		try {
			testList = "placeholder\n key = '1' '2'\n";
			is.str(testList);//update underlying string
			is.clear();//clear bad bit
			nm.read(is);
			os << "failed to detect missing string delimiter\n";
			return false;
		} catch (...) {}

		//bad string opening in list
		try {
			testList = "placeholder\n key = '1', 2\n";
			is.str(testList);//update underlying string
			is.clear();//clear bad bit
			nm.read(is);
			os << "failed to detect bad string opening in list\n";
			return false;
		} catch (...) {}

		//double quoted string
		try {
			testList = "placeholder\n key = \"1\"\n";
			is.str(testList);//update underlying string
			is.clear();//clear bad bit
			nm.read(is);
			os << "failed to detect double quoted string\n";
			return false;
		} catch (...) {}

		//unquoted string
		try {
			testList = "placeholder\n key = value\n";
			is.str(testList);//update underlying string
			is.clear();//clear bad bit
			nm.read(is);
			os << "failed to detect unquoted string\n";
			return false;
		} catch (...) {}

		//int/bool mix
		try {
			testList = "placeholder\n key = 1, .true.\n";
			is.str(testList);//update underlying string
			is.clear();//clear bad bit
			nm.read(is);
			os << "failed to detect int/bool mix\n";
			return false;
		} catch (...) {}

		//double/bool mix
		try {
			testList = "placeholder\n key = 1.0, .true.\n";
			is.str(testList);//update underlying string
			is.clear();//clear bad bit
			nm.read(is);
			os << "failed to detect double/bool mix\n";
			return false;
		} catch (...) {}

		//(string mixing is allowed)

		//if we made it this far all tests passed
		return true;
	}
}
