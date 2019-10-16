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

//@brief : check if base64 encode/decode is self consistent
//@return: true / false if tests pass / fail
bool testBase64(std::ostream& os);

int main() {
	try {
		return testBase64(std::cout) ? EXIT_SUCCESS : EXIT_FAILURE;
	} catch(std::exception& e) {
		std::cout << "caught: " << e.what();
		return EXIT_FAILURE;
	}
}

#include "util/base64.hpp"

#include <random>
#include <sstream>
#include <algorithm>
#include <limits>

//@brief : check if base64 encode/decode is self consistent
//@return: true / false if tests pass / fail
bool testBase64(std::ostream& os) {
	//make sure a few different strings have the expected encoding
	std::string test[3] = {
		"test string length % 3 == 0"  ,
		"test string length % 3 == 1  ",
		"test string length % 3 == 2 " ,
	};

	//expected encodings
	std::string enc[3] = {
		"dGVzdCBzdHJpbmcgbGVuZ3RoICUgMyA9PSAw"    ,
		"dGVzdCBzdHJpbmcgbGVuZ3RoICUgMyA9PSAxICA=",
		"dGVzdCBzdHJpbmcgbGVuZ3RoICUgMyA9PSAyIA==",
	};

	//loop over test strings
	for(size_t i = 0; i < 3; i++) {
		//make sure encoding matches expected data
		std::stringstream ss;
		size_t len = base64::encode(test[i].data(), test[i].size(), ss);//encode
		if(0 != enc[i].compare(ss.str())) {//make sure encoding matches our expected value
			os << "expected `" << test[i] << "' to encode as:\n\t`" << enc[i] << "'\nbut got\n\t`" << ss.str() << "'\n";
			return false;
		}
		if(len != ss.str().size()) {//make sure the returned size is correct
			os << "encoded length is " << ss.str().size() << " but encoding returned " << len << '\n';
			return false;
		}

		//make sure decoding matches expected data
		ss.str(std::string());
		len = base64::decode(enc[i].data(), enc[i].size(), ss);//decode
		if(0 != test[i].compare(ss.str())) {//make sure encoding matches our expected value
			os << "expected `" << enc[i] << "' to decode as:\n\t`" << test[i] << "'\nbut got\n\t`" << ss.str() << "'\n";
			return false;
		}
		if(len != ss.str().size()) {//make sure the returned size is correct
			os << "encoded length is " << ss.str().size() << " but encoding returned " << len << '\n';
			return false;
		}
	}

	//now test for round trip self consistency on random data

	//generate some random data
	char buff[64];
	std::default_random_engine gen;
	std::uniform_int_distribution<short> dist(std::numeric_limits<char>::min(), std::numeric_limits<char>::max());
	for(size_t i = 0; i < sizeof(buff); i++) buff[i] = (char)dist(gen);

	//loop over all termination types (3*n+ 0,1, or 2)
	for(size_t i = 0; i < 3; i++) {
		//encode the data
		const size_t bytes = sizeof(buff) - i;
		std::stringstream ssEnc, ssDec;
		size_t len = base64::encode(buff, bytes, ssEnc);
		if(len != ssEnc.str().size()) {
			os << "encoded length is " << ssEnc.str().size() << " but encoding returned " << len << '\n';
			return false;
		}

		//decode the data
		len = base64::decode(ssEnc.str().data(), ssEnc.str().size(), ssDec);
		if(len != ssDec.str().size()) {
			os << "encoded length is " << ssDec.str().size() << " but encoding returned " << len << '\n';
			return false;
		}

		//sanity check output vs input length
		if(len != bytes) {
			os << "round trip decoding of " << bytes << " bytes became " << len << " bytes\n";
			return false;
		}

		//make sure we got back what we started with
		if(!std::equal(buff, buff + len, ssDec.str().c_str())) {
			os << "round trip encode/decode failed for " << bytes << " bytes encoded as:\n\t`" << ssEnc.str() << "'\n";
			return false;
		}
	}

	//if we made it this far all tests passed
	os << "all tests passed\n";
	return true;
}
