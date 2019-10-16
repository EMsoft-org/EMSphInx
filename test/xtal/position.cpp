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
//      test program for functions in include/xtal/position.hpp       //
////////////////////////////////////////////////////////////////////////

#include <iostream>

namespace xtal {
	//@brief   : run all general position and generator unit tests
	//@param os: output stream to write error messages to
	//@return  : true/false if the self tests pass/fail
	bool runTests(std::ostream& os);


	//@brief   : test all general position constructors (+ static factories)
	//@param os: output stream to write error messages to
	//@return  : true/false if the self tests pass/fail
	bool testBuild(std::ostream& os);

	//@brief   : test all general position matrix operations
	//@param os: output stream to write error messages to
	//@return  : true/false if the self tests pass/fail
	bool testMat(std::ostream& os);

}


int main() {
	return xtal::runTests(std::cout) ? EXIT_SUCCESS : EXIT_FAILURE;
}


#include "xtal/position.hpp"

namespace xtal {
	//@brief   : run all general position and generator unit tests
	//@param os: output stream to write error messages to
	//@return  : true/false if the self tests pass/fail
	bool runTests(std::ostream& os) {

		//make sure conversion between 3x3 matrix and 64 matrix table is self consistent
		for(uint_fast8_t i = 0; i < 64; i++) {
			uint8_t j = xtal::GenPos::Mat3ToIdx(xtal::GenPos::IdxToMat3(i));
			if(i != j) {
				os << "self consistent matrix <==> index failed for " << i << "\n";
				return false;
			}
		}

		//make sure constructors and matrix operations work
		if(!testBuild(os)) return false;
		if(!testMat  (os)) return false;

		//make sure EMsoft parsing is self consistent for generators
		for(uint32_t i = 0; i < 64; i++) {
			GenPos p(i), q;
			char em[5] = "    ";
			auto trs = {0,4,6,8,12,16,18,20};//allowable translations
			try {
				//loop over all possible translations
				for(int8_t x : trs) {
					for(int8_t y : trs) {
						for(int8_t z : trs) {
							int8_t t[3] = {x, y, z};
							p.setTrans(t);
							p.toEMsoft(em);
							q.fromEMsoft(em);
							if(p != q) {
								os << "EMsoft conversions not self consistent\n";
								return false;
							}
						}
					}
				}
			} catch (...) {
				//not all possibilities can be converted
			}
		}

		//make sure EMsoft parsing is self consistent for origin shifts
		{
			auto trs = {0,15,18,21};//allowable translations
			GenPos p = GenPos::Identity();
			GenPos q;
			char em[5] = "    ";
			for(int8_t x : trs) {
				for(int8_t y : trs) {
					for(int8_t z : trs) {
						int8_t t[3] = {x, y, z};
						p.setTrans(t);
						p.toEMsoft(em, true);
						q.fromEMsoft(em);
						if(p != q) {
							os << "EMsoft conversions not self consistent\n";
							return false;
						}
					}
				}
			}
		}

		//if we made it this far all tests passed
		os << "all tests passed\n";
		return true;
	}

	//@brief   : test all general position constructors (+ static factories)
	//@param os: output stream to write error messages to
	//@return  : true/false if the self tests pass/fail
	bool testBuild(std::ostream& os) {
		//make sure default constructor works and check parameterless static constructors
		{
			GenPos p;
			if(p != GenPos::Identity()) {
				os << "default general position isn't identity\n";
				return false;
			}

			int8_t mat[9] = {1, 0, 0,  0, 1, 0,  0, 0, 1};
			int8_t const * m = p.getMat3();
			if(!std::equal(m, m+9, mat)) {
				os << "identity isn't identity\n";
				return false;
			}

			p = GenPos::Inversion();
			m = p.getMat3();
			mat[0] = mat[4] = mat[8] = -1;
			if(!std::equal(m, m+9, mat)) {
				os << "inversion isn't inversion\n";
				return false;
			}
		}

		//make sure mirror + 2 fold constructors work
		{
			//build possible normals
			int8_t normals[11][3] = {
				{ 1, 0, 0}, { 0, 1, 0}, { 0, 0, 1}, { 0,-1, 1},
				{ 0, 1, 1}, { 1,-1, 0}, { 1, 1, 0}, {-1, 0, 1},
				{ 1, 0, 1}, { 1, 2, 0}, { 2, 1, 0},
			};

			//loop over normals constructing matrices
			for(size_t i = 0; i < 11; i++) {
				for(bool b : {true, false}) {
					//build mirror and make sure it is correct
					GenPos p = GenPos::Mirror(normals[i], b);
					if(-2 != p.order() || !std::equal(normals[i], normals[i]+3, p.axis())) {
						os << "Mirror() constructor wrong\n";
						return false;
					}

					//build 2 fold and make sure it is correct
					p = GenPos::Two(normals[i], b);
					if( 2 != p.order() || !std::equal(normals[i], normals[i]+3, p.axis())) {
						os << "Two() constructor wrong\n";
						return false;
					}
				}
			}
		}

		//make sure 3 fold constructor works
		{
			int8_t normals[4][3] = {
				{ 1, 1, 1},
				{-1, 1, 1},
				{ 1,-1, 1},
				{-1,-1, 1},
			};

			for(size_t i = 0; i < 4; i++) {
				//test 3 @ n
				GenPos p = GenPos::Three(normals[i], false);
				if( 3 != p.order() || !std::equal(normals[i], normals[i]+3, p.axis())) {
					os << "Three() constructor wrong\n";
					return false;
				}

				//test -3 @ n
				p = GenPos::Three(normals[i], true);
				if(-3 != p.order() || !std::equal(normals[i], normals[i]+3, p.axis())) {
					os << "Three() constructor wrong\n";
					return false;
				}
			}
		}

		//make sure z rotation constructor works
		{
			int8_t norm[3] = {0,0,1};
			for(int n : {2, 3, 4, 6, -2, -3, -4, -6}) {
				GenPos p = GenPos::Z(n);
				if(n != p.order() || !std::equal(norm, norm+3, p.axis())) {
					os << "Z(" << n << ") constructor wrong\n";
					return false;
				}
			}
			//handle +/-1 specially since they have degenerate rotation axis
			if(GenPos::Z(1) != GenPos::Identity() ||GenPos::Z(-1) != GenPos::Inversion()) {
				os << "Z(+/-1) constructor wrong\n";
				return false;
			} 
		}

		//if we made it this far all constructors work
		return true;

	}

	//@brief   : test all general position matrix operations
	//@param os: output stream to write error messages to
	//@return  : true/false if the self tests pass/fail
	bool testMat(std::ostream& os) {
		//test set / get
		{
			//3x3 matrix part
			GenPos p;
			for(uint32_t i = 0; i < 64; i++) {//loop over possible 3x3 matrices
				int8_t const * m = GenPos(i).getMat3();//get the 3x3 part
				p.setMat3(m);//tell p to use that matrix part
				int8_t const * n = GenPos(i).getMat3();//make the 3x3 part from p
				if(!std::equal(m, m+9, n)) {
					os << "set/getMat3() inconsistent\n";
					return false;
				}
			}

			//translation part
			for(int8_t x = 0; x < 24; x++) {
				for(int8_t y = 0; y < 24; y++) {
					for(int8_t z = 0; z < 24; z++) {
						//test fractional 24ths setting
						int8_t t[3] = {x, y, z};
						p.setTrans(t);
						if(!std::equal(t, t+3, p.getTrans())) {
							os << "set/getTrans() inconsistent\n";
							return false;
						}

						//test real setting
						double q[3];
						double r[3] = {double(x)/24, double(y)/24, double(z)/24};
						p.removeTrans();
						p.setTransReal(r);
						p.getTransReal(q);
						if(!std::equal(t, t+3, p.getTrans()) || !std::equal(q, q+3, r)) {
							os << "set/getTransReal() inconsistent\n";
							return false;
						}
					}
				}
			}
		}

		//check property queries and order
		for(uint32_t i = 0; i < 64; i++) {//loop over possible 3x3 matrices
			GenPos p(i);//build matrix
			GenPos q(p);//copy
			const int order = p.order();//get rotational order
			GenPos ident = GenPos::Identity();//most matrices^order should be identity
			if(order < 0 && ( (-order) % 2) == 1) ident = GenPos::Inversion();//odd roto inversions should be inversion
			for(int j = 1; j < std::abs(order); j++) q = q*p;//compute p^order
			if(q != ident ) {//make sure we get the expected result
				os << "mat^order != identity\n";
				return false;
			}

			if(order * p.det() < 0) {
				os << "sign(order) !== sign(det())\n";
				return false;
			}
		}

		//make sure that inverse works
		for(uint32_t i = 0; i < 64; i++) {//loop over possible 3x3 matrices
			int8_t t[3] = {3, 8, 18};//random translation
			GenPos p(i);//build matrix from 3x3 part
			p.setTrans(t);//set translation
			p *= p.inverse();
			if(GenPos::Identity() != p) {
				os << "p * p.inverse() != identity\n";
				return false;
			}
		}

		//make sure origin shift is equivalent to transforming with identity | t
		for(uint32_t i = 0; i < 64; i++) {//loop over possible 3x3 matrices
			//save translation + identity | t
			int8_t t[3] = {3, 8, 18};//random translation
			GenPos q = GenPos::Identity();
			q.setTrans(t);

			//check for equivalence
			GenPos p(i);
			int8_t t2[3] = {4, 12, 9};//different random translation
			p.setTrans(t2);
			if(p.shiftOrigin(t) != p.transform(q)) {
				os << "shiftOrigin is not special case of transform\n";
				return false;
			}
		}

		//make sure that transform is invertible
		for(uint32_t i = 1; i < 64; i++) {//loop over possible 3x3 matrices except for identity
			//build a random matrix
			int8_t t[3] = {3, 8, 18};//random translation
			GenPos p(i);//build matrix from 3x3 part
			p.setTrans(t);//set translation

			//compute like bounds (don't multiply cubic with hexagonal matrices)
			const uint32_t start = i < 48 ?  1 : 48;
			const uint32_t end   = i < 48 ? 48 : 64;
			for(uint32_t j = start; j < end; j++) {//loop over possible 3x3 matrices in the same group
				//build another random matrix
				int8_t t2[3] = {4, 12, 9};//different random translation
				GenPos q(j);
				q.setTrans(t2);

				//transform one by the other
				GenPos r = p.transform(q);
				if(p == r || p != r.transform(q.inverse())) {//make sure a transform was actually applied and that it was inverted
					os << "transform is not reversible\n";
					return false;
				}
			}
		}

		//if we made it this far all matrix operations work
		return true;
	}
}
