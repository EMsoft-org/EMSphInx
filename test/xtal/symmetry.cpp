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
//      test program for functions in include/xtal/symmetry.hpp       //
////////////////////////////////////////////////////////////////////////

#include <iostream>

namespace xtal {

	//@brief   : run all symmetry unit tests
	//@param os: output stream to write error messages to
	//@return  : true/false if the self tests pass/fail
	template <typename Real> bool runTests(std::ostream& os);

	//@brief   : test point group constructors
	//@param os: output stream to write error messages to
	//@return  : true/false if the tests pass/fail
	bool testBuild(std::ostream& os);

	//@brief   : test point group symmetry against space groups
	//@param os: output stream to write error messages to
	//@return  : true/false if the tests pass/fail
	template <typename Real> bool testSgSym(std::ostream& os);

	//@brief   : test point group relationships
	//@param os: output stream to write error messages to
	//@return  : true/false if the tests pass/fail
	template <typename Real> bool testRel(std::ostream& os);

	//@brief   : test point group conversions
	//@param os: output stream to write error messages to
	//@return  : true/false if the tests pass/fail
	bool testConv(std::ostream& os);

	//@brief   : test symmetry operations of rotational point groups
	//@param os: output stream to write error messages to
	//@return  : true/false if the tests pass/fail
	template <typename Real> bool testRot(std::ostream& os);
}

int main() {
	//select output stream
	std::ostream& os = std::cout;

	// //run unit tests
	const bool passed =  xtal::runTests<float >(os)
	                  && xtal::runTests<double>(os);

	//return result
	return passed ? EXIT_SUCCESS : EXIT_FAILURE;
}

#include <iomanip>
#include <limits>
#include <set>
#include <vector>
#include <string>
#include <random>

#include "xtal/symmetry.hpp"
#include "xtal/position.hpp"
#include "xtal/rotations.hpp"
#include "xtal/hm.hpp"

//@brief  : square lambert projection from unit square to unit hemisphere
//@param X: x coordinate in unit square (0,1)
//@param Y: y coordinate in unit square (0,1)
//@param x: location to write x coordinate on unit sphere
//@param y: location to write y coordinate on unit sphere
//@param z: location to write z coordinate on unit sphere
template <typename Real>
void squareToSphere(Real const& X, Real const& Y, Real& x, Real& y, Real& z) {
	static const Real kPi_4 = Real(0.7853981633974483096156608458199);//pi/4
	const Real sX = Real(2) * X - 1;//[0,1] -> [-1, 1]
	const Real sY = Real(2) * Y - 1;//[0,1] -> [-1, 1]
	const Real aX = std::abs(sX);
	const Real aY = std::abs(sY);
	const Real vMax = std::max<Real>(aX, aY);
	if(vMax <= std::numeric_limits<Real>::epsilon()) {
		x = y = Real(0);
		z = Real(1);
	} else {
		if(vMax > Real(1) + std::numeric_limits<Real>::epsilon()) throw std::runtime_error("point doesn't lie in square (0,0) -> (1,1)");
		if(aX <= aY) {
			const Real q  = sY * std::sqrt(Real(2) - sY * sY);
			const Real qq = kPi_4 * sX / sY;
			x = q * std::sin(qq);
			y = q * std::cos(qq);
		} else {
			const Real q = sX * std::sqrt(Real(2) - sX * sX);
			const Real qq = kPi_4 * sY / sX;
			x = q * std::cos(qq);
			y = q * std::sin(qq);
		}
		z = Real(1) - vMax * vMax;
		const Real mag = std::sqrt(x*x + y*y + z*z);
		x /= mag; y /= mag; z /= mag;
	}
}

namespace xtal {
	//enumerate list of point group names
	const std::vector<std::string> names = {
		    "1",    "-1",   "121",   "112",   "1m1",   "11m", "12/m1", "112/m",
		  "222",   "mm2",   "mmm",     "4",  "-4"  ,   "4/m",   "422",   "4mm",
		 "-42m",  "-4m2", "4/mmm",     "3",    "-3",   "321",   "312",   "3m1",
		  "31m",  "-3m1",  "-31m",     "6",    "-6",   "6/m",   "622",   "6mm",
		 "-62m",  "-6m2", "6/mmm",    "23",    "m3",   "432",  "-43m",   "m3m",
	};

	//@brief : get a list of all point groups
	//@return: list of groups
	std::vector<PointGroup> getGroups() {
		//build list of point groups
		std::vector<PointGroup> groups;
		for(const std::string n : names) groups.push_back(PointGroup(n));
		groups.push_back(PointGroup::BuildOrtho45("222r"));
		groups.push_back(PointGroup::BuildOrtho45("mm2r"));
		groups.push_back(PointGroup::BuildOrtho45("mmmr"));
		return groups;
	}

	//@brief   : test point group constructors
	//@param os: output stream to write error messages to
	//@return  : true/false if the tests pass/fail
	bool testBuild(std::ostream& os) {
		os << "testing space group number construction\n";
		for(uint_fast8_t i = 1; i <= 230; i++) {//loop over space groups
			//get generators and split into translation / no translation
			std::vector<GenPos> gen = GenPos::CloseSet(HermannMaguin(i).generators(NULL, true));//get generators
			std::set<GenPos> noTrn;
			std::vector<GenPos> cen;
			for(GenPos& p : gen) {
				GenPos n = p;//copy position
				n.removeTrans();//remove translation component
				noTrn.insert(n);//save translation free position
				if(p.hasTranslation() && GenPos::Identity() == n) cen.push_back(p);//save any identity translations
			}

			//get the point group for this number and then the primitive space group
			PointGroup pg(i);
			size_t sg = pg.symmorphic();
			if(0 == sg) {
				os << "inconsistent symmorphic() for PointGroup(" << i << ")\n";
				return false;
			}

			//make sure primitive group matches symmorphic(P)
			std::vector<GenPos> genP = GenPos::CloseSet(HermannMaguin(sg).generators(NULL, true));//get generators
			std::vector<GenPos> genRed(noTrn.begin(), noTrn.end());//convert translation reduced operators to vector
			if(genP != genRed) {
				os << "inconsistent symmorphic() generators for PointGroup(" << i << ")\n";
				return false;
			}
		}

		//test by name construction
		os << "testing name construction\n";
		for(const std::string n : names) {
			if(PointGroup(n).name() != n) {
				os << "PointGroup(str).name() != str for str == " << n << '\n';
				return false; 
			}
		}

		//test special orthorhombic constructor
		os << "testing special construction\n";
		for(std::string n : {"222r", "mm2r", "mmmr"}) {
			if(PointGroup::BuildOrtho45(n).name() != n) {
				os << "BuildOrtho45(str).name() != str for str == " << n << '\n';
				return false; 
			}
		}

		//if we made it this far everything passed
		return true;
	}

	//@brief   : test point symmetry against space groups
	//@param os: output stream to write error messages to
	//@return  : true/false if the tests pass/fail
	template <typename Real> bool testSgSym(std::ostream& os) {
		os << "testing that symmorphic() is correct\n";
		std::vector<PointGroup> groups = getGroups();
		for(const PointGroup& pg : groups) {//loop over point groups
			for(char c : {'P', 'A', 'B', 'C', 'F', 'I', 'R'}) {//loop over lattice centering types
				uint_fast8_t sg = pg.symmorphic(c);//get symmorphic group for centering
				if(sg != 0 && NULL == pg.symmorphicTrns<Real>()) {//there is a corresponding space group
					HermannMaguin hm(sg);//get the space group
					std::string name = hm.shortSym().to_string();//get the space group name
					if(name.front() != c) {//make sure centering was correct
						os << "symmorphic space group for point group " << pg.name() << " with " << c << " centering (" << sg << "): " << name << '\n';
						return false;
					}
				}
			}
		}

		os << "testing symmetry operations against primitive space group general positions\n";
		std::set<PointGroup> rotSg, laueSg;
		for(const PointGroup& pg : groups) {//loop over point groups
			//check if we have a matching primitive group to work with
			uint_fast8_t sg = pg.symmorphic('P');
			if(0 == sg) {
				os << "no primitive space group for " << pg.name() << '\n';
				return false;
			}

			//change cell if needed
			HermannMaguin hm(sg);
			bool hasStd = true;//has a corresponding ITA space group
			std::string altMono = "b";
			if(NULL != pg.symmorphicTrns<Real>()) {//manually handle alternate monoclinic setting
				if("112" == pg.name() || "11m" == pg.name() || "112/m" == pg.name()) hm = hm.changeMonoCell(1, "c");
				else hasStd = false;//222r type
			}

			//make sure we have the right space group
			std::string name = hm.shortSym().to_string();//get the space group name
			if(!hasStd) name+= 'r';//handle 222r types
			name.erase(std::remove(name.begin(), name.end(), ' '), name.end());//remove white space
			if('P' != name[0] || 0 != name.compare(1, std::string::npos, pg.name())) {
				os << "primitive space group " << name << " doesn't match point group " << pg.name() << '\n';
				return false;
			}

			//space group matches, get the general positions
			std::vector<GenPos> gen = GenPos::CloseSet(hm.generators());//we want hexagonal instead of rhombohedral generators where possible
			if(gen.size() != pg.order()) {
				os << pg.name() << " order doesn't match space group\n";
				return false;
			}

			//check inversion flag
			const bool sgInv = gen.end() != std::find(gen.begin(), gen.end(), GenPos::Inversion());
			if(sgInv != pg.inversion()) {
				os << pg.name() << " inversion flag doesn't match space group\n";
				return false;
			}

			//loop over general positions categorizing
			bool hasMirZ = false, hasInv  = false;
			bool hasMirE = false, hasMirY = false;
			int maxZ = 1;
			int minZ = 1;
			std::vector<GenPos> rots, mirs;
			const int8_t eigVec[15][3] = {//these are all the possible eigenvectors for the 64 crystallographic 3x3 general position matrices
				{ 1, 0, 0}, { 0, 1, 0}, { 0, 0, 1}, { 0, 1, 1}, { 1, 0, 1},
				{ 1, 1, 0}, { 0,-1, 1}, {-1, 0, 1}, { 1,-1, 0}, { 1, 1, 1},
				{ 1,-1, 1}, {-1,-1, 1}, {-1, 1, 1}, { 2, 1, 0}, { 1, 2, 0},
			};
			int maxRot[15] = {1};//higher order rotation (or rotoinversion) detected at each axis
			for(GenPos& p : gen) {//loop over general positions
				//put into correct bin
				const int order = p.order();//get rotational order
				if     (order >   0) rots.push_back(p);//rotation (or identity)
				else if(order == -1) hasInv = true;
				else if(order == -2) mirs.push_back(p);
				else if(order <  -2) {}//skip rotoinverters
				else {
					os << "unclassifiable general position in space group for " << pg.name() << '\n';
					return false;
				}

				//determine rotation axis
				size_t idx;
				int8_t const * ax = p.axis();
				for(idx = 0; idx < 15; idx++) {
					if(std::equal(ax, ax + 3, eigVec[idx])) break;
				}
				if(15 == idx) {
					os << "unclassifiable general position axis in space group for " << pg.name() << '\n';
					return false;
				}

				//accumulate maximum order in each direction
				if(std::abs(order) > std::abs(maxRot[idx])) {//new largest magnitude
					maxRot[idx] = order;
				} else if(std::abs(order) == std::abs(maxRot[idx])) {//equal to largest magnitude
					if(order % 2 == order > 0 ? 0 : 1) maxRot[idx] = order;//prefer proper inversions for even orders and rotoinversinos for odd orders
				}

				//check for z elements specially
				if(std::equal(eigVec[2], eigVec[2]+3, ax)) {//this is an operation about the z axis
					minZ = std::min(order, minZ);//save largest z order
					maxZ = std::max(order, maxZ);//save largest z order
					if(-2 == p.order()) hasMirZ = true;//save mirror
				}

				//check for x/y mirror specially
				if(-2 == p.order()) {//mirror
					if(0 == ax[2]) {//z == 0
						hasMirE = true;//there is at least 1 equatorial mirror
					}
				}
			}

			//save purely rotational groups and laue groups
			if(rots.size() == gen.size()) {
				size_t rCount = rotSg.size();
				rotSg.insert(pg);
				if(rCount == rotSg.size()) {
					os << "duplicate purely rotational group\n";
					return false;
				}
			} else if(hasInv) {
				size_t lCount = laueSg.size();
				laueSg.insert(pg);
				if(lCount == laueSg.size()) {
					os << "duplicate laue group\n";
					return false;
				}
			}

			//now that we've parsed out the general positions, check symmetry attributes

			//check inversion
			if(hasInv != pg.inversion()) {
				os << pg.name() << " inversion flag doesn't match space group\n";
				return false;
			}

			//check z mirror
			if(hasMirZ != pg.zMirror()) {
				os << pg.name() << " z mirror flag doesn't match space group\n";
				return false;
			}

			//check zRot
			if(maxZ != pg.zRot()) {
				os << pg.name() << " z rotational order doesn't match space group\n";
				return false;
			}

			//check operation counts
			if(mirs.size() != pg.numMirror()) {
				os << pg.name() << " mirror count doesn't match space group\n";
				return false;
			}
			if(rots.size() != pg.numRotOps()) {
				os << pg.name() << " rotation count doesn't match space group\n";
				return false;
			}

			//check that enantiomorphism() is consistent with inversion() / numMirror()
			if(pg.enantiomorphism() != ( mirs.empty() && !pg.inversion() && -4 != minZ ) ) {
				os << pg.name() << " enantiomorphism inconsistent with inversion() / numMirror()\n";
				os << mirs.size() << ' ' << pg.inversion() << '\n';
				os << mirs.empty() << ' ' << !pg.inversion() << '\n';
				os << pg.enantiomorphism() << '\n';
				return false;
			}

			//make sure rotational components are a closed set
			if(GenPos::CloseSet(rots).size() != rots.size()) {
				os << pg.name() << " rotational elements aren't a closed set\n";
				return false;
			}

			//convert rotations to quaternions and check that they are the same
			Real om[9];
			Quat<Real> qu;
			Quat<Real> const * ops = (Quat<Real> const*) pg.rotOps<Real>();
			const bool hex = 0 == maxZ % 3;
			const Real minDot = Real(1) - std::numeric_limits<Real>::epsilon() * 10;//threshold for 2 quaternions to be the same
			for(GenPos& p : rots) {
				//get 3x3 matrix as real type
				if(hex) p.getMat3HexCart(om);
				else std::copy(p.getMat3(), p.getMat3()+9, om);
				
				//transform if needed and convert to quaternion
				if(!hasStd) {
					Real const* a = pg.symmorphicTrns<Real>();
					const Real atOm[9] = {//a^T * om
						a[0] *   om[0] + a[3] *   om[3] + a[6] *   om[6],   a[0] *   om[1] + a[3] *   om[4] + a[6] *   om[7],   a[0] *   om[2] + a[3] *   om[5] + a[6] *   om[8],
						a[1] *   om[0] + a[4] *   om[3] + a[7] *   om[6],   a[1] *   om[1] + a[4] *   om[4] + a[7] *   om[7],   a[1] *   om[2] + a[4] *   om[5] + a[7] *   om[8],
						a[2] *   om[0] + a[5] *   om[3] + a[8] *   om[6],   a[2] *   om[1] + a[5] *   om[4] + a[8] *   om[7],   a[2] *   om[2] + a[5] *   om[5] + a[8] *   om[8],
					};
					const Real atOmA[9] = {//a^T * om * a
						atOm[0] * a[0] + atOm[1] * a[3] + atOm[2] * a[6],   atOm[0] * a[1] + atOm[1] * a[4] + atOm[2] * a[7],   atOm[0] * a[2] + atOm[1] * a[5] + atOm[2] * a[8],
						atOm[3] * a[0] + atOm[4] * a[3] + atOm[5] * a[6],   atOm[3] * a[1] + atOm[4] * a[4] + atOm[5] * a[7],   atOm[3] * a[2] + atOm[4] * a[5] + atOm[5] * a[8],
						atOm[6] * a[0] + atOm[7] * a[3] + atOm[8] * a[6],   atOm[6] * a[1] + atOm[7] * a[4] + atOm[8] * a[7],   atOm[6] * a[2] + atOm[7] * a[5] + atOm[8] * a[8],
					};
					std::copy(atOmA, atOmA + 9, om);
				}
				om2qu(om, qu.data());

				//loop over ops finding the best match
				Real maxDot = 0;
				for(size_t i = 0; i < rots.size(); i++) maxDot = std::max(maxDot, std::fabs(qu.dot(ops[i])));
				if(maxDot <= minDot) {
					os << "failed to find all space group rotations in rotOps() for " << pg.name() << '\n';
					return false;
				}
			}

			//check that mirror planes are the same
			Real const * norms = pg.mirrors<Real>();
			for(GenPos& p : mirs) {
				//get mirror plane normal as unit vector
				int8_t const * ax = p.axis();//get mirror plane normal
				Real n[3] = {Real(ax[0]), Real(ax[1]), Real(ax[2])};//convert to real
				if(hex) {//correct for hex frame if needed
					n[0]  = n[0] - n[1] / 2;
					n[1] *= std::sqrt(Real(3)) / 2;
				}
				Real mag = std::sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);//get magnitude
				n[0] /= mag; n[1] /= mag; n[2] /= mag;//normalize plane normal

				//transform if needed
				if(!hasStd) {
					Real const* a = pg.symmorphicTrns<Real>();
					const Real np[3] = {//a * n
						a[0] * n[0] + a[1] * n[1] + a[2] * n[2],
						a[3] * n[0] + a[4] * n[1] + a[5] * n[2],
						a[6] * n[0] + a[7] * n[1] + a[8] * n[2],
					};
					std::copy(np, np + 3, n);
				}
				if(std::fabs(n[1]) >= 0.99) hasMirY = true;

				//loop over ops finding the best match
				Real maxDot = 0;
				for(size_t i = 0; i < mirs.size(); i++) maxDot = std::max(maxDot, std::fabs(n[0]*norms[3*i+0] + n[1]*norms[3*i+1] + n[2]*norms[3*i+2]));
				if(maxDot <= minDot) {
					os << "failed to find all space group mirrors in rotOps() for " << pg.name() << '\n';
					return false;
				}
			}

			//check mmType
			if(hasMirE) {
				if(0 == pg.mmType()) {//z rotation + equatorial mirrors should have >0 mm type
					os << pg.name() << " mmType doesn't match space group (is 0)\n";
					return false;
				} else {
					if(1 == pg.mmType() && hasMirY) {
						//unrotated mm type
					} else if(2 == pg.mmType() && !hasMirY) {
						//rotated mm type
					} else {
						os << pg.name() << " mmType doesn't match space group (mismatch)\n";
						return false;
					}
				}
			} else {
				if(0 != pg.mmType()) {//no equatorial mirrors should have no mm type
					os << pg.name() << " mmType doesn't match space group (not 0)\n";
					return false;
				}
			}

			//convert principal directions to real unit vectors
			Real eigVecReal[15][3] = {//these are all the possible eigenvectors for the 64 crystallographic 3x3 general position matrices
				{ 1, 0, 0}, { 0, 1, 0}, { 0, 0, 1}, { 0, 1, 1}, { 1, 0, 1},
				{ 1, 1, 0}, { 0,-1, 1}, {-1, 0, 1}, { 1,-1, 0}, { 1, 1, 1},
				{ 1,-1, 1}, {-1,-1, 1}, {-1, 1, 1}, { 2, 1, 0}, { 1, 2, 0},
			};
			for(size_t i = 0; i < 15; i++) {
				if(hex) {//correct for hex frame if needed
					eigVecReal[i][0]  = eigVecReal[i][0] - eigVecReal[i][1] / 2;
					eigVecReal[i][1] *= std::sqrt(Real(3)) / 2;
				}
				Real mag = std::sqrt(eigVecReal[i][0]*eigVecReal[i][0] +
				                     eigVecReal[i][1]*eigVecReal[i][1] + 
				                     eigVecReal[i][2]*eigVecReal[i][2]);//get magnitude
				eigVecReal[i][0] /= mag; eigVecReal[i][1] /= mag; eigVecReal[i][2] /= mag;//normalize plane normal

				//transform if needed
				if(!hasStd) {
					Real const* a = pg.symmorphicTrns<Real>();
					const Real np[3] = {//a * n
						a[0] * eigVecReal[i][0] + a[1] * eigVecReal[i][1] + a[2] * eigVecReal[i][2],
						a[3] * eigVecReal[i][0] + a[4] * eigVecReal[i][1] + a[5] * eigVecReal[i][2],
						a[6] * eigVecReal[i][0] + a[7] * eigVecReal[i][1] + a[8] * eigVecReal[i][2],
					};
					std::copy(np, np + 3, eigVecReal[i]);
				}
			}

			//check that principal rotational symmetry operators are correct
			size_t numOps = 0;
			Real const * ax = pg.rotAxis<Real>();//get rotation axis
			for(size_t i = 0; i < 15; i++) {
				if( maxRot[i] > 1 || maxRot[i] < -2) {//there is a rotation to find
					++numOps;
					Real maxDot = 0;
					for(size_t j = 0; j < pg.numRotAxis(); j++) {//loop over rotation axis of point group
						if(Real(maxRot[i]) != ax[4*j]) continue;//rotational order doesn't match
						const Real dot = eigVecReal[i][0] * ax[4*j+1]
						               + eigVecReal[i][1] * ax[4*j+2]
						               + eigVecReal[i][2] * ax[4*j+3];
						maxDot = std::max( maxDot, std::fabs(dot) );
					}
					if(maxDot < minDot) {
						os << "rotAxis() doesn't match detected rotations for " << pg.name() << '\n';
						return false;
					}
				} 
			}

			//make sure the number of operators matches
			if(numOps != pg.numRotAxis()) {
				os << "numRotAxis() doesn't match detected rotations for " << pg.name() << '\n';
				return false;
			}

			//loop over generators convering to 3x3 real matrix
			std::vector< std::vector<Real> > realPos;//this is supper wasteful but thats ok
			for(GenPos& p : gen) {
				//get 3x3 matrix as real type
				std::vector<Real> om(9);
				if(hex) p.getMat3HexCart(om.data());
				else std::copy(p.getMat3(), p.getMat3()+9, om.data());
				
				//transform if needed and convert to quaternion
				if(!hasStd) {
					Real const* a = pg.symmorphicTrns<Real>();
					const Real atOm[9] = {//a^T * om
						a[0] *   om[0] + a[3] *   om[3] + a[6] *   om[6],   a[0] *   om[1] + a[3] *   om[4] + a[6] *   om[7],   a[0] *   om[2] + a[3] *   om[5] + a[6] *   om[8],
						a[1] *   om[0] + a[4] *   om[3] + a[7] *   om[6],   a[1] *   om[1] + a[4] *   om[4] + a[7] *   om[7],   a[1] *   om[2] + a[4] *   om[5] + a[7] *   om[8],
						a[2] *   om[0] + a[5] *   om[3] + a[8] *   om[6],   a[2] *   om[1] + a[5] *   om[4] + a[8] *   om[7],   a[2] *   om[2] + a[5] *   om[5] + a[8] *   om[8],
					};
					const Real atOmA[9] = {//a^T * om * a
						atOm[0] * a[0] + atOm[1] * a[3] + atOm[2] * a[6],   atOm[0] * a[1] + atOm[1] * a[4] + atOm[2] * a[7],   atOm[0] * a[2] + atOm[1] * a[5] + atOm[2] * a[8],
						atOm[3] * a[0] + atOm[4] * a[3] + atOm[5] * a[6],   atOm[3] * a[1] + atOm[4] * a[4] + atOm[5] * a[7],   atOm[3] * a[2] + atOm[4] * a[5] + atOm[5] * a[8],
						atOm[6] * a[0] + atOm[7] * a[3] + atOm[8] * a[6],   atOm[6] * a[1] + atOm[7] * a[4] + atOm[8] * a[7],   atOm[6] * a[2] + atOm[7] * a[5] + atOm[8] * a[8],
					};
					std::copy(atOmA, atOmA + 9, om.data());
				}
				realPos.push_back(om);
			}

			//loop over a grid of all unit directions doing some tests
			std::vector<Real> lin(501);
			for(size_t i = 0; i < lin.size(); i++) lin[i] = Real(i) / (lin.size() - 1);// [0,1]
			size_t fsCount = 0;
			for(Real& x : lin) {
				for(Real& y : lin) {
					//get the unit direction in the north hemisphere
					Real nh[3];
					squareToSphere(x, y, nh[0], nh[1], nh[2]);
					Real sh[3] = {nh[0], nh[1], -nh[2]};

					//get fundamental sector direction and check if already inide
					Real nhFs[3], shFs[3];
					const bool nhIn = pg.fsDir(nh, nhFs);
					const bool shIn = pg.fsDir(sh, shFs);
					if(nhIn) ++fsCount;
					if(shIn) ++fsCount;

					//make sure reduction was a symmetric equivalent
					Real maxDotNh = 0, maxDotSh = 0;
					for(std::vector<Real>& om : realPos) {//loop over general positions
						//compute om * nh
						Real omNh[3] = {
							om[0] * nh[0] + om[1] * nh[1] + om[2] * nh[2],
							om[3] * nh[0] + om[4] * nh[1] + om[5] * nh[2],
							om[6] * nh[0] + om[7] * nh[1] + om[8] * nh[2],
						};

						//compute om * sh
						Real omSh[3] = {
							om[0] * sh[0] + om[1] * sh[1] + om[2] * sh[2],
							om[3] * sh[0] + om[4] * sh[1] + om[5] * sh[2],
							om[6] * sh[0] + om[7] * sh[1] + om[8] * sh[2],
						};

						//check for best match
						maxDotNh = std::max(maxDotNh, omNh[0] * nhFs[0] + omNh[1] * nhFs[1] + omNh[2] * nhFs[2]);
						maxDotSh = std::max(maxDotSh, omSh[0] * shFs[0] + omSh[1] * shFs[1] + omSh[2] * shFs[2]);
					}
					if(maxDotNh < minDot || maxDotSh < minDot) {
						os << "fundamental sector direction isn't genpos * n for " << pg.name() << '\n';
						return false;
					}
				}
			}

			//make sure fraction of directions in fundamental sector matches order
			const Real fsFrac = Real(lin.size() * lin.size() * 2) / fsCount;
			if(std::round(fsFrac) != pg.order()) {
				os << "fundmental sector fraction doesn't match order for " << pg.name() << '\n';
				os << fsFrac << ' ' << (int) pg.order() << '\n';
				return false;
			}
		}
		
		//this is a cludge until alternate settings is wrapped up
		rotSg .insert(PointGroup("112"));
		rotSg .insert(PointGroup::BuildOrtho45("222r"));
		laueSg.insert(PointGroup("112/m"));
		laueSg.insert(PointGroup::BuildOrtho45("mmmr"));

		//make sure rotational point group is actually purely rotational
		for(const PointGroup& pg : groups) {//loop over point groups
			if(pg.rotationGroup().numMirror() > 0 || pg.rotationGroup().inversion()) {
				os << "rotation group of " << pg.name() << " isn't a purely rotation group\n";
				os << pg.rotationGroup().name() << '\n';
				return false;
			}
			if(1 != rotSg.count(pg.rotationGroup())) {
				os << "rotation group of " << pg.name() << " not accumlated from space groups\n";
				os << pg.rotationGroup().name() << '\n';
				return false;
			}
		}
		
		//make sure the laue group is actualy inversion symmetric
		for(const PointGroup& pg : groups) {//loop over point groups
			if(!pg.laueGroup().inversion()) {
				os << "laue group of " << pg.name() << " isn't a centrosymmetric\n";
				os << pg.laueGroup().name() << '\n';
				return false;
			}
			if(1 != laueSg.count(pg.laueGroup())) {
				os << "laue group of " << pg.name() << " not accumlated from space groups\n";
				os << pg.laueGroup().name() << '\n';
				return false;
			}
		}

		//if we made it this far everything passed
		return true;
	}

	//@brief   : test point group relationships
	//@param os: output stream to write error messages to
	//@return  : true/false if the tests pass/fail
	template <typename Real> bool testRel(std::ostream& os) {
		os << "testing point group relations\n";
		std::vector<PointGroup> groups = getGroups();

		//make sure rotation point group is a subset of group rotations
		for(const PointGroup& pg : groups) {//loop over point groups
			//make sure number of operators is the same
			size_t num  = pg.numRotOps();
			if(pg.rotationGroup().numRotOps() != num) {
				os << "rotation group of " << pg.name() << " has different rotation count\n";
				return false;
			}

			//make sure actual operators are the same
			Real const * rots  = pg.rotOps<Real>();
			Real const * rotsR = pg.rotationGroup().rotOps<Real>();
			if(!std::equal(rots, rots + num * 4, rotsR)) {
				os << "rotation group of " << pg.name() << " has different rotations\n";
				return false;
			}
		}

		//make sure group is a subset of the laue group
		for(const PointGroup& pg : groups) {//loop over point groups
			//sanity check number of operators
			size_t numR  = pg.numRotOps();
			size_t numM  = pg.numMirror();
			size_t numRL = pg.laueGroup().numRotOps();
			size_t numML = pg.laueGroup().numMirror();
			if(numRL < numR || numML < numM) {
				os << "laue group of " << pg.name() << " has fewer operators\n";
				return false;
			}

			//make sure mirror planes are subset of laue group mirrors
			Real const * mirs  = pg.mirrors<Real>();
			Real const * mirsL = pg.laueGroup().mirrors<Real>();
			for(size_t i = 0; i < numM; i++) {//loop over mirrors
				bool found = false;
				for(size_t j = 0; j < numML; j++) {//loop over laue group mirrors looking for a match
					if(std::equal(mirs + 3 * i, mirs + 3 * i + 3, mirsL + 3 * j)) {//we found a match
						found = true;//flag
						break;//stop looking
					}
				}
				if(!found) {//check that the mirror was found
					os << "laue group of " << pg.name() << " doesn't have all mirror planes\n";
					return false;
				}
			}

			//do the same for rotations
			Real const * rots  = pg.rotOps<Real>();
			Real const * rotsL = pg.laueGroup().rotOps<Real>();
			for(size_t i = 0; i < numR; i++) {//loop over rotations
				bool found = false;
				for(size_t j = 0; j < numRL; j++) {//loop over laue group rotations looking for a match
					if(std::equal(rots + 4 * i, rots + 4 * i + 4, rotsL + 4 * j)) {//we found a match
						found = true;//flag
						break;//stop looking
					}
				}
				if(!found) {//check that the mirror was found
					os << "laue group of " << pg.name() << " doesn't have all rotations\n";
					return false;
				}
			}
		}

		//check that laueName() is the same as laueGroup().name()
		for(const PointGroup& pg : groups) {//loop over point groups
			if(pg.laueName() != pg.laueGroup().laueName()) {
				os << "laueName() of " << pg.name() << " != laueGroup().laueName()\n";
				return false;
			}
		}

		//if we made it this far everything passed
		return true;
	}

	//@brief   : test point group conversions
	//@param os: output stream to write error messages to
	//@return  : true/false if the tests pass/fail
	bool testConv(std::ostream& os) {
		os << "testing conversions\n";
		std::vector<PointGroup> groups = {
			PointGroup("-1"   ),
			PointGroup("12/m1"),
			PointGroup("mmm"  ),
			PointGroup("4/m"  ),
			PointGroup("4/mmm"),
			PointGroup("-3"   ),
			PointGroup("-3m1" ),
			PointGroup("6/m"  ),
			PointGroup("6/mmm"),
			PointGroup("m3"   ),
			PointGroup("m3m"  ),
		};

		for(const PointGroup& pg : groups) {//loop over point groups
			if(pg.laueGroup() != PointGroup::FromTSL(pg.tslNum())) {
				os << pg.name() << " round trip TSL conversion isn't self consistent\n";
				return false;
			}

			if(pg.laueGroup() != PointGroup::FromHKL(pg.hklNum())) {
				os << pg.name() << " round trip HKL conversion isn't self consistent\n";
				return false;
			}
		}

		//if we made it this far everything passed
		return true;
	}


	//@brief   : test symmetry operations of rotational point groups
	//@param os: output stream to write error messages to
	//@return  : true/false if the tests pass/fail
	template <typename Real> bool testRot(std::ostream& os) {
		//get the purely rotational groups
		os << "testing rotataional group symmetry operations:\n";
		std::vector<PointGroup> groups = getGroups();
		{
			std::set<PointGroup> rotGroups;
			for(PointGroup& pg : groups) rotGroups.insert(pg.rotationGroup());
			groups.clear();
			groups.assign(rotGroups.begin(), rotGroups.end());
		}

		//grid the side length of the cubochoric cube
		std::vector<Real> lin(151);//this needs to be pretty fine for FZ fraction to work out
		for(size_t i = 0; i < lin.size(); i++) {
			lin[i] = Real(i) / (lin.size() - 1) - Real(0.5);// [-0.5,0.5]
			lin[i] *= Constants<Real>::cuA;//grid a side of the cubochoric cube
		}

		//loop over a grid of all orientations doing some tests
		std::vector<size_t> inFz(groups.size(), 0);
		const Real minDot = Real(1) - std::numeric_limits<Real>::epsilon() * 100;//how close do quat dots need to be for equality
		for(Real& x : lin) {
			os << "\t" << x << "/" << lin.back() << "    \r";
			os.flush();
			for(Real& y : lin) {
				for(Real& z : lin) {
					//get cubochoric grid point and convert to ro
					Real cu[3] = {x, y, z};
					Real ro[4];
					cu2ro(cu, ro);

					//count number of grid point in FZ
					for(size_t i = 0; i < groups.size(); i++) {
						//check if in fundamental zone
						const bool inFzRo = groups[i].roInFz(ro);
						if(inFzRo) ++inFz[i];

						//get FZ quat
						Real qu[4], quFz[4];
						ro2qu(ro, qu);
						groups[i].fzQu(qu, quFz);

						//check that fzQu is consistent with roInFz
						const Real dot = std::fabs(quat::dot(qu, quFz));
						const bool inFzQu = dot >= minDot;
						if(inFzQu != inFzRo) {
							os << "fzQu not consistent with roInFz for " << groups[i].name() << '\n';
							return false;
						}

						//make sure we actually got an fz orientations
						Real roFz[4];
						qu2ro(quFz, roFz);//convert to FZ
						roFz[3] *= minDot;//allow for roundoff error (pull everything towards the origin a bit)
						if(!groups[i].roInFz(roFz)) {
							os << "fzQu didn't return a fundamental zone quaternion\n";
							return false;
						}
					}
				}
			}
		}
		os << '\n';

		//make sure that the order is consistent with the FZ fraction
		for(size_t i = 0; i < groups.size(); i++) {
			const Real fzOrder = Real(lin.size() * lin.size() * lin.size()) / inFz[i];
			if(groups[i].order() != (int) std::round(fzOrder)) {
				os << groups[i].name() << " order() != FZ fraction " << fzOrder << '\n';
				return false;
			}
		}

		//test disoQu and nearbyQu
		Quat<Real> qu1, qu2, delta;
		std::mt19937_64 gen(0);
		const Real minW = std::cos(Constants<Real>::pi / 12);//minimum W for hexagonal disorientations = (1+sqrt(3)) / (2*sqrt(2))
		std::uniform_real_distribution<Real>  distR(-1, 1);
		std::uniform_int_distribution<size_t> distI(0, 23);
		std::vector< Quat<Real> const * > ops;//get operators once
		std::vector< size_t > numOps;//get # operators once
		for(PointGroup& pg : groups) {
			ops.push_back((Quat<Real> const *) pg.rotOps<Real>());
			numOps.push_back(pg.numRotOps());
		}
		for(size_t i = 0; i < 10000; i++) {
			//generate a random quaternion
			qu1.w = distR(gen);
			qu1.x = distR(gen);
			qu1.y = distR(gen);
			qu1.z = distR(gen);
			qu1 = qu1.normalize();
			qu1 = qu1.expl();

			//generate another random quaternion with angle < 15 degrees (guaranteed to be a disorientation)
			do {
				delta.w = distR(gen);
				if(std::fabs(delta.w) < minW) continue;// |w| can only get smaller when we normalize
				delta.x = distR(gen);
				delta.y = distR(gen);
				delta.z = distR(gen);
				delta = delta.normalize();
			} while( std::fabs(delta.w) <= minW );
			delta = delta.expl();
			qu2 = delta.conj() * qu1;//perturb first quaternion

			//now loop over point groups testing disoQu and nearbyQu
			for(size_t j = 0; j < groups.size(); j++) {
				//rotate qu1 and qu2 by random symmetry operators
				Quat<Real> qu1P = ops[j][distI(gen) % numOps[j]] * qu1;
				Quat<Real> qu2P = ops[j][distI(gen) % numOps[j]] * qu2;

				//compute disorientation and make sure it is correct
				Quat<Real> diso;
				groups[j].disoQu(qu1P.data(), qu2P.data(), diso.data());
				if(std::fabs(diso.w - delta.w) > 0.001) {
				// if(std::fabs(diso.dot(delta)) < minDot) {//don't check dot product since axis may be inconsistent for some point groups
					os << "disoQu failed to produce applied w for " << groups[j].name() << '\n';
					return false;
				}

				//compute nearby qu and make sure it is correct
				Quat<Real> near;
				groups[j].nearbyQu(qu1.data(), qu2P.data(), near.data());
				if(std::fabs(near.dot(qu2)) < minDot) {//don't check dot product since axis may be inconsistent for some point groups
					os << "nearbyQu failed to produce applied w for " << groups[j].name() << '\n';
					return false;
				}
			}
		}

		//if we made it this far everything passed
		return true;
	}

	//@brief   : run all symmetry unit tests
	//@param os: output stream to write error messages to
	//@return  : true/false if the self tests pass/fail
	template <typename Real> bool runTests(std::ostream& os) {
		os << " * * * warning, the following functions are untested * * *\n";
		os << '\t' << "ipfColor\n";//need a good test for ipfColor (maybe a small gridding of orientation space for an unusual reference direction)
		
		return testBuild      (os) &&
		       testSgSym<Real>(os) &&
		       testRel  <Real>(os) &&
		       testConv       (os) &&
		       testRot  <Real>(os);
	}
}
