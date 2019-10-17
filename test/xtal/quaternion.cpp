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

// tests for xtal/quaternion.cpp

//@note : currently no test for quat::less

#include <iostream>

namespace xtal {
	//@brief : quaternion unit tests
	//@return: true / false if tests pass / fail
	template <typename Real> bool testQuats(std::ostream& os);
}

int main() {
	std::ostream& os = std::cout;
	return xtal::testQuats<float>(os) && xtal::testQuats<double>(os) ? EXIT_SUCCESS : EXIT_FAILURE;
}

#include <limits>
#include <cmath>
#include <algorithm>

#include "xtal/quaternion.hpp"

namespace xtal {


	//@brief : quaternion tests
	//@return: true / false if tests pass / fail
	template <typename Real> bool testQuats(std::ostream& os) {
		//make 2 random quats
		Real qu[4] = {1,  2, -3,  4};
		Real qr[4] = {5, -6,  7, -8};
		Real qv[4];
		const Real eps = std::sqrt(std::numeric_limits<Real>::epsilon());

		//check magnitude functions
		Real v;
		v = quat::norm2(qu);
		if(v != Real(30)) {
			os << "quat norm2 of (1, 2, -3, 4) != 30\n";
			return false;
		}
		v = std::fabs(quat::norm(qr) - std::sqrt(Real(174)));
		if(v > eps) {
			os << "quat norm of (5, -6, 7, -8) != sqrt(174)\n";
			return false;
		}

		//check dot product
		v = quat::dot(qu, qr);
		if(v != Real(-60)) {
			os << "quat dot of (1, 2, -3, 4), (5, -6, 7, -8) != 60\n";
			return false;
		}

		//check comparison
		if(quat::equal(qu, qr) || !quat::equal(qu, qu) || !quat::equal(qr, qr)) {
			os << "quat::equal failed\n";
			return false;
		}

		//check element wise functions
		quat::scalarAdd(qu, Real(-5), qv);
		if(Real(-4) != qv[0] || Real(-3) != qv[1] || Real(-8) != qv[2] || Real(-1) != qv[3]) {
			os << "quat::scalarAdd failed\n";
			return false;
		}
		quat::scalarSub(qv, Real(-5), qv);
		if(!quat::equal(qu, qv)) {
			os << "quat::scalarSub failed\n";
			return false;

		}
		quat::scalarDiv(qv, Real(-0.5), qv);
		if(Real(-2) != qv[0] || Real(-4) != qv[1] || Real(6) != qv[2] || Real(-8) != qv[3]) {
			os << "quat::scalarDiv failed\n";
			return false;
		}
		quat::scalarMul(qv, Real(-0.5), qv);
		if(!quat::equal(qu, qv)) {
			os << "quat::scalarMul failed\n";
			return false;
		}

		//check addition and subtraction
		quat::add(qu, qr, qv);
		if(Real(6) != qv[0] || Real(-4) != qv[1] || Real(4) != qv[2] || Real(-4) != qv[3]) {
			os << "quat::add failed\n";
			return false;
		}
		quat::sub(qv, qr, qv);
		if(!quat::equal(qu, qv)) {
			os << "quat::sub failed\n";
			return false;
		}

		//check simple urnary ops
		quat::cAbs(qu, qv);
		if(Real(1) != qv[0] || Real(2) != qv[1] || Real(3) != qv[2] || Real(4) != qv[3]) {
			os << "quat::cAbs failed\n";
			return false;
		}
		quat::conj(qu, qv);
		if(Real(1) != qv[0] || Real(-2) != qv[1] || Real(3) != qv[2] || Real(-4) != qv[3]) {
			os << "quat::conj failed\n";
			return false;
		}
		quat::neg(qv, qv);
		qv[0] = -qv[0];
		if(!quat::equal(qu, qv)) {
			os << "quat::neg failed\n";
			return false;
		}
		quat::conj(qu, qv);
		quat::neg (qv, qv);
		quat::expl(qv, qv);
		if(quat::equal(qu, qv)) {
			os << "quat::expl failed\n";
			return false;
		}

		//check multiplication
		if(1 == emsphinx::pijk) {//since pijk isn't in definition of qu / qr
			quat::conj(qu, qu);//account for original test definition w/ pijk == -1
			quat::conj(qr, qr);//account for original test definition w/ pijk == -1
			quat::mul (qu, qr, qv);
			quat::conj(qu, qu);//undo conjugation
			quat::conj(qr, qr);//undo conjugation
			quat::conj(qv, qv);//account for original test definition w/ pijk == -1
		} else {
			quat::mul(qu, qr, qv);
		}
		if(Real(70) != qv[0] || Real(8) != qv[1] || Real(0) != qv[2] || Real(16) != qv[3]) {
			os << "quat::mul failed\n";
			return false;
		}

		//check inverse
		quat::inv(qu, qv);
		quat::mul(qu, qv, qv);
		if(std::fabs(qv[0] - Real(1)) > eps) {
			os << "quat::inv failed\n";
			return false;
		}
		quat::cAbs(qv, qv);
		if(std::max(std::max(qv[1], qv[2]), qv[3]) > eps) {
			os << "quat::inv failed\n";
			return false;
		}

		//check division
		quat::mul(qu, qr, qv);//(qu * qr)
		quat::div(qv, qr, qv);//(qu * qr)
		if(!quat::equal(qu, qv)) {
			os << "quat::div failed\n";
			return false;
		}

		//check normalization
		quat::normalize(qu, qv);
		v = std::fabs(quat::norm(qv) - Real(1));
		if(v > eps) {
			os << "quat::normalize failed\n";
			return false;
		}

		//check active vector rotation
		qv[0] = Real(0.5);//120 degrees
		v = std::sqrt(Real(56) / 3);
		qv[1] = Real(-1 * pijk) / v;
		qv[2] = Real( 2 * pijk) / v;
		qv[3] = Real(-3 * pijk) / v;
		//qv is now passive rotation of 120 degrees @ {-1, 2, -3}
		//rotate {5, -6, 7}
		const Real vr[3] = {//this should be the simplified result (reguardless of pijk)
			Real( 11) / 7 - std::sqrt(Real( 6) / 7),
			Real(-36) / 7 - std::sqrt(Real(24) / 7),
			Real( 61) / 7 - std::sqrt(Real( 6) / 7),
		};
		quat::rotateVector(qv, qr, qv);
		v = std::max(std::fabs(qv[0] - vr[0]), std::max(std::fabs(qv[1] - vr[1]), std::fabs(qv[2] - vr[2])));
		if(v > eps) {
			os << "quat::rotateVector failed\n";
			return false;
		}

		//if we made it this far all tests pass
		os << "all tests passed\n";
		return true;
	}
}
