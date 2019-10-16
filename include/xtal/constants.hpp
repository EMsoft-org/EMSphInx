/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                     *
 * Copyright (c) 2019, William C. Lenthe                               *
 * All rights reserved.                                                *
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
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef _XTAL_CONSTANTS_
#define _XTAL_CONSTANTS_

#include "../constants.hpp"

#include <limits>

#define XTAL_USE_H5

namespace xtal {
	const int& pijk = emsphinx::pijk;

	//@brief: helper class to hold commonly used floating point constants
	template <typename Real>
	struct Constants {
		//make sure that the Real type has the appropriate attributes
		static_assert(std::is_floating_point<Real>::value, "Real must be a floating point type");//disallow integers
		static_assert(std::numeric_limits<Real>::has_infinity, "Real must have infinity");//must have ieee infinity (if this is relaxed inf should be set to max fp)

		//simple mathematical constants
		static const Real pi   ;//pi
		static const Real pi2  ;//pi*2
		static const Real pi_2 ;//pi/2
		static const Real r2   ;//sqrt(2)
		static const Real r3   ;//sqrt(3)
		static const Real r1_2 ;//sqrt(1/2)
		static const Real r1_3 ;//sqrt(1/3)
		static const Real r3_4 ;//sqrt(3/4) == sqrt(3) / 2

		//conversion constants
		static const Real dg2rd;//degrees to radians (pi / 180)
		static const Real rd2dg;//radians to degrees (180 / pi)

		//numerical precision constants	
		static const Real eps  ;//machine epsilon
		static const Real rEps ;//sqrt(eps)
		static const Real thr  ;//threshold for 0 (machine epsilon * small factor)
		static const Real inf  ;//infinity

		//rotations constants
		static const Real hoR  ;//radius of homochoric sphere
		static const Real hoR2 ;//radius of homochoric sphere^2
		static const Real cuA  ;//side length of cubochoric cube
		static const Real cuA_2;//side length of cubochoric cube / 2
		static const Real pjik ;//pijk as fp
		static const Real tp_8 ;//tan(pi/8) == sqrt(2) - 1
		static const Real tp_12;//tan(pi/12) == 2 - sqrt(3)
	};

}

#include <limits>
#include <cmath>

namespace xtal {
	//constants are given with 72 decimal places which is enough for IEEE octuple precision (if you ever need to support 512 bit fp numbers you'll need to recompute these...)

	//simple mathematical constants
	template <typename Real> const Real Constants<Real>::pi    = Real(3.14159265358979323846264338327950288419716939937510582097494459230781641 );
	template <typename Real> const Real Constants<Real>::pi2   = Real(6.28318530717958647692528676655900576839433879875021164194988918461563281 );
	template <typename Real> const Real Constants<Real>::pi_2  = Real(1.57079632679489661923132169163975144209858469968755291048747229615390820 );
	template <typename Real> const Real Constants<Real>::r2    = Real(1.41421356237309504880168872420969807856967187537694807317667973799073248 );
	template <typename Real> const Real Constants<Real>::r3    = Real(1.73205080756887729352744634150587236694280525381038062805580697945193302 );
	template <typename Real> const Real Constants<Real>::r1_2  = Real(0.707106781186547524400844362104849039284835937688474036588339868995366239);
	template <typename Real> const Real Constants<Real>::r1_3  = Real(0.577350269189625764509148780501957455647601751270126876018602326483977672);
	template <typename Real> const Real Constants<Real>::r3_4  = Real(0.866025403784438646763723170752936183471402626905190314027903489725966508);

	//conversion constants
	template <typename Real> const Real Constants<Real>::dg2rd = Real(0.0174532925199432957692369076848861271344287188854172545609719144017100911);
	template <typename Real> const Real Constants<Real>::rd2dg = Real(57.2957795130823208767981548141051703324054724665643215491602438612028471  );

	//numerical precision constants	
	template <typename Real> const Real Constants<Real>::eps   = std::numeric_limits<Real>::epsilon();
	template <typename Real> const Real Constants<Real>::rEps  = std::sqrt(Constants<Real>::eps);
	template <typename Real> const Real Constants<Real>::thr   = Constants<Real>::eps * Real(10);
	template <typename Real> const Real Constants<Real>::inf   = std::numeric_limits<Real>::infinity();

	//rotations constants
	template <typename Real> const Real Constants<Real>::hoR   = Real(1.33067003949146879092560751376547158439759677721615056766447114777522727);//(pi * 3 / 4) ^ (1/3)
	template <typename Real> const Real Constants<Real>::hoR2  = Real(1.77068275400022711161806356533586473059468238635461597721935515693649391);//(pi * 3 / 4) ^ (2/3)
	template <typename Real> const Real Constants<Real>::cuA   = Real(2.14502939711102560007744410094123559748666736547155699029161726142592328);//pi ^ (2/3)
	template <typename Real> const Real Constants<Real>::cuA_2 = Real(1.07251469855551280003872205047061779874333368273577849514580863071296164);//pi ^ (2/3)
	template <typename Real> const Real Constants<Real>::pjik  = Real(pijk);
	template <typename Real> const Real Constants<Real>::tp_8  = Real(0.41421356237309504880168872420969807856967187537694807317667973799073248);
	template <typename Real> const Real Constants<Real>::tp_12 = Real(0.26794919243112270647255365849412763305719474618961937194419302054806698);
}

#endif// _XTAL_CONSTANTS_
