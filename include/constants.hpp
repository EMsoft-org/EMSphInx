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

#ifndef _CONSTANTS_H_
#define _CONSTANTS_H_

#ifndef EMSPHINX_GIT_BRANCH
	static_assert(false, "git branch isn't defined (please -DEMSPHINX_GIT_BRANCH=___ or redefine empshinx::Version");
#endif

#ifndef EMSPHINX_GIT_HASH
	static_assert(false, "git hash isn't defined (please -DEMSPHINX_GIT_HASH=___ or redefine empshinx::Version");
#endif

#define EM_STRINGIFY(x) #x
#define EM_2_STR(x) EM_STRINGIFY(x)

#include <string>

namespace emsphinx {
	//@remarks
	//crystallographic orientations in this code are described as passive rotations from the sample to the crystal frame
	//this is the most common choice for EBSD codes and crystallography codes in the US and consistent with EMsoft
	//if your use case requires another choice you'll need to adjust the output accordingly

	//axis angle convention switch
	//euler angles are always treated in the passive sense such that
	//eu[3] = {alpha, beta, gamma} is a rotation of the reference frame by:
	// -first alpha about the z axis
	// -next beta about the y' axis
	// -finally gamma bout the z'' axis
	//pijk (+/-1) is \hat{i} * \hat{j} * \hat{k} (the complex components of the quaternion)
	//for pijk = +1 a passive rotation of 120@[1,1,1] -> {0.5, +0.5 * \hat{i}, +0.5 * \hat{j}, +0.5 * \hat{k}}
	//for pijk = -1 a passive rotation of 120@[1,1,1] -> {0.5, -0.5 * \hat{i}, -0.5 * \hat{j}, -0.5 * \hat{k}}
	const int pijk = +1;
	static_assert(pijk == +1 || pijk == -1, "pijk must be +/-1");

	//@brief: helper class to hold commonly used floating point constants
	template <typename Real>
	struct Constants {
		static_assert(std::is_floating_point<Real>::value, "Real must be a floating point type");//disallow integers
		static const Real pi  ;//pi
		static const Real pi2 ;//pi*2
		static const Real pi_2;//pi/2
	};

	static const std::string GitBranch = EM_2_STR(EMSPHINX_GIT_BRANCH);
	static const std::string GitHash   = EM_2_STR(EMSPHINX_GIT_HASH);
	static const std::string Version   = GitBranch + ':' + GitHash;
}

namespace emsphinx {
	//constants are given with 72 decimal places which is enough for IEEE octuple precision (if you ever need to support 512 bit fp numbers you'll need to recompute these...)
	template <typename Real> const Real Constants<Real>::pi   = Real(3.14159265358979323846264338327950288419716939937510582097494459230781641);
	template <typename Real> const Real Constants<Real>::pi2  = Real(6.28318530717958647692528676655900576839433879875021164194988918461563281);
	template <typename Real> const Real Constants<Real>::pi_2 = Real(1.57079632679489661923132169163975144209858469968755291048747229615390820 );
}

#endif//_CONSTANTS_H_
