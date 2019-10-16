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

#ifndef _phase_h_
#define _phase_h_

#include "constants.hpp"

#ifdef XTAL_USE_H5
	#include "H5Cpp.h"
#endif

#include "symmetry.hpp"

namespace xtal {
	template <typename Real>
	struct Phase {
		static_assert(std::is_floating_point<Real>::value, "Phase must be templated on floating point type");

		Real        lat[6];//lattice constants: a,b,c,alpha,beta,gamma with a,b,c in nm and alpha,beta,gamma in degrees
		std::string name  ;//phase name
		PointGroup  pg    ;//point group

		//@brief: default constructor makes an empty phase
		Phase() : name("unknown"), pg(1) {std::fill(lat, lat + 6, Real(0));}

	#ifdef XTAL_USE_H5
		//@brief    : read phase data from an EMsoft master pattern HDF file
		//@param grp: folder in h5 file to read from (e.g. "/CrystalData")
		void readMaster(H5::Group grp);

		//@brief    : read phase data from an H5 ebsd file
		//@param grp: folder in h5 file to read from (e.g. "/Scan Name/EBSD/Header/Phase/3")
		void readEBSD(H5::Group grp);

		//@brief    : write phase data to an H5 ebsd file
		//@param grp: folder in h5 file to write to (e.g. "/Scan Name/EBSD/Header/Phase/3")
		void writeEBSD(H5::Group grp) const;
	#endif//XTAL_USE_H5
	};
}

////////////////////////////////////////////////////////////////////////
//                          Implementations                           //
////////////////////////////////////////////////////////////////////////

namespace xtal {

#ifdef XTAL_USE_H5
	//@brief    : read phase data from an EMsoft master pattern HDF file
	//@param grp: folder in h5 file to read from
	template <typename Real>
	void Phase<Real>::readMaster(H5::Group grp) {
		float   latParam[6];
		int32_t spaceGroup ;
		try {
			grp.openDataSet("LatticeParameters").read( latParam  , H5::PredType::NATIVE_FLOAT);
		} catch (H5::Exception&) {
			//the lattice parameters dataset may not exists for merged master patterns
			std::fill(latParam, latParam + 6, 0.0f);//
		}
		grp.openDataSet("SpaceGroupNumber" ).read(&spaceGroup, H5::PredType::NATIVE_INT32);
		std::copy(latParam, latParam + 6, lat);
		pg = PointGroup(spaceGroup);
	}

	//@brief    : read phase data from an H5 ebsd file
	//@param grp: folder in h5 file to read from (e.g. "/Scan Name/EBSD/Header/Phase/3")
	template <typename Real>
	void Phase<Real>::readEBSD(H5::Group grp) {
		float   latParam[6];
		int32_t symNum     ;
		grp.openDataSet("Lattice Constant a"    ).read(&latParam[0], H5::PredType::NATIVE_FLOAT);
		grp.openDataSet("Lattice Constant b"    ).read(&latParam[1], H5::PredType::NATIVE_FLOAT);
		grp.openDataSet("Lattice Constant c"    ).read(&latParam[2], H5::PredType::NATIVE_FLOAT);
		grp.openDataSet("Lattice Constant alpha").read(&latParam[3], H5::PredType::NATIVE_FLOAT);
		grp.openDataSet("Lattice Constant beta" ).read(&latParam[4], H5::PredType::NATIVE_FLOAT);
		grp.openDataSet("Lattice Constant gamma").read(&latParam[5], H5::PredType::NATIVE_FLOAT);
		grp.openDataSet("Symmetry"              ).read(&symNum     , H5::PredType::NATIVE_INT32);
		grp.openDataSet("MaterialName"          ).read( name       , H5::StrType(0, H5T_VARIABLE), H5::DataSpace(H5S_SCALAR));
		for(size_t i = 0; i < 3; i++) latParam[i] /= 10;//angstrom -> nm
		std::copy(latParam, latParam + 6, lat);
		pg = PointGroup::FromTSL(symNum);
	}

	//@brief    : write phase data to an H5 ebsd file
	//@param grp: folder in h5 file to write to (e.g. "/Scan Name/EBSD/Header/Phase/3")
	template <typename Real>
	void Phase<Real>::writeEBSD(H5::Group grp) const {
		float latParam[6];
		for (size_t i = 0; i < 3; i++) latParam[i] = (float)lat[i] * 10;//nm -> angstrom
		for (size_t i = 3; i < 6; i++) latParam[i] = (float)lat[i];
		int32_t symNum = pg.tslNum();
		hsize_t dim[1] = {1};
	    H5::DataSpace dSpace(1, dim);
		grp.createDataSet("Lattice Constant a"    , H5::PredType::NATIVE_FLOAT, dSpace).write(&latParam[0], H5::PredType::NATIVE_FLOAT);
		grp.createDataSet("Lattice Constant b"    , H5::PredType::NATIVE_FLOAT, dSpace).write(&latParam[1], H5::PredType::NATIVE_FLOAT);
		grp.createDataSet("Lattice Constant c"    , H5::PredType::NATIVE_FLOAT, dSpace).write(&latParam[2], H5::PredType::NATIVE_FLOAT);
		grp.createDataSet("Lattice Constant alpha", H5::PredType::NATIVE_FLOAT, dSpace).write(&latParam[3], H5::PredType::NATIVE_FLOAT);
		grp.createDataSet("Lattice Constant beta" , H5::PredType::NATIVE_FLOAT, dSpace).write(&latParam[4], H5::PredType::NATIVE_FLOAT);
		grp.createDataSet("Lattice Constant gamma", H5::PredType::NATIVE_FLOAT, dSpace).write(&latParam[5], H5::PredType::NATIVE_FLOAT);
		grp.createDataSet("Symmetry"              , H5::PredType::NATIVE_INT32, dSpace).write(&symNum     , H5::PredType::NATIVE_INT32);
		grp.createDataSet("MaterialName"          , H5::StrType(0, H5T_VARIABLE), H5::DataSpace(H5S_SCALAR)).write(name, H5::StrType(0, H5T_VARIABLE), H5::DataSpace(H5S_SCALAR));
	}
#endif//XTAL_USE_H5
}

#endif//_phase_h_
