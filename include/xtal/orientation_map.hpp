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

#ifndef _OrientationMap_H_
#define _OrientationMap_H_

#include <string>
#include <vector>

#include "phase.hpp"
#include "quaternion.hpp"

#define XTAL_USE_TSL
#define XTAL_USE_HKL

#ifdef XTAL_USE_TSL
	#include "vendor/tsl.hpp"
#endif

#ifdef XTAL_USE_HKL
	#include "vendor/hkl.hpp"
#endif

namespace ebsd {
	//a subset of the complete geometry  (we're lucky if the vendors give us this much)
	template <typename Real = float>
	struct Calibration {
		//@brief: enumeration of normalized pattern center types
		enum class Vendor {
			Bruker,
			EDAX  ,
			Oxford,
		};

		//experimental conditions
		Real sTlt ;//sample tilt in degrees
		Real cTlt ;//camera tilt in degrees
		Real wd   ;//sample working distance in mm
		Real kv   ;//accelerating voltage in kV

		//pattern center
		//@note: xStar is the same for TSL, Oxford and Bruker
		//       Oxford yStar is is yStar(TSL) * detector width / detector height
		//       Bruker yStar is 1 - yStar(Oxford)
		//       zStar is the same for TSL and Oxford, Bruker zStar is zStar(TSL) * detector height / detector width
		Vendor ven;//pattern center vendor
		Real ratio;//detector aspect ratio (width / height)
		Real xStar;//normalized pattern center x
		Real yStar;//normalized pattern center y
		Real zStar;//normalized pattern center z

		//@brief: construct an uninitialized calibration
		Calibration() : sTlt(NAN), cTlt(NAN), wd(0), kv(0), xStar(NAN), yStar(NAN), zStar(NAN) {}

		//@brief  : convert from current vendor to new vendor
		//@param v: new vendor
		//@note   : correctly updates pattern center using aspect ratio
		void setVendor(const Vendor v);
	};
}

namespace xtal {
	//a 2D orientation map on a square grid
	template <typename Real>
	struct OrientationMap {
		//scan conditions
		uint32_t                   width  ;//scan width  in pixels
		uint32_t                   height ;//scan height in pixels
		Real                       xStep  ;//width  of pixel in microns
		Real                       yStep  ;//height of pixel in microns
		ebsd::Calibration<Real>    calib  ;//detector geometry

		//meta data
		std::string                owner  ;//name of person who collected the data
		std::string                name   ;//description of sample, project, etc

		//possible phases
		std::vector< Phase<Real> > phsList;//list of phases

		//pixel information
		std::vector< uint_fast8_t> phase  ;//index of phase at each point
		std::vector< Quat<Real>  > qu     ;//orientation at each point (as quaternions [w,x,y,z])
		std::vector< Real        > metric ;//scalar indexing quality measure at each point (e.g. CI for tsl)
		std::vector< Real        > imQual ;//scalar pattern quality measure at each point (e.g. IQ)

		//@brief   : create an empty orientation map
		//@param w : width of scan in pixels
		//@param h : height of scan in pixels
		//@param dx: pixel width in microns
		//@param dy: pixel height in microns
		OrientationMap(const uint32_t w, const uint32_t h, const Real dx = 1, const Real dy = 1);

	#ifdef XTAL_USE_TSL
		//@brief    : extract an orientation map from a tsl scan
		//@param tsl: tsl scan to extract from
		OrientationMap(const tsl::OrientationMap& om);

		//@brief : convert an orientation map to a tsl scan format
		//@return: tsl scan
		tsl::OrientationMap toTSL() const;
	#endif

	#ifdef XTAL_USE_HKL
		//@brief    : extract an orientation map from an hkl scan
		//@param hkl: hkl scan to extract from
		OrientationMap(const hkl::OrientationMap& om);

		//@brief : convert an orientation map to a hkl scan format
		//@return: hkl scan
		hkl::OrientationMap toHKL() const;
	#endif

		//@brief         : read an orientation map from a file
		//@param fileName: name of file to read from
		//@param aux     : auxilary information (for hdf5 datasets this is the path to the dataset)
		OrientationMap(const std::string fileName, const std::string aux = "") {read(fileName, aux);}

		//@brief         : read an orientation map from a file
		//@param fileName: name of file to read from
		//@param aux     : auxilary information (for hdf5 datasets this is the scan name)
		void read(const std::string fileName, const std::string aux = "");

		//@brief         : write an orientation map to a file
		//@param fileName: name of file to write to
		void write(const std::string fileName) const;

		//@brief    : write an orientation map to an hdf5 group
		//@param grp: h5 group to write to
		void writeH5(H5::Group grp) const;

		//@brief: build an ipf color map
		//@param rgb   : location to write red, green, blue (, alpha) values
		//@param refDir: reference direction (sample frame)
		//@param h2r: hsl2rgb like coloring function to use with h, s, and l in [0,1] and output as [0,1] rgb: void(Real const * const hsl, Real * const rgb)
		//@param alpha : true to write an RGBA image, false for RGB
		void ipfColor(uint8_t * const rgb, Real const * const refDir, std::function<void(Real const*const, Real*const)> h2r, const bool alpha = false) const;

		//@brief: build an ipf color map
		//@param refDir: reference direction (sample frame)
		//@param h2r: hsl2rgb like coloring function to use with h, s, and l in [0,1] and output as [0,1] rgb: void(Real const * const hsl, Real * const rgb)
		//@param alpha : true to write an RGBA image, false for RGB	
		//@return      : red, green, blue (, alpha) values
		std::vector<uint8_t> ipfColor(Real const * const refDir, std::function<void(Real const*const, Real*const)> h2r, const bool alpha = false) const;
	};
}

#include <stdexcept>

#include "rotations.hpp"
#include "constants.hpp"

namespace ebsd {

	//@brief  : convert from current vendor to new vendor
	//@param v: new vendor
	//@note   : correctly updates pattern center using aspect ratio
	template <typename Real>
	void Calibration<Real>::setVendor(const Calibration<Real>::Vendor v) {
		if(v == ven) return;//already correct
		//update y/z star (x stars are all the same)
		if       (Vendor::Bruker == ven && Vendor::EDAX   == v) {// Bruker --> EDAX
			yStar = (Real(1) - yStar ) / ratio;
			zStar = zStar / ratio;
		} else if(Vendor::Bruker == ven && Vendor::Oxford == v) {// Bruker --> Oxford
			yStar = Real(1) - yStar;
			zStar = zStar / ratio;
		} else if(Vendor::EDAX   == ven && Vendor::Oxford == v) {// EDAX   --> Oxford
			yStar = yStar * ratio;
			//TSL and Oxford zStars are the same
		} else if(Vendor::EDAX   == ven && Vendor::Bruker == v) {// EDAX   --> Bruker
			yStar = Real(1) - yStar * ratio;
			zStar = zStar * ratio;
		} else if(Vendor::Oxford == ven && Vendor::Bruker == v) {// Oxford --> Bruker
			yStar = Real(1) - yStar;
			zStar = zStar * ratio;
		} else if(Vendor::Oxford == ven && Vendor::EDAX   == v) {// Oxford --> EDAX
			yStar = yStar / ratio;
			//TSL and Oxford zStars are the same
		}
		ven = v;//save new vendor
	}

}

namespace xtal {
	namespace detail {
		//@brief         : check if a file name is a hdf5 type
		//@param fileName: name of file to check
		//@return        : true if extension is hdf, h5, or hdf5
		bool isH5(std::string fileName) {
			size_t pos = fileName.find_last_of(".");//find the last '.' in the name
			if(std::string::npos == pos) return false;//handle files with no extension
			std::string ext = fileName.substr(pos+1);//extract the file extension
			std::transform(ext.begin(), ext.end(), ext.begin(), [](char c){return std::tolower(c);});//convert to lowercase

			if     (0 == ext.compare("hdf" )) return true;
			else if(0 == ext.compare("hdf5")) return true;
			else if(0 == ext.compare("h5"  )) return true;
			else return false;
		}
	}

	//@brief   : create an empty orientation map
	//@param w : width of scan in pixels
	//@param h : height of scan in pixels
	//@param dx: pixel width in microns
	//@param dy: pixel height in microns
	template <typename Real>
	OrientationMap<Real>::OrientationMap(const uint32_t w, const uint32_t h, const Real dx, const Real dy)
		: width (w ), height(h ), 
		  xStep (dx), yStep (dy),
		  owner ("unknown"),
		  name  (""       ),
		  phase (w * h, 0                 ),
		  qu    (w * h, Quat<Real>::Zero()),
		  metric(w * h, Real(0)           ),
		  imQual(w * h, Real(0)           ) {}

#ifdef XTAL_USE_TSL
	//@brief    : extract an orientation map from a tsl scan
	//@param tsl: tsl scan to extract from
	template <typename Real>
	OrientationMap<Real>::OrientationMap(const tsl::OrientationMap& om) {
		//hex isn't allowed (could convert to square in the future)
		if(tsl::GridType::Square != om.gridType) throw std::runtime_error("hexagonal grids aren't supported");

		//extract header data
		width       = om.nColsOdd       ;//or om.nColsEven (should be the same for square grids)
		height      = om.nRows          ;
		xStep       = om.xStep          ;
		yStep       = om.yStep          ;
		calib.kv    = NAN               ;
		calib.sTlt  = om.sampTlt        ;
		calib.cTlt  = om.camTlt         ;
		calib.wd    = om.workingDistance;
		calib.ven   = ebsd::Calibration<Real>::Vendor::EDAX;
		calib.xStar = om.xStar          ;
		calib.yStar = om.yStar          ;
		calib.zStar = om.zStar          ;
		owner       = om.operatorName   ;
		name        = om.sampleId       ;
		if(!om.scanId.empty()) name += ": " + om.scanId;

		//extract phases
		for(const tsl::Phase& p : om.phaseList) {
			Phase<Real> xp;
			std::copy(p.lat, p.lat+6, xp.lat);
			for(size_t i = 0; i < 3; i++) xp.lat[i] /= 10;//angstroms -> nm
			xp.name = p.name;//currently we loose formula and info data
			xp.pg   = PointGroup::FromTSL(p.sym);
			//currently we lose hkl families, elastic constants, and cetegories
			phsList.push_back(xp);
		}

		//extract phases, orientations, and quality measure
		phase .resize(width * height);
		qu    .resize(width * height);
		metric.resize(width * height);
		imQual.resize(width * height);
		std::copy(om.phase.cbegin(), om.phase.cend(), phase .begin());
		std::copy(om.ci   .cbegin(), om.ci   .cend(), metric.begin());//could also use fit
		std::copy(om.iq   .cbegin(), om.iq   .cend(), imQual.begin());
		float const * pEu = &om.eu[0];
		Real eu[3];
		for(Quat<Real>& q : qu) {
			std::copy(pEu, pEu + 3, eu);//convert eu to Real

			//undo tsl convention
			// eu[0] = std::fmod(eu[0]+Constants<Real>::pi_2, Constants<Real>::pi2);

			pEu += 3;//increment euler angle for next iteration
			eu2qu(eu, q.data());//convert to quaternion
		}
	}

	//@brief : convert an orientation map to a tsl scan format
	//@return: tsl scan
	template <typename Real>
	tsl::OrientationMap OrientationMap<Real>::toTSL() const {
		//create map
		tsl::OrientationMap om;
		om.gridType = tsl::GridType::Square;

		//copy camera calibration
		ebsd::Calibration<Real> c = calib;
		c.setVendor(ebsd::Calibration<Real>::Vendor::EDAX);

		om.sampTlt         = (float) c.sTlt ;
		om.camTlt          = (float) c.cTlt ;
		om.xStar           = (float) c.xStar;
		om.yStar           = (float) c.yStar;
		om.zStar           = (float) c.zStar;
		om.workingDistance = (float) c.wd   ;

		//copy header info
		om.xStep        = (float) xStep ;
		om.yStep        = (float) yStep ;
		om.nColsOdd     =         width ;
		om.nColsEven    =         width ;
		om.nRows        =         height;
		om.operatorName =         owner ;
		om.sampleId     =         name  ;

		//copy phase data
		om.phaseList.resize(phsList.size());
		for(size_t i = 0; i < phsList.size(); i++) {
			om.phaseList[i].num = i+1;
			om.phaseList[i].name = phsList[i].name;
			om.phaseList[i].sym  = phsList[i].pg.tslNum();
			for(size_t j = 0; j < 6; j++) om.phaseList[i].lat[j] = (float)phsList[i].lat[j];
			for(size_t j = 0; j < 3; j++) om.phaseList[i].lat[j] *= 10;//nm to angstroms
		}

		//allocate scan data and copy
		om.allocate(8);//we dont need SEM or fit
		std::copy     (phase .cbegin(), phase .cend(), om.phase.begin());
		std::transform(metric.cbegin(), metric.cend(), om.ci   .begin(), [](const Real& r){return (float)r;});
		std::transform(imQual.cbegin(), imQual.cend(), om.iq   .begin(), [](const Real& r){return (float)r;});

		//convert quats to eulers
		Real eu[3];
		float * pEu = &om.eu[0];
		for(const Quat<Real>& q : qu) {
			qu2eu(q.data(), eu);//convert to euler angles
			for(size_t i = 0; i < 3; i++) pEu[i] = (float)eu[i];//convert eu to float
			pEu += 3;//increment euler angle for next iteration
		}

		//build x/y
		for(size_t i = 0; i < width; i++) om.x[i] = (float)(xStep * i);
		for(size_t j = 0; j < height; j++) {
			std::copy(om.x.begin(), om.x.begin() + width, om.x.begin() + j * width);//copy x from first row
			std::fill(om.y.begin() + j * width, om.y.begin() + (j+1) * width, (float)(yStep * j));//fill y with current value
		}
		return om;
	}
#endif//XTAL_USE_TSL

#ifdef XTAL_USE_HKL

	//@brief    : extract an orientation map from an hkl scan
	//@param hkl: hkl scan to extract from
	template <typename Real>
	OrientationMap<Real>::OrientationMap(const hkl::OrientationMap& om) {
		//hex isn't allowed (could convert to square in the future)
		if(0 != om.jobMode.compare("Grid")) throw std::runtime_error("hexagonal grids aren't supported");

		//extract header data
		width      = om.xCells ;
		height     = om.yCells ;
		xStep      = om.xStep  ;
		yStep      = om.yStep  ;
		calib.sTlt = om.angle  ;
		calib.kv   = om.kv     ;
		calib.ven  = ebsd::Calibration<Real>::Vendor::Oxford;
		owner      = om.author ;
		name       = om.project;

		//extract phases
		for(const hkl::Phase& p : om.phaseList) {
			Phase<Real> xp;
			std::copy(p.lat, p.lat+6, xp.lat);
			for(size_t i = 0; i < 3; i++) xp.lat[i] /= 10;//angstroms -> nm
			xp.name = p.name;
			xp.pg   = PointGroup(p.space);//convert from space group to point group (some loss of data)
			phsList.push_back(xp);
		}

		//extract phases, orientations, and quality measure
		phase .resize(width * height);
		qu    .resize(width * height);
		metric.resize(width * height);
		imQual.resize(width * height);
		std::copy    (om.phase.cbegin(), om.phase.cend(), phase .begin());
		std::for_each(phase   . begin(), phase   . end(), [](uint_fast8_t& i){--i;});//hkl phases are 1 indexed
		std::copy    (om.err  .cbegin(), om.err  .cend(), metric.begin());//could also use MAD
		std::copy    (om.bc   .cbegin(), om.bc   .cend(), imQual.begin());//could also use BS
		float const * pEu = &om.eu[0];
		Real eu[3];
		for(Quat<Real>& q : qu) {
			std::copy(pEu, pEu + 3, eu);//convert eu to Real
			for(size_t i = 0; i < 3; i++) eu[i] *= Constants<Real>::dg2rd;//convert from degrees to radians
			pEu += 3;//increment euler angle for next iteration
			eu2qu(eu, q.data());//convert to quaternion
		}
	}

	//@brief : convert an orientation map to a hkl scan format
	//@return: hkl scan
	template <typename Real>
	hkl::OrientationMap OrientationMap<Real>::toHKL() const {
		//create map
		hkl::OrientationMap om;
		om.jobMode = "Grid";

		//copy header info
		om.xCells  =         width     ;
		om.yCells  =         height    ;
		om.zCells  =         1         ;
		om.xStep   = (float) xStep     ;
		om.yStep   = (float) yStep     ;
		om.zStep   =         0         ;
		om.angle   = (float) calib.sTlt;
		om.kv      = (float) calib.kv  ;
		om.author  =         owner     ;
		om.project =         name      ;
		om.euStr   = "Euler angles refer to Sample Coordinate system (CS0)!";

		//copy phase data
		for(const Phase<Real>& xp : phsList) {
			hkl::Phase p;
			for(size_t i = 0; i < 6; i++) p.lat[i] = (float)xp.lat[i];
			for(size_t i = 0; i < 3; i++) p.lat[i] *= 10;//nm to angstroms
			p.name  = xp.name;
			p.laue  = xp.pg.hklNum    ();
			p.space = xp.pg.symmorphic();//get lowest symmetry space group of point group
			om.phaseList.push_back(p);
		}

		//allocate scan data and copy
		om.allocate( hkl::OrientationMap::CTF_Phase
			       | hkl::OrientationMap::CTF_Error
			       | hkl::OrientationMap::CTF_Euler
			       | hkl::OrientationMap::CTF_X
			       | hkl::OrientationMap::CTF_Y);
		std::transform(phase .cbegin(), phase .cend(), om.phase.begin(), [](const uint_fast8_t& i){return        i+1;});//0->1 indexed
		std::transform(metric.cbegin(), metric.cend(), om.err  .begin(), [](const Real        & r){return (float)r  ;});//real->float
		std::transform(imQual.cbegin(), imQual.cend(), om.bc   .begin(), [](const Real        & r){return (float)r  ;});//real->float

		//convert quats to eulers
		Real eu[3];
		float * pEu = &om.eu[0];
		for(const Quat<Real>& q : qu) {
			qu2eu(q.data(), eu);//convert to euler angles
			for(size_t i = 0; i < 3; i++) pEu[i] = (float) (eu[i] * Constants<Real>::rd2dg);//convert from radians to degrees
			pEu += 3;//increment euler angle for next iteration
		}

		//build x/y
		for(size_t i = 0; i < width; i++) om.x[i] = (float)(xStep * (i+1));
		for(size_t j = 0; j < height; j++) {
			std::copy(om.x.begin(), om.x.begin() + width, om.x.begin() + j * width);//copy x from first row
			std::fill(om.y.begin() + j * width, om.y.begin() + (j+1) * width, (float)(yStep * (j+1)));//fill y with current value
		}
		return om;
	}
#endif

	//@brief         : read an orientation map from a file
	//@param fileName: name of file to read from
	//@param aux     : auxilary information (for hdf5 datasets this is the scan name)
	template <typename Real>
	void OrientationMap<Real>::read(const std::string fileName, const std::string aux) {
		//first try to read as h5 dot product file
		if(H5::H5File::isHdf5(fileName)) {
			//open file and check if it has an "EMheader" group
			H5::H5File file(fileName.c_str(), H5F_ACC_RDONLY);
			int idx = 0;
			bool exists = false;
			file.iterateElems(".", &idx, [](int, const char* nm, void* pB)->int{
				if(0 == std::string("EMheader").compare(nm)) *(bool*)pB = true;
				return 0;
			}, &exists);

			//check for EMsoft format before vendors
			if(exists) {//this is an EMsoft file
				H5::Group grp = file.openGroup(aux + "/EBSD");//open the scan in question
				H5::Group hdr = grp.openGroup("Header");
				H5::Group clb = hdr.openGroup("Pattern Center Calibration");

				//read header
				float vf;
				hdr.openDataSet("nColumns"        ).read(&width , H5::PredType::NATIVE_UINT32 , H5::DataSpace(H5S_SCALAR));
				hdr.openDataSet("nRows"           ).read(&height, H5::PredType::NATIVE_UINT32 , H5::DataSpace(H5S_SCALAR));
				hdr.openDataSet("Step X"          ).read(&vf    , H5::PredType::NATIVE_FLOAT  , H5::DataSpace(H5S_SCALAR)); xStep       = (Real) vf;
				hdr.openDataSet("Step Y"          ).read(&vf    , H5::PredType::NATIVE_FLOAT  , H5::DataSpace(H5S_SCALAR)); yStep       = (Real) vf;
				hdr.openDataSet("Operator"        ).read( owner , H5::StrType(0, H5T_VARIABLE), H5::DataSpace(H5S_SCALAR));
				hdr.openDataSet("Scan ID"         ).read( name  , H5::StrType(0, H5T_VARIABLE), H5::DataSpace(H5S_SCALAR));
				hdr.openDataSet("Sample Tilt"     ).read(&vf    , H5::PredType::NATIVE_FLOAT  , H5::DataSpace(H5S_SCALAR)); calib.sTlt  = (Real) vf;
				hdr.openDataSet("Working Distance").read(&vf    , H5::PredType::NATIVE_FLOAT  , H5::DataSpace(H5S_SCALAR)); calib.wd    = (Real) vf;
				clb.openDataSet("x-star"          ).read(&vf    , H5::PredType::NATIVE_FLOAT  , H5::DataSpace(H5S_SCALAR)); calib.xStar = (Real) vf;
				clb.openDataSet("y-star"          ).read(&vf    , H5::PredType::NATIVE_FLOAT  , H5::DataSpace(H5S_SCALAR)); calib.yStar = (Real) vf;
				clb.openDataSet("z-star"          ).read(&vf    , H5::PredType::NATIVE_FLOAT  , H5::DataSpace(H5S_SCALAR)); calib.zStar = (Real) vf;

				//read phases
				hdr = hdr.openGroup("Phase");
				const size_t numPhase = (size_t)hdr.getNumObjs();
				phsList.resize(numPhase);
				for(size_t i = 0; i < numPhase; i++) {
					std::stringstream ss;
					ss << i + 1;
					phsList[i].readEBSD(hdr.openGroup(ss.str()));
				}

				//open data directory
				grp = grp.openGroup("Data");

				//read phase data
				phase .resize(width * height);
				{
					std::vector<uint8_t> buff(width * height);
					grp.openDataSet("Phase").read(buff.data(), H5::PredType::NATIVE_UINT8);
					std::copy(buff.cbegin(), buff.cend(), phase.begin());
				}

				//read orientations and quality
				qu    .resize(width * height);
				metric.resize(width * height);
				imQual.resize(width * height);
				{
					//read orientations
					Real eu[3];
					std::vector<float> buff(width * height * 3);
					float * const p0 = buff.data()                     ;
					float * const p1 = buff.data() + width * height    ;
					float * const p2 = buff.data() + width * height * 2;
					grp.openDataSet("Phi1").read(p0, H5::PredType::NATIVE_FLOAT);
					grp.openDataSet("Phi" ).read(p1, H5::PredType::NATIVE_FLOAT);
					grp.openDataSet("Phi2").read(p2, H5::PredType::NATIVE_FLOAT);
					for(size_t i = 0; i < phase.size(); i++) {
						eu[0] = (Real) p0[i];
						eu[1] = (Real) p1[i];
						eu[2] = (Real) p2[i];
						eu2qu(eu, qu[i].data());
					}

					//read indexing quality
					grp.openDataSet("CI").read(p0, H5::PredType::NATIVE_FLOAT);
					for(size_t i = 0; i < phase.size(); i++) metric[i] = (Real) p0[i];

					//read image quality
					grp.openDataSet("IQ").read(p0, H5::PredType::NATIVE_FLOAT);
					for(size_t i = 0; i < phase.size(); i++) imQual[i] = (Real) p0[i];
				}
				return;

			}
		}

	#ifdef XTAL_USE_TSL
		//first try tsl readers
		if(tsl::OrientationMap::CanRead(fileName)) {
			*this = OrientationMap(tsl::OrientationMap(fileName, aux));
			return;
		}
	#endif

		//next try hkl readers
	#ifdef XTAL_USE_HKL
		if(hkl::OrientationMap::CanRead(fileName)) {
			*this = OrientationMap(hkl::OrientationMap(fileName));
			return;
		}
	#endif

		//if we haven't read the file yet no reader was found
		throw std::runtime_error("couldn't find a reader for " + fileName);
	}


	//@brief         : read an orientation map from a file
	//@param fileName: name of file to read from
	//@brief         : write an orientation map to a file
	//@param fileName: name of file to write to
	template <typename Real>
	void OrientationMap<Real>::write(const std::string fileName) const {
		//first try to write as an h5 file
		if(detail::isH5(fileName)) {
			std::string folder = "Scan 1";
			H5::H5File file(fileName.c_str(), H5F_ACC_TRUNC);
			writeH5(file.createGroup(folder));
			return;
		}

	#ifdef XTAL_USE_TSL
		//first try tsl writers
		try {
			toTSL().write(fileName);
			return;
		} catch(...) {
		}
	#endif

		//next try hkl writers
	#ifdef XTAL_USE_HKL
		try {
			toHKL().write(fileName);
			return;
		} catch(...) {
		}
	#endif

		//if we haven't read the file yet no reader was found
		throw std::runtime_error("couldn't find a writer for " + fileName);
	}

	//@brief    : write an orientation map to an hdf5 group
	//@param grp: h5 group to write to
	template <typename Real>
	void OrientationMap<Real>::writeH5(H5::Group grp) const {
		//make primary folders
		H5::Group ebsd = grp .createGroup("EBSD"  );
		H5::Group hdr  = ebsd.createGroup("Header");
		H5::Group data = ebsd.createGroup("Data"  );

		//write header
		float vf;
		hdr.createDataSet("nColumns"        , H5::PredType::NATIVE_UINT32 , H5::DataSpace(H5S_SCALAR)).write(&width , H5::PredType::NATIVE_UINT32);
		hdr.createDataSet("nRows"           , H5::PredType::NATIVE_UINT32 , H5::DataSpace(H5S_SCALAR)).write(&height, H5::PredType::NATIVE_UINT32);
		vf = (float) xStep;
		hdr.createDataSet("Step X"          , H5::PredType::NATIVE_FLOAT  , H5::DataSpace(H5S_SCALAR)).write(&vf    , H5::PredType::NATIVE_FLOAT );
		vf = (float) yStep;
		hdr.createDataSet("Step Y"          , H5::PredType::NATIVE_FLOAT  , H5::DataSpace(H5S_SCALAR)).write(&vf    , H5::PredType::NATIVE_FLOAT );
		hdr.createDataSet("Operator"        , H5::StrType(0, H5T_VARIABLE), H5::DataSpace(H5S_SCALAR)).write( owner ,H5::StrType(0, H5T_VARIABLE));
		hdr.createDataSet("Scan ID"         , H5::StrType(0, H5T_VARIABLE), H5::DataSpace(H5S_SCALAR)).write( name  ,H5::StrType(0, H5T_VARIABLE));
		vf = (float) calib.sTlt;
		hdr.createDataSet("Sample Tilt"     , H5::PredType::NATIVE_FLOAT  , H5::DataSpace(H5S_SCALAR)).write(&vf    , H5::PredType::NATIVE_FLOAT );
		vf = (float) calib.wd;
		hdr.createDataSet("Working Distance", H5::PredType::NATIVE_FLOAT  , H5::DataSpace(H5S_SCALAR)).write(&vf    , H5::PredType::NATIVE_FLOAT );

		H5::Group clb = hdr.createGroup("Pattern Center Calibration");
		vf = (float) calib.xStar;
		clb.createDataSet("x-star"          , H5::PredType::NATIVE_FLOAT  , H5::DataSpace(H5S_SCALAR)).write(&vf    , H5::PredType::NATIVE_FLOAT );
		vf = (float) calib.yStar;
		clb.createDataSet("y-star"          , H5::PredType::NATIVE_FLOAT  , H5::DataSpace(H5S_SCALAR)).write(&vf    , H5::PredType::NATIVE_FLOAT );
		vf = (float) calib.zStar;
		clb.createDataSet("z-star"          , H5::PredType::NATIVE_FLOAT  , H5::DataSpace(H5S_SCALAR)).write(&vf    , H5::PredType::NATIVE_FLOAT );

		std::string grid("SqrGrid");
		hdr.createDataSet("Grid Type"       , H5::StrType(0, H5T_VARIABLE), H5::DataSpace(H5S_SCALAR)).write( grid  ,H5::StrType(0, H5T_VARIABLE));
		hdr.createDataSet("Notes"           , H5::StrType(0, H5T_VARIABLE), H5::DataSpace(H5S_SCALAR));
		hdr.createDataSet("Sample ID"       , H5::StrType(0, H5T_VARIABLE), H5::DataSpace(H5S_SCALAR));

		//write phases to header
		hdr = hdr.createGroup("Phase");
		// phsList.resize(numPhase);
		for(size_t i = 0; i < phsList.size(); i++) {
			std::stringstream ss;
			ss << i + 1;
			phsList[i].writeEBSD(hdr.createGroup(ss.str()));
		}

		//write data

		const hsize_t num = width * height;
		H5::DataSpace space(1, &num);
		if(std::is_same<std::uint_fast8_t,uint8_t>::value) {
			data.createDataSet("Phase", H5::PredType::NATIVE_UINT8, space).write(phase.data(), H5::PredType::NATIVE_UINT8);
		} else {
			std::vector<uint8_t> buff(phase.cbegin(), phase.cend());
			data.createDataSet("Phase", H5::PredType::NATIVE_UINT8, space).write(buff .data(), H5::PredType::NATIVE_UINT8);
		}

		//convert orientations to planar euler angles
		std::vector<float> buff(width * height * 3);
		float * const p0 = buff.data()                     ;
		float * const p1 = buff.data() + width * height    ;
		float * const p2 = buff.data() + width * height * 2;
		Real eu[3];
		for(size_t i = 0; i < phase.size(); i++) {
			qu2eu(qu[i].data(), eu);
			p0[i] = (float) eu[0];
			p1[i] = (float) eu[1];
			p2[i] = (float) eu[2];
		}

		//write orientations
		data.createDataSet("Phi1", H5::PredType::NATIVE_FLOAT, space).write(p0, H5::PredType::NATIVE_FLOAT);
		data.createDataSet("Phi" , H5::PredType::NATIVE_FLOAT, space).write(p1, H5::PredType::NATIVE_FLOAT);
		data.createDataSet("Phi2", H5::PredType::NATIVE_FLOAT, space).write(p2, H5::PredType::NATIVE_FLOAT);

		//write quality
		if(std::is_same<Real, float>::value) {
			data.createDataSet("Metric", H5::PredType::NATIVE_FLOAT, space).write(metric.data(), H5::PredType::NATIVE_FLOAT);
			data.createDataSet("IQ"    , H5::PredType::NATIVE_FLOAT, space).write(imQual.data(), H5::PredType::NATIVE_FLOAT);
		} else {
			std::transform(metric.cbegin(), metric.cend(), buff.begin(), [](const Real& v){return (float)v;});
			data.createDataSet("Metric", H5::PredType::NATIVE_FLOAT, space).write(buff  .data(), H5::PredType::NATIVE_FLOAT);
			std::transform(imQual.cbegin(), imQual.cend(), buff.begin(), [](const Real& v){return (float)v;});
			data.createDataSet("IQ"    , H5::PredType::NATIVE_FLOAT, space).write(buff  .data(), H5::PredType::NATIVE_FLOAT);
		}
	}

	//@brief: build an ipf color map
	//@param rgb   : location to write red, green, blue (, alpha) values
	//@param refDir: reference direction (sample frame)
	//@param h2r: hsl2rgb like coloring function to use with h, s, and l in [0,1] and output as [0,1] rgb: void(Real const * const hsl, Real * const rgb)
	//@param alpha : true to write an RGBA image, false for RGB
	template <typename Real>
	void OrientationMap<Real>::ipfColor(uint8_t * const rgb, Real const * const refDir, std::function<void(Real const*const, Real*const)> h2r, const bool alpha) const {
		//normalize reference direction
		const Real mag = std::sqrt(std::inner_product(refDir, refDir+3, refDir, Real(0)));
		const Real n[3] = {refDir[0] / mag, refDir[1] / mag, refDir[2] / mag};

		//build table of point groups (could just use phase[i].pg but this is more flexible for potential future optimizations)
		std::vector<PointGroup> pgs;
		for(const Phase<Real>& p : phsList) pgs.push_back(p.pg);

		//loop over pixels
		Real nx[3], color[4];
		color[3] = 255.0;
		const size_t NUM = alpha ? 4 : 3;
		for(size_t i = 0; i < qu.size(); i++) {
			if(phase[i] < phsList.size()) {
				qu[i].rotateVector(n, nx);//get direction in crystal frame
				pgs[phase[i]].ipfColor(nx, color, h2r);//get ipf color for crystal direction
				for(size_t j = 0; j < 3; j++) rgb[i * 3 + j] = (uint8_t) std::round(color[j] * 255);//[0,1] -> [0,255]
			}
		}
	}

	//@brief: build an ipf color map
	//@param refDir: reference direction (sample frame)
	//@param h2r: hsl2rgb like coloring function to use with h, s, and l in [0,1] and output as [0,1] rgb: void(Real const * const hsl, Real * const rgb)
	//@param alpha : true to write an RGBA image, false for RGB
	//@return      : red, green, blue (, alpha) values
	template <typename Real>
	std::vector<uint8_t> OrientationMap<Real>::ipfColor(Real const * const refDir, std::function<void(Real const*const, Real*const)> h2r, const bool alpha) const {
		const std::vector<uint8_t> rgb(qu.size() * (alpha ? 4 : 3));
		ipfColor((uint8_t * const)rgb.data(), refDir, h2r, alpha);
		return rgb;
	}
}

#endif//_OrientationMap_H_
