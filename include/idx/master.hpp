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

#ifndef _MASTER_H_
#define _MASTER_H_

#include "xtal/phase.hpp"
#include "xtal/quaternion.hpp"
#include "sht/square_sht.hpp"
#include "sht_file.hpp"

namespace emsphinx {

	//@brief: base class for master pattern and master spectra
	template <typename Real>
	class MasterData {
		protected:
			std::vector< xtal::Quat<Real>   > pSm;//list of pseudo-symmetric misorientations
			xtal::Phase<Real>                 phs;//crystal phase
			Real                              sig;//sample tilt
			Real                              kv ;//accelerating voltage

		public:
			//@brief : get a copy of the phase
			//@return: phase
			xtal::Phase<Real> phase() const {return phs;}

			//@brief : get a copy of the phase's point group
			//@return: symmetry group
			xtal::PointGroup pointGroup() const {return phs.pg;}

			//@brief: get pseudo-symmetric misorientations
			const std::vector< xtal::Quat<Real> >& pseudoSym() const {return pSm;}

			//@brief   : add a pseudo-symmetric misorientation
			//@param qp: pseudo-symmetric misorientation to add as quaternion
			void addPseudoSym(const xtal::Quat<Real>& qp) {pSm.push_back(qp);}

			//@brief   : add a list of pseudo-symmetric misorientation from a file
			//@param fn: name of angle file to add misorientations from
			void addPseudoSym(std::string fn);

			//@brief : get tilt
			//@return: sample tilt in degrees
			Real getSig() const {return sig;}

			//@brief : get accelerating voltage
			//@return: accelerating voltage in kV
			Real getKv() const {return kv;}
	};

	//helper class to hold master pattern files
	template <typename Real>
	struct MasterPattern : public MasterData<Real> {
		std::vector<Real>  nh ;//north hemisphere of master pattern as square lambert projection
		std::vector<Real>  sh ;//north hemisphere of master pattern as square lambert projection
		uint16_t           dim;//side length of master pattern
		square::Layout     lyt;//grid type
		
		//@brief: construct an empty master pattern
		MasterPattern(const size_t d = 0) : dim((uint16_t)d), lyt(square::Layout::Lambert), nh(d*d), sh(d*d) {}

		//@brief         : construct a master pattern from an EMsoft output file
		//@param fileName: name of EMsoft hdf5 file to read from
		MasterPattern(std::string fileName) {read(fileName);}

		//@brief : check if a master pattern can be rescaled
		//@return: true/false if the layout of the pattern can/can't be rescaled
		bool canRescale() const {return square::Layout::Lambert == lyt;}

		//@brief    : rescale a master pattern
		//@param nDm: new side length of master pattern
		//@return   : this
		//@note     : throws an exception if !canRescale()
		MasterPattern& resize(const size_t nDm);

		//@brief    : convert from a square lambert to a square legendre grid
		//@param nDm: new side length of master pattern
		//@return   : this
		//@note     : just bilinearly interpolates between grids
		MasterPattern& toLegendre(const size_t nDm);
		MasterPattern& toLegendre() {return toLegendre(dim);}

		//@brief    : convert from a square legendre to a square lambert grid
		//@param nDm: new side length of master pattern
		//@return   : this
		//@note     : just bilinearly interpolates between grids
		MasterPattern& toLambert(const size_t nDm);
		MasterPattern& toLambert() {return toLambert(dim);}

		//@brief         : read master patterns from an EMsoft output file
		//@param fileName: name of EMsoft hdf5 file to read from
		void read(std::string fileName);

		//@brief  : give a master pattern n fold symmetry about the z axis
		//@param n: rotational order
		//@param m: 0 for no mirror
		//          1 for mirror at phi = i * 180 / n (e.g. 3m1 or -6m2)
		//          2 for mirror at phi = i * 180 / n + 90 / n (i.e. -42m, 31m, -62m, -43m)
		//@note   : uses first 1/n of ring (axact for n == 2 or 4, approximate otherwise)
		//@note   : doesn't update crystal structure
		void makeNFold(const size_t n, const int m = 0);

		//@brief  : give a master pattern a mirror plane at the quator
		//@note   : doesn't update crystal structure
		void makeZMir() {sh = nh;}

		//@brief  : give a master pattern inversion symmetry
		//@note   : doesn't update crystal structure
		void makeInvSym();

		//@brief  : copy the equator of a master pattern from the north to south hemisphere
		void matchEquator();
	};

	//helper struct to hold the spherical harmonic transform of a master pattern
	template <typename Real>
	class MasterSpectra : public MasterData<Real> {
		uint16_t                          mBw;//maximum bandwidth of harmonic transform
		std::vector< std::complex<Real> > alm;//harmonic transform coefficients with a^l_m at alm[bw * m + l]

		public:
			//@brief: empty spectra
			MasterSpectra() {}

			//@brief    : construct a spectra from a master pattern
			//@param mp : master pattern
			//@param bw : desired bandwidth
			//@param nrm: should the pattern be normalized (mean 0, stdev 1) before computing the SHT
			//@return   : spectra of square legendre pattern
			//@note     : dim must be odd, bandwidth is dim - 2
			MasterSpectra(MasterPattern<Real> mp, const uint16_t bw, const bool nrm = true);

			//@brief         : construct a master pattern from an EMsoft output file
			//@param fileName: name of EMsoft hdf5 file to read from
			MasterSpectra(std::string fileName) {read(fileName);}

			//@brief : get a pointer to the spectra with (m,l) stored: (0,0), (0,1), (0,2), ..., (0,bw-1), (1,0), (1,1), (1,2), ..., (bw-1,0), (bw-1,1), (bw-1,bw-1)
			//@return: pointer with a^l_m at ptr[bw * m + l]
			std::complex<Real>       * data()       {return alm.data();}
			std::complex<Real> const * data() const {return alm.data();}

			//@brief   : repack the spectra into a new maximum bandwidth (0 padded for increasing size)
			//@param bw: bandwidth of new spectra
			//@return  : this spectra (repacked)
			MasterSpectra& resize(const uint16_t bw);

			//@brief : get rotational symmetry of the master pattern about the z axis
			//@return: rotational symmetry
			size_t nFold() const {return MasterData<Real>::pointGroup().zRot();}

			//@brief : check if there is a mirror plane at the equator of the master pattern
			//@return: true/false if there is/isn't a mirror plane
			bool mirror() const {return MasterData<Real>::pointGroup().zMirror();}

			//@brief : check if there is a inversion symmetry
			//@return: true/false if there is/isn't inversion symmetry
			bool invSym() const {return MasterData<Real>::pointGroup().inversion();}

			//@brief : get the maximum bandwidth
			//@return: max bandwidth
			size_t getBw() const {return mBw;}

			//@brief: remove DC value from spectra (make average value 0)
			void removeDC() {alm[0] = std::complex<Real>(0);}

			//@brief         : read master patterns from an EMSphInx harmonics file
			//@param fileName: name of EMSphInx spx file to read from
			void read(std::string fileName);

	};

}//emsphinx

////////////////////////////////////////////////////////////////////////
//                          Implementations                           //
////////////////////////////////////////////////////////////////////////

#include "H5Cpp.h"
#include "util/image.hpp"
#include "xtal/vendor/emsoft.hpp"

namespace emsphinx {

	//@brief   : add a list of pseudo-symmetric misorientation from a file
	//@param fn: name of angle file to add misorientations from
	template <typename Real>
	void MasterData<Real>::addPseudoSym(std::string fn) {
		emsoft::AngleFile<double> af(fn);
		if(xtal::Rotation::Quaternion != af.getType()) throw std::runtime_error("only quaternion angle files are supported for psuedo symmetry");
		xtal::Quat<double> const* qu = (xtal::Quat<double> const*)af.getAngles();
		const xtal::Quat<double> ident(1, 0, 0, 0);
		for(size_t i = 0; i < af.getNum(); i++) {
			if(!(ident == qu[i])) {
				addPseudoSym(qu[i]);
			}
		}
	}

	////////////////////////////////////////////////////////////////////////
	//                           MasterPattern                            //
	////////////////////////////////////////////////////////////////////////

	//@brief         : read master patterns from an EMsoft output file
	//@param fileName: name of EMsoft hdf5 file to read from
	template <typename Real>
	void MasterPattern<Real>::read(std::string fileName) {
		//first open the h5 file
		H5::H5File file = H5::H5File(fileName.c_str(), H5F_ACC_RDONLY);//read only access
		std::vector<hsize_t> dims, dims2;//vector to hold dimensions of datasets

		//read accelerating voltage and sample tilt
		double dVal;
		file.openDataSet("/NMLparameters/MCCLNameList/sig" ).read(&dVal, H5::PredType::NATIVE_DOUBLE, H5::DataSpace(H5S_SCALAR)); MasterData<Real>::sig = (Real) dVal;
		file.openDataSet("/NMLparameters/MCCLNameList/EkeV").read(&dVal, H5::PredType::NATIVE_DOUBLE, H5::DataSpace(H5S_SCALAR)); MasterData<Real>::kv  = (Real) dVal;

	/*
		//sanity check atom data
		dims.resize(2, 0);
		H5::DataSet xtal = file.openDataSet("CrystalData/AtomData");
		dims.resize(xtal.getSpace().getSimpleExtentNdims());//allocate space for AtomData (should be 2)
		xtal.getSpace().getSimpleExtentDims(dims.data());//read extent in each dimension
		if(2 != dims.size() || 5 != dims[0]) throw std::runtime_error("unexpected AtomData layout");
	*/
		//construct phase from lattice parameters + space group
		MasterData<Real>::phs.readMaster(file.openGroup("/CrystalData"));
		// const size_t numAtoms = dims[1];

		//read monte carlo simulation results
		H5::DataSet energy = file.openDataSet("EMData/MCOpenCL/accum_e");
		dims.resize(energy.getSpace().getSimpleExtentNdims());//allocate space for accum_e (should be 3)
		energy.getSpace().getSimpleExtentDims(dims.data());//read extent in each dimension
		if(3 != dims.size()) throw std::runtime_error("unexpected accum_e layout");
		const size_t slicePts = (size_t)(dims[0] * dims[1]);
		std::vector<uint32_t> accumE(slicePts * (size_t)dims[2]);//allocate space
		energy.read(accumE.data(), H5::PredType::NATIVE_UINT32);//read accum_e

		//use results to determine energy averaging
		std::vector<uint32_t> eCounts((size_t)dims[2], 0);
		for(size_t i = 0; i < slicePts; i++) std::transform(eCounts.cbegin(), eCounts.cend(), accumE.cbegin() + i * eCounts.size(), eCounts.begin(), std::plus<uint32_t>());
		std::vector<Real> weights(eCounts.cbegin(), eCounts.cend());
		const Real sum = std::accumulate(weights.cbegin(), weights.cend(), Real(0));
		std::for_each(weights.begin(), weights.end(), [sum](Real& i){i /= sum;});

		//check if we're in an EBSD or ECP file
		int type = 0, idx = 0;
		file.openGroup("EMData").iterateElems(".", &idx, [](int, const char* nm, void* pT)->int{
			if(0 == std::string("EBSDmaster").compare(nm)) *(int*)pT = 1;
			if(0 == std::string("ECPmaster" ).compare(nm)) *(int*)pT = 2;
			return 0;
		}, &type);

		std::string mode;
		if(1 == type) {//EBSD
			mode = "EBSD";
		} else if(2 == type) {//ECP
			mode = "ECP";
		} else {
			throw std::runtime_error("couldn't determine if master pattern was EBSD or ECP");
		}

		//open master pattern arrays and sanity check dimensions
		const std::string prefex = "EMData/" + mode + "master/";
		H5::DataSet nhData = file.openDataSet(prefex + "mLPNH");//northern hemisphere as square lambert
		H5::DataSet shData = file.openDataSet(prefex + "mLPSH");//southern hemisphere as square lambert
		dims .resize(nhData.getSpace().getSimpleExtentNdims());//allocate space for master patterns (should be 4)
		dims2.resize(shData.getSpace().getSimpleExtentNdims());//allocate space for master patterns (should be 4)
		nhData.getSpace().getSimpleExtentDims(dims .data());//read extent in each dimension
		shData.getSpace().getSimpleExtentDims(dims2.data());//read extent in each dimension
		if(dims != dims2) throw std::runtime_error("master pattern hemispehres are different shapes");
		if(1 == type) {//EBSD
			//{atom, energy, x, y}
			if(4 != dims.size()) throw std::runtime_error("unexpected master pattern layout");
		} else if(2 == type) {//ECP
			//{atom, x, y}
			if(3 != dims.size()) throw std::runtime_error("unexpected master pattern layout");
			dims.insert(dims.begin() + 1, 1);//{atom, 1, x, y}
			weights = std::vector<Real>(1,1);//energy weights
		} else {
			throw std::runtime_error("couldn't determine if master pattern was EBSD or ECP");
		}
		const size_t numAtoms = (size_t)dims[0];
		// if(dims[0] != numAtoms) throw std::runtime_error("unexpected master pattern layout");//make sure number of patterns matches atom count
		if((size_t)dims[1] != weights.size()) throw std::runtime_error("unexpected master pattern layout");//make sure energy bins matches MC output

		//save size / allocate data
		dim = (uint16_t)dims[2];//extract master pattern side length
		nh  = std::vector<Real>(size_t(dim) * size_t(dim), 0);
		sh  = std::vector<Real>(size_t(dim) * size_t(dim), 0);
		lyt = square::Layout::Lambert;

		//now read entire master pattern
		const size_t hemPts = (size_t)(dims[2] * dims[3]);//points for one hemisphere
		const size_t atmPts = (size_t)(dims[1] * hemPts );//points for all energy bins of a single atom
		std::vector<float> nhPat((size_t)dims[0] * atmPts), shPat((size_t)dims[0] * atmPts);//allocate space
		nhData.read(nhPat.data(), H5::PredType::NATIVE_FLOAT);//read north hemisphere
		shData.read(shPat.data(), H5::PredType::NATIVE_FLOAT);//read south hemisphere

		//add together all atoms for each energy
		for(size_t i = 1; i < numAtoms; i++) {//loop over each atom array accumulating
			std::transform(nhPat.cbegin(), nhPat.cbegin() + atmPts, nhPat.cbegin() + i * atmPts, nhPat.begin(), std::plus<float>());
			std::transform(shPat.cbegin(), shPat.cbegin() + atmPts, shPat.cbegin() + i * atmPts, shPat.begin(), std::plus<float>());
		}

		//finaly compute energy weighted average
		for(size_t i = 0; i < weights.size(); i++) {
			const Real w = weights[i];
			auto func = [w](const Real& a, const float& b){return a + w * b;};
			std::transform(nh.cbegin(), nh.cend(), nhPat.cbegin() + hemPts * i, nh.begin(), func);
			std::transform(sh.cbegin(), sh.cend(), shPat.cbegin() + hemPts * i, sh.begin(), func);
		}
	}

	//@brief    : rescale a master pattern
	//@param nDm: new side length of master pattern
	//@return   : this
	//@note     : throws an exception if !canRescale()
	template <typename Real>
	MasterPattern<Real>& MasterPattern<Real>::resize(const size_t nDm) {
		if(!canRescale()) throw std::runtime_error("can't rescale non-lambert master patterns");
		if(nDm == dim) return *this;

		//rescale hemispheres into new storage
		std::vector<Real> nhScaled(nDm*nDm), shScaled(nDm*nDm);//allocate space for scaled patterns
		image::Rescaler<Real> sclr(dim, dim, Real(nDm) / dim, fft::flag::Plan::Estimate);//build rescaler, we're only rescaling as a prestep so don't waste time planning
		sclr.scale(nh.data(), nhScaled.data(), false);//rescale north hemisphere
		sclr.scale(sh.data(), shScaled.data(), false);//rescale south hemisphere

		//correct fftw rescaling
		const Real factor = Real(0.5) / nhScaled.size();
		std::for_each(nhScaled.begin(), nhScaled.end(), [factor](Real& v){v *= factor;});
		std::for_each(shScaled.begin(), shScaled.end(), [factor](Real& v){v *= factor;});

		//update members
		nh.swap(nhScaled);
		sh.swap(shScaled);
		dim = (uint16_t) nDm;
		return *this;
	}

	//@brief    : convert from a square lambert to a square legendre grid
	//@param nDm: new side length of master pattern
	//@note     : just bilinearly interpolates between grids
	//@return   : this
	template <typename Real>
	MasterPattern<Real>& MasterPattern<Real>::toLegendre(const size_t nDm) {
		//sanity check grid type
		if(square::Layout::Legendre == lyt) return *this;//already legendre
		else if(square::Layout::Lambert != lyt) throw std::logic_error("unknown square grid type");

		//rescale larger for better interpolation
		const size_t dimScaled = (size_t)std::round(std::sqrt(Real(2)) * nDm);
		resize(dimScaled);//rescale pattern if needed

		//compute normals of square legendre grid
		const size_t nPts = nDm * nDm;
		std::vector<Real> xyz(nPts * 3);//allocate space
		square::legendre::normals(nDm, xyz.data());//compute normals for northern hemisphere

		//interpolate master pattern on legendre grid points
		std::vector<Real> lgNh(nPts), lgSh(nPts);//allocate space for north/south hemispheres of legendre grids
		for(size_t i = 0; i < nPts; i++) {
			//square lambert project to [0, dim - 1]
			Real XY[2];
			Real const * const pt = xyz.data() + 3 * i;
			square::lambert::sphereToSquare(pt[0], pt[1], pt[2], XY[0], XY[1]);

			//bilinear interpolate from master pattern
			image::BiPix<Real> p;
			p.bilinearCoeff(XY[0], XY[1], dimScaled, dimScaled);
			lgNh[i] = p.interpolate(nh.data());
			lgSh[i] = p.interpolate(sh.data());
		}

		//update members
		nh.swap(lgNh);
		sh.swap(lgSh);
		lyt = square::Layout::Legendre;
		dim = (uint16_t)nDm;
		return *this;
	}

	//@brief    : convert from a square legendre to a square lambert grid
	//@param nDm: new side length of master pattern
	//@note     : just bilinearly interpolates between grids
	//@return   : this
	template <typename Real>
	MasterPattern<Real>& MasterPattern<Real>::toLambert(const size_t nDm) {
		//sanity check grid type
		if(square::Layout::Lambert == lyt) return *this;//already lambert
		else if(square::Layout::Legendre != lyt) throw std::logic_error("unknown square grid type");
		const bool even = 0 == (dim % 2);

		//compute cosines of legendre grid once
		std::vector<Real> cLat = square::cosLats<Real>(dim, square::Layout::Legendre);
		std::vector<Real> xyz = square::normals<Real>(dim, square::Layout::Legendre);

		//allocate space for north/south hemispheres of legendre grids
		const size_t nPts = nDm * nDm;
		std::vector<Real> lmNh(nPts), lmSh(nPts);

		//loop over new grid points bilinearly interpolating
		Real XY[2];
		for(size_t j = 0; j < nDm; j++) {//loop over rows of square lambert grid
			XY[1] = Real(j) / (nDm - 1);//fractional y position
			const Real aY = std::fabs(XY[1] - Real(0.5));//convert from [0,1] to | [-0.5,0.5] |
			for(size_t i = 0; i < nDm; i++) {//loop over
				XY[0] = Real(i) / (nDm - 1);//fractional x position

				//now determing indices of bounding points on legendre grid
				Real n[3];
				size_t inds[4];
				square::lambert::squareToSphere(XY[0], XY[1], n[0], n[1], n[2]);
				square::legendre::boundingInds(dim, cLat.data(), n, inds);

				//next determine the 3 nearest points
				const Real dots[4] = {
					std::inner_product(n, n+3, xyz.data() + 3 * inds[0], Real(0)),
					std::inner_product(n, n+3, xyz.data() + 3 * inds[1], Real(0)),
					std::inner_product(n, n+3, xyz.data() + 3 * inds[2], Real(0)),
					std::inner_product(n, n+3, xyz.data() + 3 * inds[3], Real(0)),
				};

				//for now just nearest neighbor interpolate
				//in the future this should probably switch to barycentric between 3 closests points
				const size_t idx = std::distance(dots, std::max_element(dots, dots+4));
				lmNh[j*nDm+i] = nh[inds[idx]];
				lmSh[j*nDm+i] = sh[inds[idx]];
			}
		}

		//update members
		nh.swap(lmNh);
		sh.swap(lmSh);
		lyt = square::Layout::Lambert;
		dim = (uint16_t)nDm;
		return *this;
	}

	//@brief  : give a master pattern n fold symmetry about the z axis
	//@param n: rotational order
	//@param m: 0 for no mirror
	//          1 for mirror at phi = i * 180 / n (e.g. 3m1 or -6m2)
	//          2 for mirror at phi = i * 180 / n + 90 / n (i.e. -42m, 31m, -62m, -43m)
	//@note   : uses first 1/n of ring (axact for n == 2 or 4, approximate otherwise)
	//@note   : doesn't update crystal structure
	template <typename Real>
	void MasterPattern<Real>::makeNFold(const size_t n, const int m) {
		if(0 == dim % 2) throw std::runtime_error("must be odd side length");
		std::function<void(Real * const)> doHemi = [this,n,m](Real*const hemi) {
			const size_t rings = (this->dim + 1) / 2;
			for(size_t i = 1; i < rings; i++) {
				std::vector<Real> work(4 * (this->dim - 1));//large enough to hold equator
				const size_t count = square::readRing(this->dim, i, hemi, work.data());//extract nth ring
				const double repeatNum = double(count) / n;
				const size_t iRepeatNum = (size_t)repeatNum;
				const size_t iRepeatNum14 = (size_t)(repeatNum * 0.25);
				const size_t iRepeatNum24 = (size_t)(repeatNum * 0.50);
				const size_t iRepeatNum34 = (size_t)(repeatNum * 0.75);

				if(0 == m) {
					//no mirror
				} else if(1 == m) {
					//mirror at phi = i * 180 / n
					std::reverse_copy(work.begin()               , work.begin() + iRepeatNum24, work.begin() + iRepeatNum24);
				} else if(2 == m) {
					//mirror at phi = i * 180 / n + 90 / n
					std::reverse_copy(work.begin()               , work.begin() + iRepeatNum14, work.begin() + iRepeatNum14);
					std::reverse_copy(work.begin() + iRepeatNum24, work.begin() + iRepeatNum34, work.begin() + iRepeatNum34);
				} else throw std::runtime_error("invalid m");

				for(size_t j = 1; j < n; j++) {
					const size_t start = (size_t)std::round(repeatNum * j);
					std::copy(work.begin(), work.begin() + iRepeatNum, work.begin() + start);
				}
				square::writeRing(this->dim, i, hemi, work.data());//replace nth ring
			}
		};
		doHemi(nh.data());
		doHemi(sh.data());
	}

	//@brief  : give a master pattern inversion symmetry
	//@note   : doesn't update crystal structure
	template <typename Real>
	void MasterPattern<Real>::makeInvSym() {
		for(size_t j = 0; j < dim; j++) {
			for(size_t i = 0; i < dim; i++) {
				sh[(dim - 1 - j) * dim + (dim - 1 - i)] = nh[dim * j + i];
			}
		}
	}

	//@brief  : copy the equator of a master pattern from the north to south hemisphere
	template <typename Real>
	void MasterPattern<Real>::matchEquator() {
		std::copy(nh. begin(), nh. begin() + dim, sh. begin());
		for(size_t i = 1; i < (size_t)(dim-1); i++) {
			sh[dim* i     ] = nh[dim* i     ];
			sh[dim*(i+1)-1] = nh[dim*(i+1)-1];
		}
		std::copy(nh.rbegin(), nh.rbegin() + dim, sh.rbegin());
	}

	////////////////////////////////////////////////////////////////////////
	//                           MasterSpectra                            //
	////////////////////////////////////////////////////////////////////////

	//@param mp : master pattern
	//@param bw : desired bandwidth
	//@param nrm: should the pattern be normalized (mean 0, stdev 1) before computing the SHT
	//@return   : spectra of square legendre pattern
	//@note     : dim must be odd, bandwidth is dim - 2
	template <typename Real>
	MasterSpectra<Real>::MasterSpectra(MasterPattern<Real> mp, const uint16_t bw, const bool nrm) {
		const uint16_t dimLg = bw + 2 + (bw % 2 == 0 ? 1 : 0);//compute side length of legendre grid for target bandwidth (must be odd)
		if(square::Layout::Legendre != mp.lyt) {//convert grid type if needed
			mp.toLegendre(dimLg);
		}

		const size_t dim = mp.dim;
		if(nrm) {
			//compute pixels sizes relative to average area
			if(dim < dimLg) throw std::runtime_error("insufficient grid resolution for requested bandwidth");
			std::vector<Real> weights(dim * dim);
			std::vector<Real> omega = square::solidAngles<Real>(dim, square::Layout::Legendre);
			for(size_t i = 0; i < weights.size(); i++) weights[i] = omega[square::ringNum(dim, i)];

			//correct weights for double counting of equator
			for(size_t i = 0; i < dim; i++) weights[i                  ] /= 2;
			for(size_t i = 0; i < dim; i++) weights[i * dim            ] /= 2;
			for(size_t i = 0; i < dim; i++) weights[i * dim +  dim - 1 ] /= 2;
			for(size_t i = 0; i < dim; i++) weights[i + dim * (dim - 1)] /= 2;
			const Real totW = std::accumulate(weights.cbegin(), weights.cend(), Real(0));

			//make area weighted average 0
			const Real mean = ( std::inner_product(weights.cbegin(), weights.cend(), mp.nh.cbegin(), Real(0))
			                +   std::inner_product(weights.cbegin(), weights.cend(), mp.sh.cbegin(), Real(0)) ) / totW ;
			std::for_each(mp.nh.begin(), mp.nh.end(), [mean](Real& v){v -= mean;});
			std::for_each(mp.sh.begin(), mp.sh.end(), [mean](Real& v){v -= mean;});

			//make the area weighted standard deviation 1
			std::function<Real(const Real&, const Real&)> multOp = [mean](const Real& v, const Real& w){return v * v * w;};//lambda for inner product multiplication
			const Real sumNh = std::inner_product(mp.nh.cbegin(), mp.nh.cend(), weights.cbegin(), Real(0), std::plus<Real>(), multOp);
			const Real sumSh = std::inner_product(mp.sh.cbegin(), mp.sh.cend(), weights.cbegin(), Real(0), std::plus<Real>(), multOp);
			const Real stdev = std::sqrt((sumNh + sumSh) / (totW * 2));
			std::for_each(mp.nh.begin(), mp.nh.end(), [stdev](Real& v){v /= stdev;});
			std::for_each(mp.sh.begin(), mp.sh.end(), [stdev](Real& v){v /= stdev;});
		}

		//compute SHT
		if(0 == dim % 2) throw std::runtime_error("only odd side lengths are supported");
		MasterData<Real>::phs = mp.phase ();
		MasterData<Real>::kv  = mp.getKv ();
		MasterData<Real>::sig = mp.getSig();
		mBw = bw;
		alm.resize(size_t(mBw) * size_t(mBw));
		square::DiscreteSHT<Real>::Legendre(dim).analyze(mp.nh.data(), mp.sh.data(), alm.data(), mBw, mBw);
		// alm[0] = std::complex<Real>(0);//make sure the mean is truly zero
	}

	//@brief   : repack the spectra into a new maximum bandwidth (0 padded for increasing size)
	//@param bw: bandwidth of new spectra
	//@return  : repacked spectra
	template <typename Real>
	MasterSpectra<Real>& MasterSpectra<Real>::resize(const uint16_t bw) {
		if(bw > mBw) {//zero pad up
			alm.resize(bw * bw, 0);
			for(size_t m = bw - 1; m < bw; m--) {
				std::copy(alm.begin() + m * mBw, alm.begin() + m * mBw + mBw, alm.begin() + m * bw);//move
				std::fill(alm.begin() + m * bw + mBw, alm.begin() + m * bw + bw, std::complex<Real>(0));//zero pad
			}
		} else if(bw < mBw) {//crop down
			for(size_t m = 1; m < bw; m++) std::copy(alm.begin() + m * mBw, alm.begin() + m * mBw + bw, alm.begin() + m * bw);//move
			alm.resize(bw * bw, 0);
		}
		mBw = bw;
		return *this;
	}

	//@brief         : read master patterns from an EMSphInx harmonics file
	//@param fileName: name of EMSphInx spx file to read from
	template <typename Real>
	void MasterSpectra<Real>::read(std::string fileName) {

		//read file into spx structure
		sht::File file;
		std::ifstream is(fileName, std::ios::in | std::ios::binary);
		file.read(is);

		//extract relevant header data
		MasterData<Real>::kv  = (Real) file.header.beamEnergy  ();
		MasterData<Real>::sig = (Real) file.header.primaryAngle();
		mBw = file.harmonics.bw     ();
		MasterData<Real>::phs.pg = xtal::PointGroup(file.material.sgEff  ());
		MasterData<Real>::phs.name.clear();
		MasterData<Real>::pSm.clear();
		if(1 == file.material.numXtal()) {
			for(size_t i = 0; i < 6; i++) MasterData<Real>::phs.lat[i] = (Real) file.material.xtals[0].lat()[i];
		}
		
		//resize data + uncompress
		alm.resize(mBw * mBw, 0);
		file.harmonics.unpackHarm(file.harmonics.alm.data(), alm.data());
	}

}//emsphinx

#endif//_MASTER_H_
