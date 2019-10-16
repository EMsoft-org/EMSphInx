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


#include <vector>

#include "H5Cpp.h"
#include "modality/ebsd/detector.hpp"

struct Dictionary {
	std::vector<float>              eu ;//euler angles
	std::vector<float>              pat;//patterns
	size_t                          pix;//pixels per ebsd pattern
	size_t                          num;//number of orientations in dictionary
	emsphinx::ebsd::Geometry<float> geo;//geometry

	//@brief  : get a pointer to the start of an euler agnle
	//@param i: euler angle number to get pointer to
	//@return : pointer to euler angle start
	float const * getEu (const size_t i) const {return eu .data() + i * 3  ;}

	//@brief  : get a pointer to the start of a pattern
	//@param i: pattern number to get pointer to
	//@return : pointer to pattern start
	float const * getPat(const size_t i) const {return pat.data() + i * pix;}

	//@brief  : get the number of patterns in the dictionary
	//@return : total number of patterns
	size_t numPat() const {return num;}

	//@brief         : read a dictionary
	//@param fileName: name of dictionary file to read
	//@param maxPat  : maximum number of patterns to read (0 for all)
	Dictionary(std::string fileName, const size_t maxPat = 0);

	//@brief    : get detector geometry
	//@param tlt: sample tilt
	//@return   : geometry
	template <typename Real> emsphinx::ebsd::Geometry<Real> getGeom(Real tlt) const;
};

#include <iostream>

//@brief   : write a namelist template file
//@param os: ostream to write namelist template to
//@return  : os
std::ostream& writeDictNml(std::ostream& os) {
	os << " &Placeholder\n";
	os << "!#################################################################\n";
	os << "! Indexing Parameters\n";
	os << "!#################################################################\n";
	os << "\n";
	os << "! spherical harmonic bandwidth to be used (2*bw-1 should be a product of small primes for speed)\n";
	os << "! some reasonable values are: 53, 63, 68, 74, 88, 95, 113, 123, 158\n";
	os << " bw        = 68,\n";
	os << "\n";
	os << "! should newton's method based orientation refinement be used?\n";
	os << " refine    = .TRUE.,\n";
	os << "\n";
	os << "! number of work threads\n";
	os << "! 0 to multithread with an automatic number of threads\n";
	os << "! 1 for serial threading\n";
	os << "! N to multithread with N threads\n";
	os << " nthread   = 0,\n";
	os << "\n";
	os << "! number of patterns to index per work itme (ignored for single threading)\n";
	os << "! should be large enough to make the task significant compared to thread overhead\n";
	os << "! should be small enough to enable enough work items for load balancing\n";
	os << "! 0 to estimate a reasonable value based on speed\n";
	os << " batchsize = 0,\n";
	os << "\n";
	os << "\n";
	os << "!#################################################################\n";
	os << "! Camera Calibration\n";
	os << "!#################################################################\n";
	os << "\n";
	os << "! sample tilt angle from horizontal [degrees]\n";
 	os << " sig      = 70.0,\n";
	os << "\n";
	os << "\n";
	os << "!#################################################################\n";
	os << "! Input Data\n";
	os << "!#################################################################\n";
	os << "\n";
	os << "! input path, empty for current working directory\n";
	os << " ipath       = '',\n";
	os << "\n";
	os << "! master pattern with phases to index (relative to ipath)\n";
	os <<  " masterfile = 'master.h5',\n";
	os << "\n";
	os << "! dictionary file to index (relative to ipath)\n";
	os <<  " dictfile   = 'dict.h5',\n";
	os << "\n";
	os << "! number of patterns to consider (0 for all)\n";
	os <<  " numpat     = 0\n";
	os << " /\n";
	return os;
}

//@brief   : write program usage instructions
//@param os: ostream to write instructions
//@return  : os
std::ostream& writeInstructions(std::ostream& os, char* name) {
	os << "useage: ";
	os << "\tindex using a nml file : " << name << "input.nml\n";
	os << "\tgenerate a template nml: " << name << "-t\n";
	return os;
}

#include <fstream>

#include "modality/ebsd/idx.hpp"
#include "xtal/rotations.hpp"
#include "util/threadpool.hpp"
#include "util/timer.hpp"
#include "util/nml.hpp"
#include "idx/indexer.hpp"

int main(int argc, char *argv[]) {
	////////////////////////////////////////////////////////////////////////
	//                          Parse Arguments                           //
	////////////////////////////////////////////////////////////////////////

	//check argument count
	if(2 != argc) {
		writeInstructions(std::cout, argv[0]);
		return EXIT_FAILURE;
	}

	//check for template request
	const std::string nmlName(argv[1]);
	if(0 == nmlName.compare("-t")) {
		std::ofstream os(std::string(argv[0]) + ".nml");//create programname.nml
		writeDictNml(os);//write template to nml file
		return EXIT_SUCCESS;
	}

	//read nml and parse
	nml::NameList nameList(nmlName);
	typedef double Real;//should we use float, double, or long double for calculations?

	//parse indexing parameters
	const size_t bw        = (size_t) nameList.getInt ("bw"       );//what bandwidth should be used, if 2*bw-1 is product of small primes it is a good candidate for speed (fft is significant fraction of time): 32,38,41,53,63,68,74,88,95,113,123,158
	const bool   refine    = (size_t) nameList.getBool("refine"   );//should newton refinement be used
	const size_t nThread   = (size_t) nameList.getInt ("nthread"  );
	      size_t batchSize = (size_t) nameList.getInt ("batchsize");//number of patterns per work item (should be large enough that the task is significant but small enough that there are enough jobs for load balancing)

	//parse geometry
	const double sampleTilt = nameList.getDouble("sig");

	//parse input files
	std::string ipath, scanName, patName;
	try {        ipath      =         nameList.getString("ipath"     );} catch (...) {}//if ipath isn't found we'll just use cwd
	std::string  masterFile = ipath + nameList.getString("masterfile");
	std::string  dictFile   = ipath + nameList.getString("dictfile"  );
	const size_t maxPat     =         nameList.getInt   ("numpat"    );

	//check for unused inputs
	if(!nameList.fullyParsed()) std::cout << "\nwarning - some namelist parameters weren't used: " << nameList.unusedTokens() << "\n" << std::endl;

	////////////////////////////////////////////////////////////////////////
	//                            Read Inputs                             //
	////////////////////////////////////////////////////////////////////////

	//read master pattern, get symmetry, and read dictionary
	std::vector< emsphinx::MasterSpectra<Real> > phases;
	emsphinx::MasterSpectra<Real> spec(emsphinx::MasterPattern<Real>(masterFile), (uint16_t)bw);//compute SHT of master patterns once
	xtal::PointGroup sym(spec.pointGroup());
	Dictionary dict(dictFile, maxPat);

	//build image processor
	emsphinx::ebsd::Geometry<Real> geom = dict.getGeom<Real>(sampleTilt);
	std::unique_ptr< emsphinx::ebsd::PatternProcessor<Real> > prc(new emsphinx::ebsd::PatternProcessor<Real>());
	prc->setSize(geom.w, geom.h, -1, false, 0);//no circular mask, no background subtraction, no AHE

	//build back projector
	const size_t gridDim = bw + (bw % 2 == 0 ? 3 : 2);
	std::array<Real, 4> quNp = geom.northPoleQuat();
	std::unique_ptr< emsphinx::ebsd::BackProjector<Real> > prj(new emsphinx::ebsd::BackProjector<Real>(geom, gridDim, std::sqrt(Real(2)), quNp.data()));

	//build unnormalized spherical cross correlators
	std::vector< std::unique_ptr< emsphinx::sphere::PhaseCorrelator<Real> > > corrs;
	std::shared_ptr< std::vector< std::complex<Real> > > flm = std::make_shared< std::vector< std::complex<Real> > >(spec.data(), spec.data() + bw * bw);//copy harmonics into shared pointer
	std::unique_ptr< emsphinx::sphere::UnNormalizedCorrelator<Real> > pCorr(new emsphinx::sphere::UnNormalizedCorrelator<Real>((int)bw, flm, spec.mirror(), spec.nFold() ) );//build correlator
	corrs.push_back( std::move( pCorr ) );

	//determine threading parameters and build indexers
	const size_t threadCount = nThread == 0 ? ThreadPool::Concurrency() : nThread;
	if(0 == batchSize) batchSize = emsphinx::Indexer<Real>::BatchEstimate(bw, threadCount, dict.numPat());
	ThreadPool pool(threadCount);//pool

	//build a single indexer for each thread
	std::unique_ptr<emsphinx::Indexer<Real> > idx(new emsphinx::Indexer<Real>((int)bw, prc->clone(), prj->clone(), corrs));//make a single indexer (good for 1 thread)
	std::vector< std::shared_ptr< emsphinx::Indexer<Real> > > indexers;
	indexers.push_back(std::move(idx));//move original indexer onto 
	for(size_t i = 1; i < threadCount; i++) indexers.push_back( std::move(indexers.front()->clone()) );//duplicate n-1 times

	//allocate space to hold indexing result and build work item
	std::vector<std::string> exceptStr(dict.numPat());
	std::vector< emsphinx::Result<Real> > res(dict.numPat());//allocate space for results
	std::function<void(size_t,size_t,size_t)> workItem = [&](const size_t start, const size_t end, const size_t idx){//work function
		for(size_t i = start; i < end; i++) {
			try {
				indexers[idx]->indexImage(dict.getPat(i), &res[i], 1, refine);//index corresponding pattern
			} catch (std::exception& e) {
				exceptStr[i] = e.what();
			}
		}
	};

	////////////////////////////////////////////////////////////////////////
	//                            Do Indexing                             //
	////////////////////////////////////////////////////////////////////////

	//parallel index
	Timer t;
	size_t batches = dict.numPat() / batchSize;//how many batches are needed
	if(batches * batchSize < dict.numPat()) ++batches;//extra batch for leftovers
	for(size_t i = 0; i < batches; i++) {//loop over batches
		const size_t start = i * batchSize;//first pattern
		const size_t end = std::min(start + batchSize, dict.numPat());//last pattern
		pool.schedule(std::bind(workItem, start, end, std::placeholders::_1));//queue indexing
	}
	pool.waitAll();//wait for work to finish

	//print results (+serial index if needed)
	Real eu[3], qu[4], diso[4];
	for(size_t i = 0; i < dict.numPat(); i++) {
		Real * qr = res[i].qu;
		// for(size_t j = 1; j < 4; j++) qr[j] = -qr[j];////////////////////////////!!!!!!TEMPROARY TRY NEGATE
		std::copy(dict.getEu(i), dict.getEu(i) + 3, eu);//get answer
		xtal::eu2qu(eu, qu);//convert to quaternion
		if(0 == res[i].phase) {//result found
			sym.disoQu(qu, qr, diso);//compute disorientation
			Real delta = std::acos(diso[0]) * 360.0 / emsphinx::Constants<Real>::pi;//get disorientation angle
			std::cout << '\t' << i << ": " << delta << '\n';//print error
			if(delta > 2) {//print full result for debugging in case of large error
				std::cout << "\t\t" << qu  [0] << ' ' << qu  [1] << ' ' << qu  [2] << ' ' << qu  [3] << '\n';
				std::cout << "\t\t" << qr  [0] << ' ' << qr  [1] << ' ' << qr  [2] << ' ' << qr  [3] << '\n';
				std::cout << "\t\t" << diso[0] << ' ' << diso[1] << ' ' << diso[2] << ' ' << diso[3] << '\n';
			}
		} else {//error during refinement
			std::cout << '\t' << i << ": " << exceptStr[i] << '\n';
		}
	}

	std::cout << t.poll() << '\n';
	return 0;
}

////////////////////////////////////////////////////////////////////////
//                         Dictionary Details                         //
////////////////////////////////////////////////////////////////////////

//@brief         : read a dictionary
//@param fileName: name of dictionary file to read
Dictionary::Dictionary(std::string fileName, const size_t maxPat) {
	//first open the h5 file and read size
	hsize_t dims[3];//2 or 3 if makeditctionary was true/false
	H5::H5File file = H5::H5File(fileName.c_str(), H5F_ACC_RDONLY);//read only access
	H5::DataSet eulers   = file.openDataSet("EMData/EBSD/EulerAngles");
	H5::DataSet patterns = file.openDataSet("EMData/EBSD/EBSDPatterns");
	patterns.getSpace().getSimpleExtentDims(dims);//read extent in each dimension
	num = maxPat == 0 ? dims[0] : std::min<size_t>(dims[0], maxPat);
	const size_t nDims = patterns.getSpace().getSimpleExtentNdims();
	const bool vectorized = nDims == 2;
	if(!vectorized) {
		if(nDims != 3) throw std::runtime_error("only 2d or 3d patterns are supported");
	}
	pix = vectorized ? dims[1] : dims[1] * dims[2];

	//select hyperslabs to read only first num patterns/euler angles
	hsize_t slabOffsets[3];// hyperslab offset in memory
	hsize_t euSlabDims[2]; // size of the hyperslab in memory
	hsize_t patSlabDims[3]; // size of the hyperslab in memory
	slabOffsets[0] = slabOffsets[1] = slabOffsets[2] = 0;
	 euSlabDims[0] = num;
	 euSlabDims[1] = 3;
	patSlabDims[0] = num;
	patSlabDims[1] = pix;
	if(!vectorized) {
		patSlabDims[1] = dims[1];
		patSlabDims[2] = dims[2];
	}
	H5::DataSpace  euSpace = eulers  .getSpace();
	H5::DataSpace patSpace = patterns.getSpace();

	 euSpace.selectHyperslab(H5S_SELECT_SET,  euSlabDims, slabOffsets);
	patSpace.selectHyperslab(H5S_SELECT_SET, patSlabDims, slabOffsets);

	//allocate memory and read
	eu .resize(3   * num);
	pat.resize(pix * num);
	eulers  .read(eu .data(), H5::PredType::NATIVE_FLOAT,  euSpace,  euSpace);
	patterns.read(pat.data(), H5::PredType::NATIVE_FLOAT, patSpace, patSpace);

	//also parse geometery
	geo.readEMsoft(file.openGroup("NMLparameters/EBSDNameList"));
}

//@brief    : get detector geometry
//@param tlt: sample tilt
//@return   : geometry
template <typename Real>
emsphinx::ebsd::Geometry<Real> Dictionary::getGeom(Real tlt) const {
	emsphinx::ebsd::Geometry<Real> geom;
	geom.sampleTilt   (tlt                        );//sample tilt in degrees
	geom.cameraTilt   (geo.dTlt                   );//camera tilt in degrees
	geom.cameraSize   (geo.w    , geo.h , geo.pX  );///width, height, pixel size in microns
	geom.patternCenter(geo.cX   , geo.cY, geo.sDst);//x/y pattern center + scintilator distance
	geom.maskPattern  (geo.circ                   );
	geom.flipPattern  (geo.flip                   );
	return geom;
}
