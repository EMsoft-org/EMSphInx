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

#ifndef _EBSD_NML_
#define _EBSD_NML_

#include <vector>
#include <string>

#include "idx/roi.h"
#include "util/nml.hpp"

namespace emsphinx {

	namespace ebsd {

		//@brief: encapsulation of all the parameters needed for EBSD indexing
		struct Namelist {

			std::string              ipath          ;//input file path
			std::string              patFile        ;//pattern file [hdf, up1, up2, ebsp, data]
			std::string              patName        ;//path to dataset within pattern file (hdf only)
			std::vector<std::string> masterFiles    ;//master pattern files
			std::string              pSymFile       ;//file for psuedosymmetry

			int32_t                  patDims[2]     ;//binned width and height of ebsd patterns
			int32_t                  circRad        ;//circular mask radius
			bool                     gausBckg       ;//should a gaussian background subtraction be applied
			int32_t                  nRegions       ;//AHE tile count

			double                   delta          ;//unbinned pixel size in microns
			std::string              ven            ;//pattern center type
			double                   pctr[3]        ;//pattern center
			double                   thetac         ;//camera tilt in degrees

			int32_t                  scanDims[2]    ;//scan width and height in pixels
			double                   scanSteps[2]   ;//pixel width and height in microns
			std::string              scanFile       ;//scan file to read dimensions from
			std::string              scanName       ;//scan file to read pixel size from
			RoiSelection             roi            ;//indexing region of interest

			int32_t                  bw             ;//spherical harmonic bandwidth to index with
			bool                     normed         ;//should normalized cross correlation be used instead of unnormalized
			bool                     refine         ;//should newton's method refinement be used instead of subpixel interpolation
			int32_t                  nThread        ;
			int32_t                  batchSize      ;

			std::string              opath          ;//output file path
			std::string              dataFile       ;//output HDF file
			std::string              vendorFile     ;//output ang or ctf
			std::string              ipfName        ;//output ipf map
			std::string              qualName       ;//output confidence map (cross correlation)

			std::string              namelistFile   ;//namelist file used for building (actual file contents)

			//@brief: clear inputs files
			void clearInputs() {ipath = patFile = patName = pSymFile = ""; masterFiles.clear();}

			//@brief: clear the pattern file image processing
			void clearImPrc() {patDims[0] = patDims[1] = -1; circRad = -1; gausBckg = false; nRegions = 0;}

			//@brief: clear the pattern center
			void clearPctr() {ven = "EMsoft"; pctr[0] = pctr[1] = pctr[2] = NAN; thetac = delta = NAN;}

			//@brief: clear the scan dimensions
			void clearDims() {scanDims[0] = scanDims[1] = -1; scanSteps[0] = scanSteps[1] = NAN; scanFile = scanName = ""; roi.clear();}

			//@brief: clear indexing parameters
			void clearIdxPrm() {bw = -1; normed = false; refine = false; nThread = 0; batchSize = 0;}

			//@brief: clear outputs
			void clearOutputs() {opath = dataFile = vendorFile = ipfName = qualName = "";}

			//@brief: initialize / reset values with defaults
			void clear();

			//@brief: fill namelist with reasonable defaults
			void defaults();

			//@biref: construct an emty namelist
			Namelist() {clear();}

			//@brief    : parse indexing values from a namelist file
			//@param nml: file to parse (actual contents)
			//@return   : warning string (list of unused namelist values)
			std::string from_string(std::string nml);

			//@brief    : parse indexing values from a namelist file
			//@param nml: file to parse (namelist object)
			//@return   : warning string (list of unused namelist values)
			std::string parse_nml(nml::NameList& nml);

			//@brief    : convert a namelist to a string
			//@param nml: namelist name (for fortran)
			//@return   : namelist file string
			std::string to_string(std::string nml = "EMSphInx") const;

			//@brief     : set pattern file and parse as much information as possible
			//@param file: file to set
			//@param aux : path to h5 dataset if needed
			//@note      : for h5 files searches for a scan that the dataset belongs to
			//@note      : for up1/up2 files searches for a matching .ang
			//@note      : for ebsp files searches for a matching ctf
			void getPatDims(std::string file, std::string aux);

			//@brief    : search for a scan file based on the current scan file name
			//@param pIq: [optional] pointer to location to store image quality map from scan file
			//@param pCi: [optional] pointer to location to indexing quality map from scan file
			//@return   : true if additional information was parsed, false otherwise
			//@note     : may parse pattern center, camera tilt, scan dimensions, and pixel size
			bool findScanFile(std::vector<float>* pIq = NULL, std::vector<float>* pCi = NULL);

			//@brief    : read a scan file with a known file name (scanFile/scanName)
			//@param pIq: [optional] pointer to location to store image quality map from scan file
			//@param pCi: [optional] pointer to location to indexing quality map from scan file
			//@return   : true if additional information was parsed, false otherwise
			//@note     : may parse pattern center, camera tilt, scan dimensions, and pixel size
			bool readScanFile(std::vector<float>* pIq = NULL, std::vector<float>* pCi = NULL);

			//@brief: create the output data file and write header data to it
			void writeFileHeader() const;
		};

	}//ebsd

}//emsphinx

#include <sstream>

#include "modality/ebsd/pattern.hpp"
#include "util/sysnames.hpp"
#include "xtal/orientation_map.hpp"

namespace emsphinx {

	namespace ebsd {
		
		//@brief: initialize / reset values with defaults
		void Namelist::clear() {
			clearInputs ();
			clearImPrc  ();
			clearPctr   ();
			clearDims   ();
			clearIdxPrm ();
			clearOutputs();
		}

		void Namelist::defaults() {
			clear();
			ipath        = ""                        ;
			patFile      = "scan.h5"                 ;
			patName      = "Scan 1/EBSD/Data/Pattern";
			masterFiles  = std::vector<std::string>(1, "master.h5");
			pSymFile     = ""                        ;
			patDims[0]   = 640; patDims[1] = 480     ;
			circRad      = -1                        ;
			gausBckg     = false                     ;
			nRegions     = 10                        ;
			delta        = 50.0                      ;
			ven          = "EMsoft"                  ;
			pctr[0]      = 0.0                       ;
			pctr[1]      = 0.0                       ;
			pctr[2]      = 15000.0                   ;
			thetac       = 10.0                      ;
			scanDims[0]  = 256                       ;
			scanDims[1]  = 256                       ;
			scanSteps[0] = 1.0                       ;
			scanSteps[1] = 1.0                       ;
			roi.clear();
			bw           = 68                        ;
			normed       = true                      ;
			refine       = true                      ;
			nThread      = 0                         ;
			batchSize    = 0                         ;
			opath        = ""                        ;
			dataFile     = "SphInx_Scan.h5"          ;
			vendorFile   = "reindexed.ang"           ;
			ipfName      = "ipf.png"                 ;
			qualName     = "qual.png"                ;
		}

		//@brief    : parse indexing values from a namelist file
		//@param nml: file to parse (actual contents)
		//@return   : warning string (list of unused namelist values)
		std::string Namelist::from_string(std::string nml) {
			//start by wrapping input as stream and parsing with namelist reader
			clear();
			namelistFile = nml;//save a copy of the file
			std::istringstream iss(nml);
			nml::NameList nameList;
			nameList.read(iss); 
			return parse_nml(nameList);
		}

		//@brief    : parse indexing values from a namelist file
		//@param nml: file to parse (namelist object)
		//@return   : warning string (list of unused namelist values)
		std::string Namelist::parse_nml(nml::NameList& nml) {
			//parse inputs first
			try {     ipath       =         nml.getString ("ipath"     );} catch (...) {ipath   .clear();}//if ipath isn't found we'll just use cwd
			try {     pSymFile    =         nml.getString ("psymfile"  );} catch (...) {pSymFile.clear();}//if pSymFile isn't found no operators will be checked
			          masterFiles =         nml.getStrings("masterfile");
			          patFile     = ipath + nml.getString("patfile"   );
			const bool h5Pat = H5::H5File::isHdf5(patFile);
			if(h5Pat) patName     =         nml.getString("patdset"   );
			if(!pSymFile.empty()) patFile = ipath + patFile;
			for(std::string& str : masterFiles) str = ipath + str;

			//parse scan dimensions (before pattern center since this will overwrite values)
			try{
				scanFile = ipath + nml.getString("scandims");//try to get scan dims as a file name
				const bool h5Scn = H5::H5File::isHdf5(patFile);
				if(h5Scn) scanName = nml.getString("scanname");
				if(!readScanFile()) throw std::runtime_error("failed to read scan dimensions from " + scanFile);
				pctr[0] = pctr[1] = pctr[2] = thetac = NAN;
			} catch (std::runtime_error&) {
				//if we didn't get a file name check for scan dimensions
				std::vector<double> dims = nml.getDoubles("scandims");
				if(3 != dims.size() && 4 != dims.size()) throw std::runtime_error("expected a filename or dimensions + resolution for 'scandims' in namelist");

				//sanity check scan dimensions
				if(dims[0] != (uint32_t) dims[0] || dims[1] != (uint32_t) dims[1]) throw std::runtime_error("scan dimensinos must be non-negative integers");
				scanDims [0] = (uint32_t) dims[0];
				scanDims [1] = (uint32_t) dims[1];
				scanSteps[0] = dims[2];
				scanSteps[1] = dims.back();
			}
			roi.from_string(nml.getString("roimask"));

			//parse pattern info
			{
				std::vector<int> dims = nml.getInts("patdims");
				if(2 != dims.size()) throw std::runtime_error("patdims must be 2 elements");
				patDims[0] = dims[0]; patDims[1] = dims[1];
			}
			circRad  =          nml.getInt ("circmask");
			gausBckg =          nml.getBool("gausbckg");
			nRegions = (size_t) nml.getInt ("nregions");//adaptive histogram grid resolution

			//parse pattern center
			delta  = nml.getDouble("delta" );
			thetac = nml.getDouble("thetac");
			ven    = nml.getString ("vendor" );
			{
				std::vector<double> ctr = nml.getDoubles("pctr"   );
				if(3 != ctr.size()) throw std::runtime_error("pctr    must be 3 elements");
				pctr[0] = ctr[0]; pctr[1] = ctr[1]; pctr[2] = ctr[2]; 
			}
			if(!("EMsoft" == ven ||
			     "EDAX"   == ven ||
			     "tsl"    == ven ||
			     "Oxford" == ven ||
			     "Bruker" == ven ||
			     "tsl"    == ven )) throw std::runtime_error("unknown vendor for pattern center `" + ven + "'");

			//parse indexing parameters
			bw        = (size_t) nml.getInt ("bw"       );//what bandwidth should be used, if 2*bw-1 is product of small primes it is a good candidate for speed (fft is significant fraction of time): 32,38,41,53,63,68,74,88,95,113,123,158
			normed    = (size_t) nml.getBool("normed"   );//should normalized or unnormalized cross correlation be used
			refine    = (size_t) nml.getBool("refine"   );//should refinement be used
			nThread   = (size_t) nml.getInt ("nthread"  );
			batchSize = (size_t) nml.getInt ("batchsize");//number of patterns per work item (should be large enough that the task is significant but small enough that there are enough jobs for load balancing)

			//parse outputs
			try {opath      = nml.getString("opath"     );} catch (...) {opath     .clear();}//if ipath isn't found we'll just use cwd
			     dataFile   = nml.getString("datafile"  );
			try {vendorFile = nml.getString("vendorfile");} catch (...) {vendorFile.clear();}//if vendorfile isn't found we won't make an vendor file
			try {ipfName    = nml.getString("ipfmap"    );} catch (...) {ipfName   .clear();}//if ipfName    isn't found we won't make an ipf     map
			try {qualName   = nml.getString("qualmap"   );} catch (...) {qualName  .clear();}//if qualName   isn't found we won't make a  quality map

			//check for unused inputs
			if(!nml.fullyParsed()) return nml.unusedTokens();
			return "";
		}

		//@brief    : convert a namelist to a string
		//@param nml: namelist name (for fortran)
		//@return   : namelist file string
		std::string Namelist::to_string(std::string nml) const {
			std::ostringstream ss;
			ss << " &" << nml << "\n";//this if for the fortran version
			ss << "!#################################################################\n";
			ss << "! Input Files\n";
			ss << "!#################################################################\n";
			ss << "\n";
			ss << "! input path, empty for current working directory\n";
			ss << " ipath      = '" << ipath << "',\n";//ignored for fortran version
			ss << "\n";
			ss << "! raw pattern file (relative to ipath) [can be up1, up2, or hdf5]\n";
			ss << " patfile    = '" << patFile << "',\n";
			ss << "\n";
			ss << "! h5 path of raw pattern  (ignored for non hdf5 patfile)\n";
			ss << " patdset    = '" << patName << "',\n";
			ss << "\n";
			ss << "! master pattern with phases to index (relative to ipath)\n";
			ss << " masterfile = "; for(const std::string& str : masterFiles) ss << '\'' << str << "', "; ss << '\n';
			ss << "\n";
			ss << "! file with list of pseudo symmetric rotations to check (or '' for no psym check)\n";
			ss << " psymfile   = '" << pSymFile << "',\n";
			ss << "\n";
			ss << "\n";
			ss << "!#################################################################\n";
			ss << "! Pattern Processing\n";
			ss << "!#################################################################\n";
			ss << "\n";
			ss << "! number of CCD pixels along x and y\n";
			ss << " patdims    = " << patDims[0] << ", " << patDims[1] << ",\n";
			ss << "\n";
			ss << "! should a circular mask be applied (-1 for no mask, 0 for largest inscribed circle, >0 to specify radius in pixels)\n";
			ss << " circmask   = " << circRad << ",\n";
			ss << "\n";
			ss << "! should a 2D gaussian background be subtracted\n";
			ss << " gausbckg   = ." << (gausBckg ? "TRUE" : "FALSE") << ".,\n";
			ss << "\n";
			ss << "! how many regions should be used for adaptive histogram equalization (0 for no AHE)\n";
			ss << " nregions   = " << nRegions << ",\n";
			ss << "\n";
			ss << "\n";
			ss << "!#################################################################\n";
			ss << "! Camera Calibration\n";
			ss << "!#################################################################\n";
			ss << "\n";
			ss << "! CCD pixel size on the scintillator surface [microns]\n";
			ss << " delta      = " << delta << ",\n";
			ss << "\n";
			ss << "! pattern center coordinates and vendor\n";
			ss << "! vendor must be one of the following:\n";
			ss << "!   EMsoft, EDAX, TSL, Oxford, Bruker\n";
			ss << "! with pctr interpreted accordingly:\n";
			ss << "!   EMsoft   - pcx (pixels), pcy (pixels), scintillator distance (microns)\n";
			ss << "!   EDAX/TSL - x*, y*, z*\n";
			ss << "!   Oxford   - x*, y*, z*\n";
			ss << "!   Bruker   - x*, y*, z*\n";
			ss << "! note that vendors use different x*, y*, and z* : https://doi.org/10.1007/s40192-019-00137-4\n";
			ss << " pctr       = " << pctr[0] << ", " << pctr[1] << ", " << pctr[2] << ",\n";
			ss << " vendor     = '" << ven << "',\n";
			ss << "\n";
			ss << "! tilt angle of the camera (positive below horizontal, [degrees]\n";
			ss << " thetac     = " << thetac << ",\n";
			ss << "\n";
			ss << "\n";
			ss << "!#################################################################\n";
			ss << "! Scan Information\n";
			ss << "!#################################################################\n";
			ss << "\n";
			ss << "! dimensions of scan to index and pixel size\n";
			ss << "! x, y, step   for an x by y scan with square pixels of 'step' microns\n";
			ss << "! x, y, sx, sy for an x by y scan with rectangular pixels of 'sx' by 'sy' microns\n";
			ss << "! string to read dimensions from a scan file (*.ang, *.ctf, or *.h5)\n";
			ss << " scandims   = " << scanDims[0] << ", " << scanDims[1] << ", " << scanSteps[0] << ", " << scanSteps[1] << ",\n";//fortran version can't be string
			ss << "\n";
			ss << "! h5 path of scan data folder if scandims is an h5 file (ignored otherwise)\n";
			ss << " scanname   = '',\n";//N/A for fortran version
			ss << "\n";
			ss << "! region of interest for indexing\n";
			ss << "! 0 (or omitted) to index the entire scan\n";
			ss << "! x0, y0, dx, dy for a (dx, dy) rectangular starting at pixel (x0, y0)\n";
			ss << "! string for an ROI mask file\n";
			ss << " roimask    = '" << roi.to_string() << "',\n";//fortran version can't be string
			ss << "\n";
			ss << "!#################################################################\n";
			ss << "! Indexing Parameters\n";
			ss << "!#################################################################\n";
			ss << "\n";
			ss << "! spherical harmonic bandwidth to be used (2*bw-1 should be a product of small primes for speed)\n";
			ss << "! some reasonable values are: 53, 63, 68, 74, 88, 95, 113, 122, 123, 158, 172, 188, 203, 221, 263, 284, 313\n";
			ss << "! a nice range for parameter studies is 53, 68, 88, 113, 158, 203, 263, 338 (~a factor of 1.3 between each)\n";
			ss << "! any value is now pretty fast since the transform is zero padded to the nearest fast FFT size\n";
			ss << " bw         = " << bw << ",\n";
			ss << "\n";
			ss << "! should normalized / unnormalized spherical cross correlation be used?\n";
			ss << "! normalization is more robust for (esp. for lower symmetries) but is slower\n";
			ss << " normed     = ." << (normed ? "TRUE" : "FALSE") << ".,\n";
			ss << "\n";
			ss << "! should newton's method orientation refinement be used?\n";
			ss << "! normalization is more robust for (esp. for lower symmetries) but is slower\n";
			ss << " refine     = ." << (refine ? "TRUE" : "FALSE") << ".,\n";
			ss << "\n";
			ss << "! number of work threads\n";
			ss << "! 0 (or omitted) to multithread with an automatic number of threads\n";
			ss << "! 1 for serial threading\n";
			ss << "! N to multithread with N threads\n";
			ss << " nthread    = " << nThread<< ",\n";
			ss << "\n";
			ss << "! number of patterns to index per work itme (ignored for single threading)\n";
			ss << "! should be large enough to make the task significant compared to thread overhead\n";
			ss << "! should be small enough to enable enough work items for load balancing\n";
			ss << "! should be small enough so nthread * batchsize patterns can be held in memory\n";
			ss << "! 0 (or omitted) to estimate a reasonable value based on speed\n";
			ss << " batchsize  = " << batchSize << ",\n";
			ss << "\n";
			ss << "\n";
			ss << "!#################################################################\n";
			ss << "! Output Files\n";
			ss << "!#################################################################\n";
			ss << "\n";
			ss << "! output path, empty for current working directory\n";
			ss << " opath      = '" << opath << "',\n";//ignored for fortran version
			ss << "\n";
			ss << "! output orientation map name relative to opath [must be hdf5 type]\n";
			ss << " datafile   = '" << dataFile << "',\n"; //should be hdf
			ss << "\n";
			ss << "! output orientation map name relative to opath (or omitted for no vendor output) [can be ang or ctf]\n";
			ss << " vendorfile = '" << vendorFile << "',\n"; //should be hdf
			ss << "\n";
			ss << "! output ipf map with {0,0,1} reference direction (or omitted for no ipf map) [must be png]\n";
			ss << " ipfmap     = '" << ipfName << "',\n";
			ss << "\n";
			ss << "! output quality map with (or omitted for no quality map) [must be png]\n";
			ss << " qualmap    = '" << qualName << "'\n";
			ss << " /\n";
			return ss.str();
		}

		//@brief     : set pattern file and parse as much information as possible
		//@param file: file to set
		//@param aux : path to h5 dataset if needed
		//@note      : for h5 files searches for a scan that the dataset belongs to
		//@note      : for up1/up2 files searches for a matching .ang
		//@note      : for ebsp files searches for a matching ctf
		void Namelist::getPatDims(std::string file, std::string aux) {
			//most pattern files have some basic info
			clear();
			int w, h;
			uint64_t num;
			PatternFile::Bits bits;
			PatternFile::GetFileDims(file, w, h, bits, num, aux);
			patFile = file;
			patName = aux;
			patDims[0] = w;//may be -1 for *.data
			patDims[1] = h;//may be -1 for *.data
		}

		//@brief    : search for a scan file based on the current scan file name
		//@param pIq: [optional] pointer to location to store image quality map from scan file
		//@param pCi: [optional] pointer to location to indexing quality map from scan file
		//@return   : true if additional information was parsed, false otherwise
		//@note     : may parse pattern center, camera tilt, scan dimensions, and pixel size
		bool Namelist::findScanFile(std::vector<float>* pIq, std::vector<float>* pCi) {
			//start by getting file extension
			std::string ext = "";
			std::string name = patFile;//get a copy of the pattern file
			size_t pos = name.find_last_of(".");//find the last '.' in the file
			if(std::string::npos != pos) {
				ext = name.substr(pos+1);//extract the file extension
				std::transform(ext.begin(), ext.end(), ext.begin(), [](char c){return std::tolower(c);});//convert to lowercase
			}

			//next look for vendor file
			bool tryRead = false;
			std::string aux = patName;
			if("up1" == ext || "up2" == ext) {//edax
				ven = "EDAX";
				name[pos+1] = 'a'; name[pos+2] = 'n'; name[pos+3] = 'g';
				tryRead = fileExists(name);
			} else if(".ebsp" == ext) {//oxford
				ven = "Oxford";
				name[pos+1] = 'c'; name[pos+2] = 't'; name[pos+3] = 'f'; name.pop_back();
				tryRead = fileExists(name);
			} else if("h5" == ext || "hdf" == ext || "hdf5" == ext) {//any
				ven = PatternFile::GetVendor(name);
				if(ven.empty()) {
					ven = "EMsoft";
				} else {//we found a vendor string, hopefully this is an orientation map file
					std::string dataPath = aux     .substr(0, aux     .find_last_of("/\\"));
					std::string ebsdPath = dataPath.substr(0, dataPath.find_last_of("/\\"));
					std::string scanPath = ebsdPath.substr(0, ebsdPath.find_last_of("/\\"));//this is the folder that the scan would be in
					if(dataPath != aux && !dataPath.empty()) {//folder structure was as expected (could make this a bit more strict)
						aux = scanPath + '/';
						tryRead = true;
					}
				}
			}

			//try to read vendor file if found
			if(tryRead) {
				scanFile = name;
				scanName = aux ;
				if(!readScanFile(pIq, pCi)) {
					scanFile = scanName = "";//file was bad
				} else {
					return true;
				}
			}
			return false;
		}

		//@brief    : read a scan file with a known file name (scanFile/scanName)
		//@param pIq: [optional] pointer to location to store image quality map from scan file
		//@param pCi: [optional] pointer to location to indexing quality map from scan file
		//@return   : true if additional information was parsed, false otherwise
		//@note     : may parse pattern center, camera tilt, scan dimensions, and pixel size
		bool Namelist::readScanFile(std::vector<float>* pIq, std::vector<float>* pCi) {
			try {
				xtal::OrientationMap<float> om(scanFile, scanName);//this will throw on errors
				//if we made it to here the file had everything we were after
				if(NULL != pIq) pIq->swap(om.imQual);
				if(NULL != pCi) pCi->swap(om.metric);
				pctr     [0] = om.calib.xStar ;
				pctr     [1] = om.calib.yStar ;
				pctr     [2] = om.calib.zStar ;
				thetac       = om.calib.cTlt  ;
				scanDims [0] = om      .width ;
				scanDims [1] = om      .height;
				scanSteps[0] = om      .xStep ;
				scanSteps[1] = om      .yStep ;
				return true;
			} catch (...) {
			}
			
			//if we made it to here we were unable to parse the values of interest
			pctr[0]         =  0      ;//x pattern center (emsoft)
			pctr[1]         =  0      ;//y pattern center (emsoft)
			pctr[2]         =  0      ;//scintillator distance in microns (emsoft)
			thetac          =  0      ;//camera tilt in degrees
			scanDims[0]     = -1      ;//scan width  in pixels
			scanDims[1]     = -1      ;//scan height in pixels
			scanSteps[0]    =  0      ;//pixel width  in microns
			scanSteps[1]    =  0      ;//pixel height in microns
			return false;
		}

		//@brief: create the output data file and write header data to it
		void Namelist::writeFileHeader() const {
			//reconstruct namelist
			std::string str = namelistFile.empty() ? to_string() : namelistFile;//original string if parsed from file, to_string if built programatically
			std::istringstream iss(str);
			nml::NameList nameList;
			nameList.read(iss);
			ebsd::Namelist dummy;
			dummy.parse_nml(nameList);//read files so hdf5 writes used/unused fields correctly

			//open the file and write top level values
			//this is important to do now since we don't want a failure at the end if e.g. there aren't write privileges
			H5::H5File file(opath + dataFile, H5F_ACC_TRUNC);
			std::string manufact = "EMSphInx", version = emsphinx::Version;
			file.createDataSet("Manufacturer", H5::StrType(0, H5T_VARIABLE), H5::DataSpace(H5S_SCALAR)).write(manufact, H5::StrType(0, H5T_VARIABLE));
			file.createDataSet("Version"     , H5::StrType(0, H5T_VARIABLE), H5::DataSpace(H5S_SCALAR)).write(version , H5::StrType(0, H5T_VARIABLE));

			//write nml file (line by line + parsed values)
			std::string programName("IndexEBSD");
			nameList.writeFile(file.createGroup("NMLfiles"), programName);
			nameList.writeParameters(file.createGroup("NMLparameters").createGroup(programName));//copy parsed values

			//build program, date, host, and user strings string
			std::stringstream ss;
			time_t tm = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
			ss << std::put_time(std::localtime(&tm), "%a %b %d %Y");
			programName += " (index_ebsd.cpp)";
			std::string hostName = getComputerName();
			std::string userName = getUserName();

			//write header data
			H5::Group header = file.createGroup("EMheader");
			header.createDataSet("Date"       , H5::StrType(0, H5T_VARIABLE), H5::DataSpace(H5S_SCALAR)).write(ss.str()   , H5::StrType(0, H5T_VARIABLE));
			header.createDataSet("HostName"   , H5::StrType(0, H5T_VARIABLE), H5::DataSpace(H5S_SCALAR)).write(hostName   , H5::StrType(0, H5T_VARIABLE));
			header.createDataSet("UserName"   , H5::StrType(0, H5T_VARIABLE), H5::DataSpace(H5S_SCALAR)).write(userName   , H5::StrType(0, H5T_VARIABLE));
			header.createDataSet("ProgramName", H5::StrType(0, H5T_VARIABLE), H5::DataSpace(H5S_SCALAR)).write(programName, H5::StrType(0, H5T_VARIABLE));
			header.createDataSet("Version"    , H5::StrType(0, H5T_VARIABLE), H5::DataSpace(H5S_SCALAR)).write(version    , H5::StrType(0, H5T_VARIABLE));
		}

	}//ebsd

}//emsphinx


#endif//_EBSD_NML_
