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

#ifndef _tsl_h_
#define _tsl_h_

#include <vector>
#include <string>
#include <cstdint>

#ifdef XTAL_USE_H5
	#include "H5Cpp.h"
#endif

namespace tsl {
	struct HKLFamily {
		int32_t h,k,l    ;//hkl plane
		float   intensity;//diffraction intensity
		int8_t  useIdx   ;//use in indexing
		int8_t  showBands;//overlay bands on indexed patterns
	};

	struct Phase {
		size_t                 num   ;//phase number
		std::string            name  ;//material name
		std::string            form  ;//chemical formula
		std::string            info  ;//additional information
		int32_t                sym   ;//tsl symmetry number
		float                  lat[6];//lattice constants (a, b, c, alpha, beta, gamma in angstroms, degrees)
		std::vector<HKLFamily> hklFam;//hkl families
		float                  el[36];//elastic constants (6x6 matrix in row major order)
		std::vector<size_t>    cats  ;//categories (not sure what this is...)
	};

	enum class GridType {Unknown, Square, Hexagonal};//enumeration of grid types

	//@brief     : read a grid type from an input stream
	//@param is  : stream to read from
	//@param grid: grid type to read into
	//@return    : stream read from
	std::istream& operator>>(std::istream& is,       GridType& grid);

	//@brief     : write a grid to an output stream
	//@param os  : stream to write to
	//@param grid: grid type to write from
	//@return    : stream written to
	std::ostream& operator<<(std::ostream& os, const GridType& grid);

	//@brief         : get the type of a file
	//@param fileName: name to parse extension from
	//@return        : parsed extension
	enum class FileType {Unknown, Ang, Osc, Hdf};   //enumeration of file types
	FileType getFileType(std::string fileName);

	struct OrientationMap {
		//header information
		float               pixPerUm             ;
		float               xStar, yStar, zStar  ;//pattern center calibration
		float               sampTlt, camTlt      ;//sample/camera tilt in degrees (only in h5 files)
		float               workingDistance      ;//working distance in mm
		float               xStep   , yStep      ;//pixel size in microns
		int32_t             nColsOdd, nColsEven  ;//width in pixels (same for square grid, alternating rows for hex grid)
		int32_t             nRows                ;//height in pixels
		std::string         operatorName         ;//operator name
		std::string         sampleId             ;//same ID string
		std::string         scanId               ;//scan ID string
		GridType            gridType             ;//square/hex grid
		std::vector<Phase>  phaseList            ;//list of indexed phases (the 'phase' scan data indexes into this array)

		//scan data (all in row major order)
		std::vector<float > eu                   ;//euler angle triples for each pixel
		std::vector<float > x, y                 ;//x/y coordinate of pixel in microns
		std::vector<float > iq                   ;//image quality
		std::vector<float > ci                   ;//confidence index
		std::vector<float > sem                  ;//secondary electron signal
		std::vector<float > fit                  ;//fit
		std::vector<int8_t> phase                ;//phase ID of each pixel (indexes into phaseList)

		//@brief: construct an empty orientation map
		OrientationMap() : pixPerUm(1), 
			xStar(NAN), yStar(NAN), zStar(NAN),
			sampTlt(NAN), camTlt(NAN),
			workingDistance(NAN),
			xStep(NAN), yStep(NAN),
			gridType(GridType::Unknown) {}

		//@brief         : construct an orientation map from a file
		//@param fileName: file to read
		//@param aux     : auxilary information (for hdf5 datasets this is the path to the dataset)
		OrientationMap(std::string fileName, const std::string aux = "") : gridType(GridType::Unknown) {read(fileName, aux);}

		//@brief : check if a file can be ready by this class (based on file extension)
		//@return: true/false if the file type can/cannot be read
		static bool CanRead(std::string fileName);

		//@brief           : allocate space to hold scan data based on grid type and dimensions
		//@param tokenCount: number of arrays to use - eu, x, y, iq, ci and phase are always allocated, sem is only allocated for 9+ tokens and fit for 10+
		void allocate(const size_t tokenCount);

		//@brief         : read scan data from a TSL orientation map file
		//@param fileName: file to read (currently only .ang is supported)
		//@param aux     : auxilary information (for hdf5 datasets this is the path to the dataset)
		void read(std::string fileName, const std::string aux = "");

		//@brief         : write scan data to a TSL orientation map file
		//@param fileName: file to write (currently only .ang is supported)
		void write(std::string fileName);

#ifdef XTAL_USE_H5
		//@brief    : read data from a '.h5' file
		//@param grp: folder to read data from
		void readH5(H5::Group grp);
#endif

		private:
			//@brief         : read data from a '.ang' file
			//@param fileName: name of ang file to read
			//@return        : number of scan points read from file
			void readAng(std::string fileName);

			//@brief   : read an ang header and parse the values
			//@param is: input stream to read the header from
			//@return  : number of tokens (number of data columns)
			size_t readAngHeader(std::istream& is);

			//@brief       : read ang data using an input stream
			//@param is    : input stream set data start
			//@param tokens: number of tokens per point
			//@return      : number of points (rows) parsed
			size_t readAngData(std::istream& is, size_t tokens);

			//@brief         : write data to an '.ang' file
			//@param fileName: name of ang file to write
			void writeAng(std::string fileName);

			//@brief   : write ang header
			//@param os: output stream to write the header to
			void writeAngHeader(std::ostream& os);

			//@brief       : write ang data
			//@param os    : output stream to write data to
			void writeAngData(std::ostream& os);
	};

}

////////////////////////////////////////////////////////////////////////////////
//                           Implementation Details                           //
////////////////////////////////////////////////////////////////////////////////

#include <cstring>
#include <cctype>
#include <sstream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <iterator>

namespace tsl {

	//@brief     : read a grid type from an input stream
	//@param is  : stream to read from
	//@param grid: grid type to read into
	//@return    : stream read from
	std::istream& operator>>(std::istream& is, GridType& grid) {
		std::string name;
		is >> name;
		if     (0 == name.compare("SqrGrid")) grid = GridType::Square;
		else if(0 == name.compare("HexGrid")) grid = GridType::Hexagonal;
		else grid = GridType::Unknown;
		return is;
	}

	//@brief     : write a grid to an output stream
	//@param os  : stream to write to
	//@param grid: grid type to write from
	//@return    : stream written to
	std::ostream& operator<<(std::ostream& os, const GridType& grid) {
		switch(grid) {
			case GridType::Square   : return os << "SqrGrid";
			case GridType::Hexagonal: return os << "HexGrid";
			default: throw std::runtime_error("unknown grid type");
		}
	}

	//@brief         : get the type of a file
	//@param fileName: name to parse extension from
	//@return        : parsed extension
	FileType getFileType(std::string fileName) {
		size_t pos = fileName.find_last_of(".");//find the last '.' in the name
		if(std::string::npos == pos) return FileType::Unknown;//handle files with no extension
		std::string ext = fileName.substr(pos+1);//extract the file extension
		std::transform(ext.begin(), ext.end(), ext.begin(), [](char c){return std::tolower(c);});//convert to lowercase

		if     (0 == ext.compare("ang" )) return FileType::Ang;
		else if(0 == ext.compare("osc" )) return FileType::Osc;
		else if(0 == ext.compare("hdf" )) return FileType::Hdf;
		else if(0 == ext.compare("hdf5")) return FileType::Hdf;
		else if(0 == ext.compare("h5"  )) return FileType::Hdf;
		else return FileType::Unknown;
	}

	//@brief : check if a file can be ready by this class (based on file extension)
	//@return: true/false if the file type can/cannot be read
	bool OrientationMap::CanRead(std::string fileName) {
		switch(getFileType(fileName)) {
			case FileType::Ang    : return true ;
		#ifdef XTAL_USE_H5
			case FileType::Hdf    : return true ;
		#else
			case FileType::Hdf    : return false;
		#endif
			case FileType::Osc    : return false;
			case FileType::Unknown: return false;
		}
		return false;
	}

	//@brief           : allocate space to hold scan data based on grid type and dimensions
	//@param tokenCount: number of arrays to use - eu, x, y, iq, ci and phase are always allocated, sem is only allocated for 9+ tokens and fit for 10+
	void OrientationMap::allocate(const size_t tokenCount) {
		//compute number of pixels based on dimensions and grid type
		size_t totalPoints = 0;
		switch(gridType) {
			case GridType::Square: {
				totalPoints = std::max(nColsOdd, nColsEven) * nRows;
			} break;

			case GridType::Hexagonal: {
				totalPoints = size_t(nRows / 2) * (nColsOdd + nColsEven);
				if(1 == nRows % 2) totalPoints += nColsOdd;
			} break;

			default: throw std::runtime_error("only Square and Hexagonal grid types are supported");
		}

		//allocate arrays (filling new space with zero)
		eu   .resize(3*totalPoints);
		x    .resize(  totalPoints);
		y    .resize(  totalPoints);
		iq   .resize(  totalPoints);
		ci   .resize(  totalPoints);
		phase.resize(  totalPoints);
		if(tokenCount > 8) sem.resize(totalPoints);
		if(tokenCount > 9) fit.resize(totalPoints);
	}

	//@brief         : read scan data from a TSL orientation map file
	//@param fileName: file to read (currently only .ang is supported)
	//@param aux     : auxilary information (for hdf5 datasets this is the path to the dataset)
	void OrientationMap::read(std::string fileName, const std::string aux) {
		//read data from the file
		switch(getFileType(fileName)) {//dispatch the file to the appropraite reader based on the extension
			case FileType::Ang: readAng(fileName); break;
		#ifdef XTAL_USE_H5
			case FileType::Hdf: {
				H5::H5File file(fileName.c_str(), H5F_ACC_RDONLY);
				readH5(file.openGroup(aux));
			}
		#endif
			default: std::runtime_error("unsupported file type (currently only .ang files are supported)");
		}
	}

	//@brief         : write scan data to a TSL orientation map file
	//@param fileName: file to write (currently only .ang is supported)
	void OrientationMap::write(std::string fileName) {
		switch(getFileType(fileName)) {//dispatch the file to the appropraite writer based on the extension
			case FileType::Ang: writeAng(fileName); return;
			default: throw std::runtime_error("unsupported file type (currently only .ang files are supported)");
		}
	}

	//@brief         : read data from a '.ang' file
	//@param fileName: name of ang file to read
	//@return        : number of scan points read from file
	void OrientationMap::readAng(std::string fileName) {
		//parse the header
		std::ifstream is(fileName.c_str());//open file
		if(!is) throw std::runtime_error("ang file " + fileName + " doesn't exist");
		size_t tokenCount = readAngHeader(is);//read header and count number of tokens per point
		allocate(tokenCount);//allocate space

		//read the data
		const size_t pointsRead = readAngData(is, tokenCount);

		//check that enough data was read
		if(pointsRead < iq.size()) {//I'll compare against IQ since it is always present (and easier than comparing against euler angles)
			std::stringstream ss;
			ss << "file ended after reading " << pointsRead << " of " << iq.size() << " data points";
			throw std::runtime_error(ss.str());
		}
	}

	//@brief   : read an ang header and parse the values
	//@param is: input stream to read the header from
	//@return  : number of tokens (number of data columns)
	size_t OrientationMap::readAngHeader(std::istream& is) {
		//flags for which header tokens have been parsed
		bool readPixPerUm        = false;
		bool readXStar           = false, readYStar    = false, readZStar = false;
		bool readWorkingDistance = false;
		bool readXStep           = false, readYStep    = false;
		bool readColsOdd         = false, readColsEven = false;
		bool readRows            = false;
		bool readOperatorName    = false;
		bool readSampleId        = false;
		bool readScanId          = false;
		bool readGridType        = false;

		//flags for which phase subheader tokens have been parsed (set to true values in case there are no phases listed)
		bool   readPhaseSymmetry = true;
		size_t targetFamilies    = 0   ;
		bool   readPhaseMaterial = true;
		bool   readPhaseFormula  = true;
		bool   readPhaseInfo     = true;
		bool   readPhaseLattice  = true;
		bool   readPhaseHkl      = true;
		size_t phaseElasticCount = 6   ;
		bool readPhaseCategories = true;

		//initialize values that don't exist in angs
		sampTlt = camTlt = NAN;

		//now start parsing the header
		std::string line;//buffer to hold single line of header
		std::string token;//current header token
		while('#' == is.peek()) {//all header lines start with #, keep going until we get to data
			std::getline(is, line);//extract entire line from file
			std::istringstream iss(line.data()+1);//skip the '#'
			if(iss >> token) {//get the key word if the line isn't blank
				//get value for appropriate key
				if       (0 == token.compare("TEM_PIXperUM"   )) {iss >> pixPerUm       ; readPixPerUm        = true;
				} else if(0 == token.compare("x-star"         )) {iss >> xStar          ; readXStar           = true;
				} else if(0 == token.compare("y-star"         )) {iss >> yStar          ; readYStar           = true;
				} else if(0 == token.compare("z-star"         )) {iss >> zStar          ; readZStar           = true;

				//these are only in relatively new angs
				} else if(0 == token.compare("SampleTiltAngle"     )) {iss >> sampTlt   ;
				} else if(0 == token.compare("CameraElevationAngle")) {iss >> camTlt    ;
				} else if(0 == token.compare("CameraAzimuthalAngle")) {
				} else if(0 == token.compare("PointGroupID"        )) {

				} else if(0 == token.compare("WorkingDistance")) {iss >> workingDistance; readWorkingDistance = true;
				} else if(0 == token.compare("GRID:"          )) {iss >> gridType       ; readGridType        = true;
				} else if(0 == token.compare("XSTEP:"         )) {iss >> xStep          ; readXStep           = true;
				} else if(0 == token.compare("YSTEP:"         )) {iss >> yStep          ; readYStep           = true;
				} else if(0 == token.compare("NCOLS_ODD:"     )) {iss >> nColsOdd       ; readColsOdd         = true;
				} else if(0 == token.compare("NCOLS_EVEN:"    )) {iss >> nColsEven      ; readColsEven        = true;
				} else if(0 == token.compare("NROWS:"         )) {iss >> nRows          ; readRows            = true;
				} else if(0 == token.compare("OPERATOR:"      )) {iss >> operatorName   ; readOperatorName    = true;
				} else if(0 == token.compare("SAMPLEID:"      )) {iss >> sampleId       ; readSampleId        = true;
				} else if(0 == token.compare("SCANID:"        )) {iss >> scanId         ; readScanId          = true;
				} else if(0 == token.compare("Phase"          )) {
					//check that all attributes for previous phase were read
					std::stringstream ss;
					ss << phaseList.size();
					if(    !readPhaseMaterial ) throw std::runtime_error("ang missing material name for phase "     + ss.str());
					if(    !readPhaseFormula  ) throw std::runtime_error("ang missing formula for phase "           + ss.str());
					if(    !readPhaseInfo     ) throw std::runtime_error("ang missing info for phase "              + ss.str());
					if(    !readPhaseSymmetry ) throw std::runtime_error("ang missing symmetry for phase "          + ss.str());
					if(    !readPhaseLattice  ) throw std::runtime_error("ang missing lattice constants for phase " + ss.str());
					if(    !readPhaseHkl      ) throw std::runtime_error("ang missing hkl families for phase "      + ss.str());
					// if(6 != phaseElasticCount ) throw std::runtime_error("ang missing elastic constants for phase " + ss.str());//elastic constants are optional
					// if(    !readPhaseCategories) throw std::runtime_error("ang missing categories for phase "       + ss.str());//categories are optional
					if(    !phaseList.empty()  ) {
						if(targetFamilies < phaseList.back().hklFam.size())
							throw std::runtime_error("ang missing some hkl families for phase " + ss.str());
					}
					targetFamilies = phaseElasticCount = 0;
					readPhaseSymmetry = readPhaseMaterial = readPhaseFormula = readPhaseInfo = false;
					readPhaseLattice = readPhaseHkl = readPhaseCategories = false;

					//add a new blank phase to the list
					phaseList.resize(phaseList.size() + 1);
					iss >> phaseList.back().num;
				} else if(0 == token.compare("MaterialName"  )) {iss >> phaseList.back().name; readPhaseMaterial = true;
				} else if(0 == token.compare("Formula"       )) {iss >> phaseList.back().form; readPhaseFormula  = true;
				} else if(0 == token.compare("Info"          )) {iss >> phaseList.back().info; readPhaseInfo     = true;
				} else if(0 == token.compare("Symmetry"      )) {iss >> phaseList.back().sym ; readPhaseSymmetry = true;
				} else if(0 == token.compare("NumberFamilies")) {
					iss >> targetFamilies;//read the number of families
					phaseList.back().hklFam.reserve(targetFamilies);//allocate space for the families
					readPhaseHkl = true;
				} else if(0 == token.compare("LatticeConstants")) {
					iss >> phaseList.back().lat[0]
					    >> phaseList.back().lat[1]
					    >> phaseList.back().lat[2]
					    >> phaseList.back().lat[3]
					    >> phaseList.back().lat[4]
					    >> phaseList.back().lat[5];
					readPhaseLattice = true;
				} else if(0 == token.compare("hklFamilies")) {
					phaseList.back().hklFam.resize(phaseList.back().hklFam.size() + 1);//add new family (space was already reserved)
					iss >> phaseList.back().hklFam.back().h
					    >> phaseList.back().hklFam.back().k
					    >> phaseList.back().hklFam.back().l
					    >> phaseList.back().hklFam.back().useIdx
					    >> phaseList.back().hklFam.back().intensity
					    >> phaseList.back().hklFam.back().showBands;
				} else if(0 == token.compare("ElasticConstants")) {
					iss >> phaseList.back().el[6*phaseElasticCount + 0]
					    >> phaseList.back().el[6*phaseElasticCount + 1]
					    >> phaseList.back().el[6*phaseElasticCount + 2]
					    >> phaseList.back().el[6*phaseElasticCount + 3]
					    >> phaseList.back().el[6*phaseElasticCount + 4]
					    >> phaseList.back().el[6*phaseElasticCount + 5];
					++phaseElasticCount;
				} else if(0 == token.compare(0, 10, "Categories")) {//tsl doesn't print space between categories and first number
					//rewind to first character after categories key
					iss.seekg((size_t) iss.tellg() + 10 - token.length());
					std::copy(std::istream_iterator<size_t>(iss), std::istream_iterator<size_t>(), std::back_inserter(phaseList.back().cats));
					readPhaseCategories = true;
				} else {
					if("VERSION"          == token ||
					   "HEADER:"          == token ||
						"COLUMN_COUNT:"   == token ||
						"COLUMN_HEADERS:" == token ||
						"COLUMN_UNITS:"   == token) continue;
					if("NOTES:" == token || "COLUMN_NOTES:" == token) {
						std::string end = "# " + token + " End";
						iss >> token;
						if("Start" == token) {//multi line notes
							while(true) {
								if('#' != is.peek()) throw std::runtime_error("header ended before end of multiline note");
								std::getline(is, line);//extract entire line from file
								if(0 == line.compare(0, end.size(), end)) break;//end of notes
							}
						}
						continue;
					}

					throw std::runtime_error("unknown ang header keyword '" + token + "'");
				}
			}
		}

		//check that all values of final phase were read
		std::stringstream ss;
		ss << phaseList.size();
		if(    !readPhaseMaterial  ) throw std::runtime_error("ang missing material name for phase "     + ss.str());
		if(    !readPhaseFormula   ) throw std::runtime_error("ang missing formula for phase "           + ss.str());
		if(    !readPhaseInfo      ) throw std::runtime_error("ang missing info for phase "              + ss.str());
		if(    !readPhaseSymmetry  ) throw std::runtime_error("ang missing symmetry for phase "          + ss.str());
		if(    !readPhaseLattice   ) throw std::runtime_error("ang missing lattice constants for phase " + ss.str());
		if(    !readPhaseHkl       ) throw std::runtime_error("ang missing hkl families for phase "      + ss.str());
		// if(6 != phaseElasticCount  ) throw std::runtime_error("ang missing elastic constants for phase " + ss.str());//elastic constants are optional
		// if(    !readPhaseCategories) throw std::runtime_error("ang missing categories for phase "        + ss.str());//categories are optional
		if(    !phaseList.empty()  ) {
			if(targetFamilies < phaseList.back().hklFam.size())
				throw std::runtime_error("ang missing some hkl families for phase " + ss.str());
		}

		//make sure all header values were read
		if(!readPixPerUm       ) throw std::runtime_error("missing ang header value TEM_PIXperUM"   );
		if(!readXStar          ) throw std::runtime_error("missing ang header value x-star"         );
		if(!readYStar          ) throw std::runtime_error("missing ang header value y-star"         );
		if(!readZStar          ) throw std::runtime_error("missing ang header value z-star"         );
		if(!readWorkingDistance) throw std::runtime_error("missing ang header value WorkingDistance");
		if(!readGridType       ) throw std::runtime_error("missing ang header value GRID"           );
		if(!readXStep          ) throw std::runtime_error("missing ang header value XSTEP"          );
		if(!readYStep          ) throw std::runtime_error("missing ang header value YSTEP"          );
		if(!readColsOdd        ) throw std::runtime_error("missing ang header value NCOLS_ODD"      );
		if(!readColsEven       ) throw std::runtime_error("missing ang header value NCOLS_EVEN"     );
		if(!readRows           ) throw std::runtime_error("missing ang header value NROWS"          );
		if(!readOperatorName   ) throw std::runtime_error("missing ang header value OPERATOR"       );
		if(!readSampleId       ) throw std::runtime_error("missing ang header value SAMPLEID"       );
		if(!readScanId         ) throw std::runtime_error("missing ang header value SCANID"         );

		//extract first line of data without advancing stream
		const std::streamoff start = is.tellg();//save position of data start
		std::getline(is, line);//copy first data line
		is.seekg(start);//rewind stream to data start

		//get number of tokens
		size_t tokenCount = 0;
		std::istringstream iss(line);
		while(iss >> token) tokenCount++;//count number of values in first data line
		if(tokenCount < 8) {//make sure there are enough columns
			std::stringstream ss;
			ss << "unexpected number of ang values per point (got " << tokenCount << ", expected at least 8)";
			throw std::runtime_error(ss.str());
		}
		return tokenCount;
	}

	//@brief       : read ang data using an input stream
	//@param is    : input stream set data start
	//@param tokens: number of tokens per point
	//@return      : number of points (rows) parsed
	size_t OrientationMap::readAngData(std::istream& is, size_t tokens) {
		char line[512];//most ang files have 128 byte lines including '\n' so this should be plenty
		struct PixData{
			float  eu[3]   ;
			float  x  , y  ;
			float  iq , ci ;
			float  sem, fit;
			int8_t phase   ;

			//@brief    : lexicographic sort
			//@param rhs: other pixel data to compare to
			bool operator<(const PixData& rhs) const {return y < rhs.y ? true : x < rhs.x;}
		};

		//read all point (can be in arbitrary order)
		size_t pointsRead = 0;
		std::vector<PixData> pixelData;
		const size_t totalPoints = iq.size();
		pixelData.resize(totalPoints);
		const bool readSem = tokens > 8;
		const bool readFit = tokens > 9;
		for(size_t i = 0; i < totalPoints; i++) {
			char* data = line;//get pointer to line start
			if(!is.getline(line, sizeof(line))) {//try to get line
				pixelData.resize(i+1);//clip empty points
				break;//get next line
			}
			PixData& p = pixelData[i];

			//parse line
			p.eu[0] =         std::strtof(data, &data    );//parse first euler angle
			p.eu[1] =         std::strtof(data, &data    );//parse second euler angle
			p.eu[2] =         std::strtof(data, &data    );//parse third euler angle
			p.x     =         std::strtof(data, &data    );//parse x
			p.y     =         std::strtof(data, &data    );//parse y
			p.iq    =         std::strtof(data, &data    );//parse image quality
			p.ci    =         std::strtof(data, &data    );//parse confidence index
			p.phase = (int8_t)std::strtol(data, &data, 10);//parse phase (as base 10)
			if(readSem) {//are there 9 or more tokens?
				p.sem     =   std::strtof(data, &data);//parse SE signal
				if(readFit) {//are there 10 or more tokens?
					p.fit =   std::strtof(data, NULL);//parse fit
				}
			}
		}

		//now sort points and extract data
		std::sort(pixelData.begin(), pixelData.end());
		for(size_t i = 0; i < pixelData.size(); i++) {
			PixData& p = pixelData[i];
			eu   [3*i  ] = p.eu[0];
			eu   [3*i+1] = p.eu[1];
			eu   [3*i+2] = p.eu[2];
			x    [i]     = p.x    ;
			y    [i]     = p.y    ;
			iq   [i]     = p.iq   ;
			ci   [i]     = p.ci   ;
			phase[i]     = p.phase;
			if(readSem) {//are there 9 or more tokens?
				sem[i]   = p.sem;
				if(readFit) {//are there 10 or more tokens?
					fit[i] = p.fit;
				}
			}
		}
		return pixelData.size();
	}
	
#ifdef XTAL_USE_H5
	//@brief    : read data from a '.h5' file
	//@param grp: folder to read data from
	void OrientationMap::readH5(H5::Group grp) {
		//first get header and data folders
		H5::Group header = grp.openGroup("EBSD/Header");
		H5::Group data   = grp.openGroup("EBSD/Data"  );

		//define some useful h5 parameters once
		H5::DataSpace ds(H5S_SCALAR);//single element data space
		H5::CompType hklType( sizeof(HKLFamily) );//compound type for hkl family
		hklType.insertMember("H"                    , HOFFSET(HKLFamily, h        ), H5::PredType::NATIVE_INT32);
		hklType.insertMember("K"                    , HOFFSET(HKLFamily, k        ), H5::PredType::NATIVE_INT32);
		hklType.insertMember("L"                    , HOFFSET(HKLFamily, l        ), H5::PredType::NATIVE_INT32);
		hklType.insertMember("Diffraction Intensity", HOFFSET(HKLFamily, intensity), H5::PredType::NATIVE_FLOAT);
		hklType.insertMember("Use in Indexing"      , HOFFSET(HKLFamily, useIdx   ), H5::PredType::NATIVE_INT8 );
		hklType.insertMember("Show bands"           , HOFFSET(HKLFamily, showBands), H5::PredType::NATIVE_INT8 );

		//now parse out header contents
		std::string grid;
		H5::DataSet dOp = header.openDataSet("Operator"       );
		H5::DataSet dSa = header.openDataSet("Sample ID"      );
		H5::DataSet dSc = header.openDataSet("Scan ID"        );
		H5::DataSet dGr = header.openDataSet("Grid Type"      );
		header.openDataSet("Pattern Center Calibration/x-star").read(&xStar          , H5::PredType::NATIVE_FLOAT, ds);
		header.openDataSet("Pattern Center Calibration/y-star").read(&yStar          , H5::PredType::NATIVE_FLOAT, ds);
		header.openDataSet("Pattern Center Calibration/z-star").read(&zStar          , H5::PredType::NATIVE_FLOAT, ds);
		header.openDataSet("Sample Tilt"                      ).read(&sampTlt        , H5::PredType::NATIVE_FLOAT, ds);
		header.openDataSet("Camera Elevation Angle"           ).read(&camTlt         , H5::PredType::NATIVE_FLOAT, ds);
		header.openDataSet("Working Distance"                 ).read(&workingDistance, H5::PredType::NATIVE_FLOAT, ds);
		header.openDataSet("Step X"                           ).read(&xStep          , H5::PredType::NATIVE_FLOAT, ds);
		header.openDataSet("Step Y"                           ).read(&yStep          , H5::PredType::NATIVE_FLOAT, ds);
		header.openDataSet("nColumns"                         ).read(&nColsOdd       , H5::PredType::NATIVE_INT32, ds); nColsEven = nColsOdd;
		header.openDataSet("nRows"                            ).read(&nRows          , H5::PredType::NATIVE_INT32, ds);
		dOp.                                                    read( operatorName   , dOp.getStrType()          , ds);
		dSa.                                                    read( sampleId       , dSa.getStrType()          , ds);
		dSc.                                                    read( scanId         , dSc.getStrType()          , ds);
		dGr.                                                    read( grid           , dGr.getStrType()          , ds); std::stringstream(grid) >> gridType;

		//loop over phases
		H5::Group pGrp = header.openGroup("Phase");
		const size_t numPhase = (size_t)pGrp.getNumObjs();
		if(0 == numPhase) throw std::runtime_error("no phases in tsl HDF file");
		phaseList.resize(numPhase);
		for(size_t i = 0; i < numPhase; i++) {
			Phase& p = phaseList[i];
			std::stringstream ss;
			ss << i + 1;
			H5::Group phs = pGrp.openGroup(ss.str());
			p.num = i+1;
			int32_t famNum = 0;

			H5::DataSet dMn = phs.openDataSet("MaterialName");
			H5::DataSet dFm = phs.openDataSet("Formula"     );
			H5::DataSet dIn = phs.openDataSet("Info"        );
			dMn                                      .read( p.name  , dMn.getStrType()          , ds);
			dFm                                      .read( p.form  , dFm.getStrType()          , ds);
			dIn                                      .read( p.info  , dIn.getStrType()          , ds);
			phs.openDataSet("Symmetry"              ).read(&p.sym   , H5::PredType::NATIVE_INT32, ds);
			phs.openDataSet("Lattice Constant a"    ).read(&p.lat[0], H5::PredType::NATIVE_FLOAT, ds);
			phs.openDataSet("Lattice Constant b"    ).read(&p.lat[1], H5::PredType::NATIVE_FLOAT, ds);
			phs.openDataSet("Lattice Constant c"    ).read(&p.lat[2], H5::PredType::NATIVE_FLOAT, ds);
			phs.openDataSet("Lattice Constant alpha").read(&p.lat[3], H5::PredType::NATIVE_FLOAT, ds);
			phs.openDataSet("Lattice Constant beta" ).read(&p.lat[4], H5::PredType::NATIVE_FLOAT, ds);
			phs.openDataSet("Lattice Constant gamma").read(&p.lat[5], H5::PredType::NATIVE_FLOAT, ds);
			phs.openDataSet("NumberFamilies"        ).read(&famNum  , H5::PredType::NATIVE_INT32, ds);
			p.hklFam.resize(famNum);
			phs.openDataSet("hkl Families"          ).read(p.hklFam.data(), hklType);
		}

		//open scan data
		allocate(10);//allocate all arrays
		const size_t scanPts = iq.size();
		H5::DataSet phi1 = data.openDataSet("Phi1"      );
		H5::DataSet phi  = data.openDataSet("Phi"       );
		H5::DataSet phi2 = data.openDataSet("Phi2"      );
		H5::DataSet dX   = data.openDataSet("X Position");
		H5::DataSet dY   = data.openDataSet("Y Position");
		H5::DataSet dIq  = data.openDataSet("IQ"        );
		H5::DataSet dCi  = data.openDataSet("CI"        );
		H5::DataSet dSem = data.openDataSet("SEM Signal");
		H5::DataSet dFit = data.openDataSet("Fit"       );
		H5::DataSet dPhs = data.openDataSet("Phase"     );

		//sanity check sizes
		hsize_t dSetPts[10] = {0};
		phi1.getSpace().getSimpleExtentDims(&dSetPts[0]);
		phi .getSpace().getSimpleExtentDims(&dSetPts[1]);
		phi2.getSpace().getSimpleExtentDims(&dSetPts[2]);
		dX  .getSpace().getSimpleExtentDims(&dSetPts[3]);
		dY  .getSpace().getSimpleExtentDims(&dSetPts[4]);
		dIq .getSpace().getSimpleExtentDims(&dSetPts[5]);
		dCi .getSpace().getSimpleExtentDims(&dSetPts[6]);
		dSem.getSpace().getSimpleExtentDims(&dSetPts[7]);
		dFit.getSpace().getSimpleExtentDims(&dSetPts[8]);
		dPhs.getSpace().getSimpleExtentDims(&dSetPts[9]);
		for(size_t i = 0; i < 10; i++) if(dSetPts[i] != scanPts) throw std::runtime_error("scan size / H5 dataset mismatch reading tsl H5 file");

		//read scan data
		std::vector<float> work(scanPts);
		phi1.read(work .data(), H5::PredType::NATIVE_FLOAT);
		for(size_t i = 0; i < scanPts; i++) eu[3*i  ] = work[i];
		phi .read(work .data(), H5::PredType::NATIVE_FLOAT);
		for(size_t i = 0; i < scanPts; i++) eu[3*i+1] = work[i];
		phi2.read(work .data(), H5::PredType::NATIVE_FLOAT);
		for(size_t i = 0; i < scanPts; i++) eu[3*i+2] = work[i];
		dX  .read(x    .data(), H5::PredType::NATIVE_FLOAT);
		dY  .read(y    .data(), H5::PredType::NATIVE_FLOAT);
		dIq .read(iq   .data(), H5::PredType::NATIVE_FLOAT);
		dCi .read(ci   .data(), H5::PredType::NATIVE_FLOAT);
		dSem.read(sem  .data(), H5::PredType::NATIVE_FLOAT);
		dFit.read(fit  .data(), H5::PredType::NATIVE_FLOAT);
		dPhs.read(phase.data(), H5::PredType::NATIVE_INT8 );
	}
#endif

	//@brief         : write data to an '.ang' file
	//@param fileName: name of ang file to write
	void OrientationMap::writeAng(std::string fileName) {
		std::ofstream os(fileName.c_str());//open file
		writeAngHeader(os);//write data
		writeAngData(os);//write data
	}

	//@brief   : write ang header
	//@param os: output stream to write the header to
	void OrientationMap::writeAngHeader(std::ostream& os) {
		os << std::fixed << std::setprecision(6);
		os << "# TEM_PIXperUM     " << pixPerUm        << '\n';
		os << "# x-star           " << xStar           << '\n';
		os << "# y-star           " << yStar           << '\n';
		os << "# z-star           " << zStar           << '\n';
		os << "# WorkingDistance  " << workingDistance << '\n';
		os << "#\n";
		for(const Phase& phase : phaseList) {
			os << "# Phase            " << phase.num           << '\n';
			os << "# MaterialName     " << phase.name          << '\n';
			os << "# Formula          " << phase.form          << '\n';
			os << "# Info             " << phase.info          << '\n';
			os << "# Symmetry         " << phase.sym           << '\n';
			os << "# NumberFamilies   " << phase.hklFam.size() << '\n';
			os << std::setprecision(3);
			os << "# LatticeConstants " << phase.lat[0]        << ' ' 
			                            << phase.lat[1]        << ' ' 
			                            << phase.lat[2]        << ' ' 
			                            << phase.lat[3]        << ' ' 
			                            << phase.lat[4]        << ' ' 
			                            << phase.lat[5]        << '\n';
            os << std::setprecision(6);
			for(const HKLFamily& hkl : phase.hklFam) {
				os << "# hklFamilies      " << std::setw(3) <<       hkl.h         << ' '
				                            << std::setw(3) <<       hkl.k         << ' '
				                            << std::setw(3) <<       hkl.l         << ' '
				                                            << (int) hkl.useIdx    << ' '
				                            << std::setw(9) <<       hkl.intensity << ' '
				                                            << (int) hkl.showBands << '\n';
			}
			for(size_t i = 0; i < 6; i++) {
				os << "# ElasticConstants " << phase.el[6*i+0] << ' '
				                            << phase.el[6*i+1] << ' '
				                            << phase.el[6*i+2] << ' '
				                            << phase.el[6*i+3] << ' '
				                            << phase.el[6*i+4] << ' '
				                            << phase.el[6*i+5] << '\n';
			}
			os << "# Categories";
			for(const size_t & c : phase.cats) os << c << ' ';
			os << '\n';
			os << "#\n";
		}
		os << "# GRID:            " << gridType        << '\n';
		os << "# XSTEP:           " << xStep           << '\n';
		os << "# YSTEP:           " << yStep           << '\n';
		os << "# NCOLS_ODD:       " << nColsOdd        << '\n';
		os << "# NCOLS_EVEN:      " << nColsEven       << '\n';
		os << "# NROWS:           " << nRows           << '\n';
		os << "#\n";
		os << "# OPERATOR:        " << operatorName    << '\n';
		os << "#\n";
		os << "# SAMPLEID:        " << sampleId        << '\n';
		os << "#\n";
		os << "# SCANID:          " << scanId          << '\n';
		os << "#\n";
	}

	//@brief       : write ang data
	//@param os    : output stream to write data to
	void OrientationMap::writeAngData(std::ostream& os) {
		const bool writeSem = !sem.empty();
		const bool writeFit = !fit.empty() && writeSem;
		const size_t pts = nRows * nColsOdd;
		if(eu   .size() != pts * 3) throw std::runtime_error("incomplete eu    data for ang");
		if(x    .size() != pts    ) throw std::runtime_error("incomplete x     data for ang");
		if(y    .size() != pts    ) throw std::runtime_error("incomplete y     data for ang");
		if(iq   .size() != pts    ) throw std::runtime_error("incomplete iq    data for ang");
		if(ci   .size() != pts    ) throw std::runtime_error("incomplete ci    data for ang");
		if(phase.size() != pts    ) throw std::runtime_error("incomplete phase data for ang");
		if(writeSem) {
			if(sem.size() != pts) throw std::runtime_error("incomplete sem data for ang");
			if(writeFit) {
				if(fit.size() != pts) throw std::runtime_error("incomplete fit data for ang");
			}
		}

		for(int32_t j = 0; j < nRows; j++) {
			for(int32_t i = 0; i < nColsOdd; i++) {
				const int32_t idx = j * nRows + i;
				os << std::setprecision(5);
				os        << std::setw( 9) <<      eu   [3*idx  ];
				os << ' ' << std::setw( 9) <<      eu   [3*idx+1];
				os << ' ' << std::setw( 9) <<      eu   [3*idx+2];
				os << ' ' << std::setw(12) <<      x    [  idx  ];
				os << ' ' << std::setw(12) <<      y    [  idx  ];
				os << std::setprecision(1);
				os << ' ' << std::setw( 7) <<      iq   [  idx  ];
				os << std::setprecision(3);
				os << ' ' << std::setw( 6) <<      ci   [  idx  ];
				os << ' ' << std::setw( 2) << (int)phase[  idx  ];
				if(writeSem) {
					os << ' ' << std::setw( 6) << sem[idx];
					if(writeFit) os << ' ' << std::setw( 6) << fit[idx];
				}
				os << '\n';
			}
		}


/*
		bool evenRow = true;
		size_t pointsRead = 0;
		size_t completeRowPoints = 0;
		size_t currentCol = nColsEven - 1;
		for(size_t idx = 0; idx < pts; idx++) {
			const size_t i = completeRowPoints + currentCol;//get index of point currently being parsed
			os << std::setprecision(5);
			os        << std::setw( 9) <<      eu   [3*i  ];
			os << ' ' << std::setw( 9) <<      eu   [3*i+1];
			os << ' ' << std::setw( 9) <<      eu   [3*i+2];
			os << ' ' << std::setw(12) <<      x    [  i  ];
			os << ' ' << std::setw(12) <<      y    [  i  ];
			os << std::setprecision(1);
			os << ' ' << std::setw( 7) <<      iq   [  i  ];
			os << std::setprecision(3);
			os << ' ' << std::setw( 6) <<      ci   [  i  ];
			os << ' ' << std::setw( 2) << (int)phase[  i  ];
			if(writeSem) {
				os << ' ' << std::setw( 6) << sem[i];
				if(writeFit) os << ' ' << std::setw( 6) << fit[i];
			}
			os << '\n';
			if(0 == currentCol--) {//decrement current column and check if we've reached a new row
				completeRowPoints += evenRow ? nColsEven : nColsOdd;//update increment to start of current row
				evenRow = !evenRow;//are we currently on an even or odd row?
				currentCol = evenRow ? nColsEven - 1 : nColsOdd - 1;//get number of point in new row
			}
		}
		*/
	}
}

#endif//_tsl_h_
