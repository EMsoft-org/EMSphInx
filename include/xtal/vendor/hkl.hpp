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

#ifndef _HKL_H_
#define _HKL_H_

#include <vector>
#include <string>

namespace hkl {
	struct Phase {
		float       lat[6];//a, b, c, alpha, beta, gamma (in angstroms, degrees)
		std::string name  ;
		int32_t     laue  ;//laue group number
		int32_t     space ;//space group number
	};

	//@brief         : get the type of a file
	//@param fileName: name to parse extension from
	//@return        : parsed extension
	enum class FileType {Unknown, Ctf, Hdf};   //enumeration of file types
	FileType getFileType(std::string fileName);

	struct OrientationMap {
		std::string          project               ;//project name
		std::string          author                ;//author name
		std::string          jobMode               ;//job mode
		int32_t              xCells, yCells, zCells;//number of voxels in each direction
		float                xStep , yStep , zStep ;//voxel resolution in microns
		float                acqE1 , acqE2 , acqE3 ;
		std::string          euStr                 ;//euler angle string
		float                mag                   ;//magnification
		float                cover                 ;//coverage
		float                dev                   ;//device
		float                kv                    ;//accelerating voltage 
		float                angle                 ;//tilt angle
		float                axis                  ;//tilt axis
		std::vector<Phase>   phaseList             ;//list of phases

		std::vector<uint8_t> phase                 ;//phase
		std::vector<float  > x, y                  ;//x/y coordinate of pixel in microns
		std::vector<uint8_t> bands                 ;//bands
		std::vector<float  > err                   ;//error measurement
		std::vector<float  > eu                    ;//euler angles
		std::vector<float  > mad                   ;//MAD
		std::vector<float >  bc                    ;//BC
		std::vector<float >  bs                    ;//BS

		//@brief: construct an empty orientation map
		OrientationMap(): zCells(1) {}

		//@brief         : construct an orientation map from a file
		//@param fileName: file to read (currently only .ctf is supported)
		OrientationMap(std::string fileName) {read(fileName);}

		//@brief : check if a file can be ready by this class (based on file extension)
		//@return: true/false if the file type can/cannot be read
		static bool CanRead(std::string fileName);

		//@brief     : allocate space to hold scan data based on grid type and dimensions
		//@param cols: bitmask of columns to allocate
		enum ColumnType : uint16_t {
			CTF_Phase   = 0x0001,
			CTF_X       = 0x0002,
			CTF_Y       = 0x0004,
			CTF_Bands   = 0x0008,
			CTF_Error   = 0x0010,
			CTF_Euler1  = 0x0020,
			CTF_Euler2  = 0x0040,
			CTF_Euler3  = 0x0080,
			CTF_MAD     = 0x0100,
			CTF_BC      = 0x0200,
			CTF_BS      = 0x0400,
			CTF_ALL     = 0x07FF,//all field
			CTF_Euler   = 0x00E0,//CTF_Euler1 | CTF_Euler2 | CTF_Euler3
			CTF_Min     = 0x00E1 //minimal value for a valid CTF file = CTF_Phase | CTF_Euler
		};
		void allocate(const uint16_t cols);

		//@brief         : read scan data from an HKL orientation map file
		//@param fileName: file to read (currently only .ctf is supported)
		void read(std::string fileName);

		//@brief         : write scan data to an HKL orientation map file
		//@param fileName: file to write (currently only .ctf is supported)
		void write(std::string fileName);

		private:
			//@brief         : read data from a '.ctf' file
			//@param fileName: name of ctf file to read
			//@return        : number of scan points read from file
			void readCtf(std::string fileName);

			//@brief   : read an ctf header and parse the values
			//@param is: input stream to read the header from
			//@return  : colunns detected
			std::vector<uint16_t> readCtfHeader(std::istream& is);

			//@brief     : read ctf data using an input stream
			//@param is  : input stream set data start
			//@param cols: column headers
			//@return    : number of points (rows) parsed
			size_t readCtfData(std::istream& is, const std::vector<uint16_t>& cols);

			//@brief         : read ctf data using a memory map
			//@param fileName: name of file to read
			//@param offset  : offset to data start (in bytes)
			//@param cols    : column headers
			//@return        : number of points (rows) parsed
			// size_t readCtfDataMemMap(std::string fileName, std::streamoff offset, const std::vector<uint16_t>& cols);

			//@brief         : write data to a '.ctf' file
			//@param fileName: name of ctf file to write
			void writeCtf(std::string fileName);

			//@brief   : write ctf header
			//@param os: output stream to write the header to
			void writeCtfHeader(std::ostream& os);

			//@brief       : write ctf data
			//@param os    : output stream to write data to
			void writeCtfData(std::ostream& os);
	};
}

////////////////////////////////////////////////////////////////////////////////
//                           Implementation Details                           //
////////////////////////////////////////////////////////////////////////////////

#include <fstream>
#include <sstream>
#include <iomanip>

namespace hkl {
	//@brief         : get the type of a file
	//@param fileName: name to parse extension from
	//@return        : parsed extension
	FileType getFileType(std::string fileName) {
		size_t pos = fileName.find_last_of(".");//find the last '.' in the name
		if(std::string::npos == pos) return FileType::Unknown;//handle files with no extension
		std::string ext = fileName.substr(pos+1);//extract the file extension
		std::transform(ext.begin(), ext.end(), ext.begin(), [](char c){return std::tolower(c);});//convert to lowercase

		if     (0 == ext.compare("ctf" )) return FileType::Ctf;
		else if(0 == ext.compare("hdf" )) return FileType::Hdf;
		else if(0 == ext.compare("hdf5")) return FileType::Hdf;
		else if(0 == ext.compare("h5"  )) return FileType::Hdf;
		else return FileType::Unknown;
	}

	//@brief : check if a file can be ready by this class (based on file extension)
	//@return: true/false if the file type can/cannot be read
	bool OrientationMap::CanRead(std::string fileName) {
		switch(getFileType(fileName)) {
			case FileType::Ctf    : return true ;
			case FileType::Hdf    : return false;
			case FileType::Unknown: return false;
		}
		return false;
	}

	//@brief     : allocate space to hold scan data based on grid type and dimensions
	//@param cols: bitmask of columns to allocate
	void OrientationMap::allocate(const uint16_t cols) {
		//compute number of pixels based on dimensions and grid type
		size_t totalPoints = xCells * yCells;
		if(1 != zCells) throw std::runtime_error("only 2D scans are currently supported");

		if(!(cols & CTF_Min)) throw std::runtime_error("ctf must have at least phase + orientation");

		//allocate arrays (filling new space with zero)
		const bool hasEu = (cols & CTF_Euler1) || (cols & CTF_Euler2) || (cols & CTF_Euler3);
		if(hasEu && !(cols & CTF_Euler)) throw std::runtime_error("cannot parse partial euler angles");
		if(cols & CTF_Phase) phase.resize(  totalPoints);
		if(cols & CTF_X    ) x    .resize(  totalPoints);
		if(cols & CTF_Y    ) y    .resize(  totalPoints);
		if(cols & CTF_Bands) bands.resize(  totalPoints);
		if(cols & CTF_Error) err  .resize(  totalPoints);
		if(cols & CTF_Euler) eu   .resize(3*totalPoints);
		if(cols & CTF_MAD  ) mad  .resize(  totalPoints);
		if(cols & CTF_BC   ) bc   .resize(  totalPoints);
		if(cols & CTF_BS   ) bs   .resize(  totalPoints);
	}

	//@brief         : read scan data from an HKL orientation map file
	//@param fileName: file to read (currently only .ctf is supported)
	void OrientationMap::read(std::string fileName) {
		//read data from the file
		switch(getFileType(fileName)) {//dispatch the file to the appropraite reader based on the extension
			case FileType::Ctf: readCtf(fileName); break;
			default: std::runtime_error("unsupported file type (currently only .ctf files are supported)");
		}
	}

	//@brief         : write scan data to an HKL orientation map file
	//@param fileName: file to write (currently only .ctf is supported)
	void OrientationMap::write(std::string fileName) {
		switch(getFileType(fileName)) {//dispatch the file to the appropraite writer based on the extension
			case FileType::Ctf: writeCtf(fileName); return;
			default: throw std::runtime_error("unsupported file type (currently only .ctf files are supported)");
		}
	}

	//@brief         : read data from a '.ctf' file
	//@param fileName: name of ctf file to read
	//@return        : number of scan points read from file
	void OrientationMap::readCtf(std::string fileName) {
		//parse the header
		std::ifstream is(fileName.c_str());//open file
		if(!is) throw std::runtime_error("ctf file " + fileName + " doesn't exist");
		std::vector<uint16_t> cols = readCtfHeader(is);//read header and count number of tokens per point
		if(0 != jobMode.compare("Grid")) throw std::runtime_error("only grid job mode is supported (got '" + jobMode + "')");

		//allocate space
		uint16_t c = 0;
		for(const uint16_t& i : cols) c |= i;
		allocate(c);
		//read the data
		size_t pointsRead = 0;//nothing has been read
		static const bool UseMemMap = false;
		if(UseMemMap) {
			std::streamoff offset = is.tellg();//get offset to data start
			is.close();//close the file
			throw std::runtime_error("not yet implemented");
		} else {
			pointsRead = readCtfData(is, cols);
		}

		//check that enough data was read
		if(pointsRead < phase.size()) {//I'll compare against phase since it is always present (and easier than comparing against euler angles)
			std::stringstream ss;
			ss << "file ended after reading " << pointsRead << " of " << phase.size() << " data points";
			throw std::runtime_error(ss.str());
		}
	}

	//@brief   : read an ctf header and parse the values
	//@param is: input stream to read the header from
	//@return  : colunns detected
	std::vector<uint16_t> OrientationMap::readCtfHeader(std::istream& is) {
		std::string line;
		std::string token;
		bool readProject   = false, readAuthor    = false, readJobMode  = false;
		bool readXCells    = false, readYCells    = false, readZCells   = false;
		bool readXStep     = false, readYStep     = false, readZStep    = false;
		bool readacqE1     = false, readacqE2     = false, readacqE3    = false;
		bool readEulerLine = false;
		bool readMag       = false, readCoverage  = false, readDevice   = false;
		bool readKv        = false, readTiltAngle = false, readTiltAxis = false;
		std::vector<uint16_t> cols;

		//check for magic string
		std::getline(is, line);//extract first line from file
		if(0 != std::string(line).compare(0, 17, "Channel Text File")) throw std::runtime_error("not a valid ctf file (bad header)");//need to compare this way to handle \r\n from windows files on a mac

		//parse header
		while(cols.empty()) {
			//extract entire line from file and copy to stringstream
			if(!std::getline(is, line)) throw std::runtime_error("not a valid ctf file (no column headers)");//if we failed to extract a line we're in trouble
			std::istringstream iss(line);

			//get the key word if the line isn't blank
			if(iss >> token) {
				//get value for appropriate key
				if(0 == token.compare("Phase" ) || // check for column headers first so we find Euler(1/2/3) before "Euler angles refer to..."
				   0 == token.compare("X"     ) || 
				   0 == token.compare("Y"     ) || 
				   0 == token.compare("Bands" ) || 
				   0 == token.compare("Error" ) || 
				   0 == token.compare("Euler1") || 
				   0 == token.compare("Euler2") || 
				   0 == token.compare("Euler3") || 
				   0 == token.compare("MAD"   ) || 
				   0 == token.compare("BC"    ) || 
				   0 == token.compare("BS"    ) ) {
					//parse the list of tokens
					bool euPresent[3] = {false, false, false};
					do {
						if     (0 == token.compare("Phase" )) cols.push_back(CTF_Phase );
						else if(0 == token.compare("X"     )) cols.push_back(CTF_X     );
						else if(0 == token.compare("Y"     )) cols.push_back(CTF_Y     );
						else if(0 == token.compare("Bands" )) cols.push_back(CTF_Bands );
						else if(0 == token.compare("Error" )) cols.push_back(CTF_Error );
						else if(0 == token.compare("Euler1")){cols.push_back(CTF_Euler1); euPresent[0] = true;}
						else if(0 == token.compare("Euler2")){cols.push_back(CTF_Euler2); euPresent[1] = true;}
						else if(0 == token.compare("Euler3")){cols.push_back(CTF_Euler3); euPresent[2] = true;}
						else if(0 == token.compare("MAD"   )) cols.push_back(CTF_MAD   );
						else if(0 == token.compare("BC"    )) cols.push_back(CTF_BC    );
						else if(0 == token.compare("BS"    )) cols.push_back(CTF_BS    );
						else throw std::runtime_error("unknown ctf column header: " + token);
					} while(iss >> token);
					bool allEulers = euPresent[0] && euPresent[1] && euPresent[2];
					bool anyEulers = euPresent[0] || euPresent[1] || euPresent[2];
					if(anyEulers && !allEulers) throw std::runtime_error("only some euler angles present in column headers");
				} else if(0 == token.compare("Prj"      )) { std::getline(iss >> std::ws, project); readProject   = true;
				} else if(0 == token.compare("Author"   )) { std::getline(iss >> std::ws, author ); readAuthor    = true;
				} else if(0 == token.compare("JobMode"  )) { std::getline(iss >> std::ws, jobMode); readJobMode   = true;
				} else if(0 == token.compare("XCells"   )) {              iss >> xCells           ; readXCells    = true;
				} else if(0 == token.compare("YCells"   )) {              iss >> yCells           ; readYCells    = true;
				} else if(0 == token.compare("ZCells"   )) {              iss >> zCells           ; readZCells    = true;
				} else if(0 == token.compare("XStep"    )) {              iss >> xStep            ; readXStep     = true;
				} else if(0 == token.compare("YStep"    )) {              iss >> yStep            ; readYStep     = true;
				} else if(0 == token.compare("ZStep"    )) {              iss >> zStep            ; readZStep     = true;
				} else if(0 == token.compare("AcqE1"    )) {              iss >> acqE1            ; readacqE1     = true;
				} else if(0 == token.compare("AcqE2"    )) {              iss >> acqE1            ; readacqE2     = true;
				} else if(0 == token.compare("AcqE3"    )) {              iss >> acqE2            ; readacqE3     = true;
				} else if(0 == token.compare("Euler"    )) {
					euStr = std::string(line);
					readEulerLine = true;
					while(iss >> token) {
						if       (0 == token.compare("Mag"      )) { iss >> mag         ; readMag       = true;
						} else if(0 == token.compare("Coverage" )) { iss >> cover      ; readCoverage  = true;
						} else if(0 == token.compare("Device"   )) { iss >> dev         ; readDevice    = true;
						} else if(0 == token.compare("KV"       )) { iss >> kv          ; readKv        = true;
						} else if(0 == token.compare("TiltAngle")) { iss >>     angle   ; readTiltAngle = true;
						} else if(0 == token.compare("TiltAxis" )) { iss >>     axis    ; readTiltAxis  = true;
						}
					}
				} else if(0 == token.compare("Phases"   )) {
					size_t numPhases;
					iss >> numPhases;
					for(size_t i = 0; i < numPhases; i++) {
						//phase lines are formatted as 'a;b;c\talpha;beta;gamma\tname\tlaue_group\tspacegroup'
						std::getline(is, line);//extract entire line from file
						std::replace(line.begin(), line.end(), ';', '\t');//convert to fully tab delimited to make things easier
						iss.str(line);//change istream buffer
						iss.clear();//clear any flags from istream

						//sanity check (phase lines should start with a digit)
						iss >> std::skipws;
						if(!std::isdigit(iss.peek())) throw std::runtime_error("end of ctf phases reached before epxected number of phases was read");

						//parse phase
						Phase p;
						iss >> p.lat[0] >> p.lat[1] >> p.lat[2];//a, b, c
						iss >> p.lat[3] >> p.lat[4] >> p.lat[5];//alpha, beta, gamma
						iss >> p.name   >> p.laue   >> p.space;///name, laue number, space group number
						phaseList.push_back(std::move(p));
					}
				} else throw std::runtime_error("unknown ctf header keyword '" + token + "'");
			}
		}

		//make sure all required parameters were read
		if(!readXCells || !readYCells) throw std::runtime_error("missing ctf dimensions");
		if(!readXStep  || !readYStep ) throw std::runtime_error("missing ctf resolution");
		if(!readacqE1 || !readacqE2 || !readacqE3 ) throw std::runtime_error("missing ctf header AcqE values");
		if(!readEulerLine) throw std::runtime_error("missing ctf euler angle convention line");
		if(!readMag      ) throw std::runtime_error("missing ctf magnification"              );
		if(!readCoverage ) throw std::runtime_error("missing ctf coverage"                   );
		if(!readDevice   ) throw std::runtime_error("missing ctf device"                     );
		if(!readKv       ) throw std::runtime_error("missing ctf accelerating voltage"       );
		if(!readTiltAngle) throw std::runtime_error("missing ctf tilt angle"                 );
		if(!readTiltAxis ) throw std::runtime_error("missing ctf tilt axis"                  );

		//handle missing optional values and return column headers
		if(!readZCells) zCells = 1;
		return cols;
	}

	//@brief     : read ctf data using an input stream
	//@param is  : input stream set data start
	//@param cols: column headers
	//@return    : number of points (rows) parsed
	size_t OrientationMap::readCtfData(std::istream& is, const std::vector<uint16_t>& cols) {
		//check if the columns are in the normal layout
		const uint16_t defaultLayout[11] = {
			CTF_Phase ,
			CTF_X     ,
			CTF_Y     ,
			CTF_Bands ,
			CTF_Error ,
			CTF_Euler1,
			CTF_Euler2,
			CTF_Euler3,
			CTF_MAD   ,
			CTF_BC    ,
			CTF_BS    ,
		};
		const bool isDefault = std::equal(cols.cbegin(), cols.cend(), defaultLayout);

		char line[512];
		const size_t totalPoints = xCells * yCells * zCells;

		if(isDefault) {
			for(size_t i = 0; i < totalPoints; i++) {
				char* data = line;
				if(!is.getline(line, sizeof(line))) return i;//get next line
				phase[i    ] = (uint8_t) std::strtoul(data, &data, 10);
				x    [i    ] =           std::strtof (data, &data    );
				y    [i    ] =           std::strtof (data, &data    );
				bands[i    ] = (uint8_t) std::strtoul(data, &data, 10);
				err  [i    ] =           std::strtof (data, &data    );
				eu   [3*i  ] =           std::strtof (data, &data    );
				eu   [3*i+1] =           std::strtof (data, &data    );
				eu   [3*i+2] =           std::strtof (data, &data    );
				mad  [i    ] =           std::strtof (data, &data    );
				bc   [i    ] =           std::strtof (data, &data    );
				bs   [i    ] =           std::strtof (data, &data    );
			}
		} else {
			//get a void pointer to each item in order and compute strides
			std::vector<size_t> strides ;//size of data type in each column
			std::vector<char* > pointers;//pointers to start of data in each column
			std::vector<size_t> isInt   ;//1/0 if each column is int/float
			for(const uint16_t& c : cols) {
				switch(c) {
					case CTF_Phase : strides.push_back(sizeof(phase.front())); pointers.push_back((char*)phase .data()                         ); isInt.push_back(1); break;
					case CTF_X     : strides.push_back(sizeof(x    .front())); pointers.push_back((char*)x     .data()                         ); isInt.push_back(0); break;
					case CTF_Y     : strides.push_back(sizeof(y    .front())); pointers.push_back((char*)y     .data()                         ); isInt.push_back(0); break;
					case CTF_Bands : strides.push_back(sizeof(bands.front())); pointers.push_back((char*)bands .data()                         ); isInt.push_back(1); break;
					case CTF_Error : strides.push_back(sizeof(err  .front())); pointers.push_back((char*)err   .data()                         ); isInt.push_back(0); break;
					case CTF_Euler1: strides.push_back(sizeof(eu   .front())); pointers.push_back((char*)eu    .data()                         ); isInt.push_back(0); break;
					case CTF_Euler2: strides.push_back(sizeof(eu   .front())); pointers.push_back((char*)eu    .data() + sizeof(eu.front())    ); isInt.push_back(0); break;
					case CTF_Euler3: strides.push_back(sizeof(eu   .front())); pointers.push_back((char*)eu    .data() + sizeof(eu.front()) * 2); isInt.push_back(0); break;
					case CTF_MAD   : strides.push_back(sizeof(mad  .front())); pointers.push_back((char*)mad   .data()                         ); isInt.push_back(0); break;
					case CTF_BC    : strides.push_back(sizeof(bc   .front())); pointers.push_back((char*)bc    .data()                         ); isInt.push_back(0); break;
					case CTF_BS    : strides.push_back(sizeof(bs   .front())); pointers.push_back((char*)bs    .data()                         ); isInt.push_back(0); break;
					default: throw std::runtime_error("unknown column type");
				}
			}

			//parse data
			for(size_t i = 0; i < totalPoints; i++) {
				char* data = line;
				if(!is.getline(line, sizeof(line))) return i;//get next line
				for(size_t j = 0; j < strides.size(); j++) {//loop over columns
					char* pData = pointers[j];//get pointer to column data
					if(isInt[j])//parse
						reinterpret_cast<uint8_t*>(pData)[i] = (uint8_t)std::strtoul(data, &data, 10);
					else
						reinterpret_cast<float  *>(pData)[i] =          std::strtof (data, &data    );
					pointers[j] += strides[j];//increment pointer
				}
			}
		}
		return totalPoints;
	}

	//@brief         : write data to a '.ctf' file
	//@param fileName: name of ctf file to write
	void OrientationMap::writeCtf(std::string fileName) {
		std::ofstream os(fileName.c_str());//open file
		writeCtfHeader(os);//write data
		writeCtfData(os);//write data
	}

	//@brief   : write ctf header
	//@param os: output stream to write the header to
	void OrientationMap::writeCtfHeader(std::ostream& os) {
		os << "Channel Text File\n";
		os << "Prj "      << project << '\n';
		os << "Author\t"  << author  << '\n';
		os << "JobMode\t" << jobMode << '\n';
		os << "XCells\t"  << xCells  << '\n';
		os << "YCells\t"  << yCells  << '\n';
		if(zCells != 1) os << "ZCells\t"  << zCells  << '\n';
		os << "XStep\t"   << xStep  << '\n';
		os << "YStep\t"   << yStep  << '\n';
		if(zCells != 1) os << "ZStep\t"   << zStep  << '\n';
		os << "AcqE1\t"   << acqE1  << '\n';
		os << "AcqE2\t"   << acqE2  << '\n';
		os << "AcqE3\t"   << acqE3  << '\n';
		os << "Euler angles refer to Sample Coordinate system (CS0)! ";
		os << "Mag\t"       << mag   << ' ';
		os << "Coverage\t"  << cover << ' ';
		os << "Device\t"    << dev   << ' ';
		os << "KV\t"        << kv    << ' ';
		os << "TiltAngle\t" << angle << ' ';
		os << "TiltAxis\t"  << axis  << '\n';
		os << "Phases\t"    << phaseList.size() << '\n';
		for(const Phase& p : phaseList) {
			os << p.lat[0] << ';'  << p.lat[1] << ';' << p.lat[2] << '\t';
			os << p.lat[3] << ';'  << p.lat[4] << ';' << p.lat[5] << '\t';
			os << p.name   << '\t' << p.laue  << '\t' << p.space  << '\n';
		}
		                   os << "Phase\t" ;
		if(!x    .empty()) os << "X\t"     ;
		if(!y    .empty()) os << "Y\t"     ;
		if(!bands.empty()) os << "Bands\t" ;
		if(!err  .empty()) os << "Error\t" ;
		                   os << "Euler1\t";
		                   os << "Euler2\t";
		                   os << "Euler3\t";
		if(!mad  .empty()) os << "MAD\t"   ;
		if(!bc   .empty()) os << "BC\t"    ;
		if(!bs   .empty()) os << "BS"      ;
		os << '\n';
	}

	//@brief       : write ctf data
	//@param os    : output stream to write data to
	void OrientationMap::writeCtfData(std::ostream& os) {
		os << std::fixed;
		for(size_t i = 0; i < phase.size(); i++) {
			                   os << (int) phase [i] << '\t';
			                   os << std::setprecision(3);
			if(!x    .empty()) os <<       x     [i] << '\t';
			if(!y    .empty()) os <<       y     [i] << '\t';
			if(!bands.empty()) os << (int) bands [i] << '\t';
			                   os << std::setprecision(6);
			if(!err  .empty()) os <<       err   [i] << '\t';
			                   os << std::setprecision(3);
			                   os <<       eu[3*i  ] << '\t';
			                   os <<       eu[3*i+1] << '\t';
			                   os <<       eu[3*i+2] << '\t';
			                   os << std::setprecision(0);
			if(!mad  .empty()) os <<       mad   [i] << '\t';
			if(!bc   .empty()) os <<       bc    [i] << '\t';
			if(!bs   .empty()) os <<       bs    [i];
			os << '\n';
		}
	}
}



/*
		//convert from degrees to radians and 1->0 phase indexing
		const T factor = M_PI / 180.0;
		std::for_each(eu.begin(), eu.end(), [factor](T& e){e *= factor;});
		std::for_each(phase.begin(), phase.end(), [](size_t& i){--i;});

		//correct for euler convention

		T* pEulers = eu.data();
		const T pi = T(2) * std::acos(T(0));
		const T twoPi = T(2) * std::acos(T(0));
		for(size_t i = 0; i < eu.size() / 3; i++)
			pEulers[3*i] = std::fmod(pEulers[3*i]+pi, twoPi);
		*/


#endif//_HKL_H_
