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

#ifndef _EBSD_PATTERN_H_
#define _EBSD_PATTERN_H_

#include <memory>
#include <vector>
#include <string>
#include <istream>
#include <mutex>

#include "idx/base.hpp"

namespace emsphinx {

	namespace ebsd {

		//@brief: abstract base class to hold patterns for indexing
		class PatternFile : public ImageSource {
			public:
		
				//@brief : get patterns held in file
				//@return: number of patterns
				size_t numPat() const {return num;}

				//@brief : check if the patterns need to be vertically flipped
				//@return: true/false if the patterns do/don't need to be flipped
				bool flipY() const {return flp;}

				//@brief   : set pattern shape
				//@param px: width of pattern in pixels
				//@param py: height of pattern in pixels
				//@param bt: bit depth of pixels
				void setShape(const size_t px, const size_t py, const Bits bt);

				//@brief   : set number of patterns
				//@param np: width of pattern in pixels
				void setNum(const size_t np) {num = np;}

				//@brief    : extract the next batch of patterns into a buffer
				//@param out: location to write
				//@param cnt: maximum number of patterns to extract
				//@return   : vector of the indices of each pattern extracted (e.g. {0,2,1,3} for the first 4 patterns but out of order)
				//@note     : implementations should be thread safe (only for other extract calls)
				virtual std::vector<size_t> extract(char * const out, const size_t cnt) const = 0;

				//@brief     : get pattern info from a file without reading entire thing
				//@param name: file name to get info from
				//@param w   : location to write width  (-1 if the file have width  information)
				//@param h   : location to write height (-1 if the file have height information)
				//@param bit : location to write bitdepth
				//@param num : location to write number (or file size in bytes if width/height are unknown)
				//@param aux : pattern file path for h5 files
				//@return    : true if the file can be read (based on extension), false otherwise
				static bool GetFileDims(const std::string name, int& w, int& h, Bits& bit, uint64_t& num, std::string aux = "");

				//@brief     : read experimental patterns from a file
				//@param name: name of pattern file to read (appropriate reader will be selected from the file name)
				//@param aux : auxiliary information (for hdf5 datasets this is the path to the dataset)
				//@param px  : width of patterns in pixels (0 to determine automatically)
				//@param py  : height of patterns in pixels (0 to determine automatically)
				//@return    : shared pointer to a pattern file
				static std::shared_ptr<PatternFile> Read(const std::string name, const std::string aux = "", const size_t px = 0, const size_t py = 0);

				//@brief     : search an h5 file for pattern datasets
				//@param name: name of pattern file to check (must be hdf5 type)
				//@return    : paths to suitable datasets (empty if none found), these are suitable for the 'aux' argument of the Read function
				static std::vector<std::string> SearchH5(const std::string name);

				//@brief     : get the vendor string from an hdf file
				//@param name: file to get vendor string from
				//@return    : vendor string or "" if none were found
				static std::string GetVendor(std::string name);

				//@brief    : read experimental patterns from individual image files
				//@param fmt: expression of pattern files to read (used with printf formatting to create files)
				//@param px : width of scan (how many patterns)
				//@param py : height of scan (how many patterns)
				//@return   : shared pointer to a pattern file
				static std::shared_ptr<PatternFile> FromImages(const std::string fmt, const size_t px, const size_t py);

			protected:
				size_t             num ;//number of patterns
				size_t             byt ;//bytes per pattern (w * h * bytes/pix)
				bool               flp ;//do the patterns need to be vertically flipped
				mutable std::mutex mut ;//mutex for thread safe access to extraction
		};

		//@brief: abstract intermediate class for cases where all patterns are stored contiguously
		class ContigousPatternFile : public PatternFile {
			protected:
				mutable size_t       idx ;//index of next pattern to extract
				mutable char const * ptr ;//pointer to next pattern to extract

			public:
				//@brief    : extract the next batch of patterns into a buffer
				//@param ptr: location to write
				//@param cnt: maximum number of patterns to extract
				//@return   : vector of the indices of each pattern extracted (e.g. {0,2,1,3} for the first 4 patterns but out of order)
				//@note     : implementations should be thread safe (only for other extract calls)
				std::vector<size_t> extract(char * const out, const size_t cnt) const;

				//@brief: virtual destructor (this is an abstract intermediate class)
				virtual ~ContigousPatternFile() = 0;
		};

		//@bierf: intermediate class for cases where all patterns are NOT stored contiguously (e.g. individual image file per patterb)
		class StreamedPatternFile : public PatternFile {
			protected:
				//all streamed pattern files need an input source
				mutable size_t idx;//index of next pattern to extract
				std::istream&  is ;

			public:
				//@brief    : extract the next batch of patterns into a buffer
				//@param ptr: location to write
				//@param cnt: maximum number of patterns to extract
				//@return   : vector of the indices of each pattern extracted (e.g. {0,2,1,3} for the first 4 patterns but out of order)
				//@note     : implementations should be thread safe (only for other extract calls)
				std::vector<size_t> extract(char * const out, const size_t cnt) const;

				//@brief  : construct a streamed pattern file
				//@param s: input stream to pull patterns from
				StreamedPatternFile(std::istream& s) : is(s) {}

				//@brief: virtual destructor (this is an abstract intermediate class)
				virtual ~StreamedPatternFile() = 0;
		};

		//@bierf: abstract intermediate class for cases where all patterns are NOT stored contiguously (e.g. individual image file per patter or chunked hdf5)
		class ChunkedPatternFile : public PatternFile {
			public:
				//@brief: virtual destructor (this is an abstract intermediate class)
				virtual ~ChunkedPatternFile() = 0;
		};

	}//ebsd

}//emsphinx

////////////////////////////////////////////////////////////////////////
//                        Non-abstract Classes                        //
////////////////////////////////////////////////////////////////////////

#include <fstream>

#include "util/sysnames.hpp"//fileSize
#include "util/bmp.hpp"
#include "H5Cpp.h"

namespace emsphinx {

	namespace ebsd {

		//@brief: class to hold all patterns in memory
		class BufferedPatternFile : public ContigousPatternFile {
			// std::vector<char> buff;//pattern buffer

			public:
			std::vector<char> buff;//pattern buffer
				//@brief: allocate buffer using current pattern size + number
				//@note : buffer is unitialized
				void allocate() {buff = std::vector<char>(imBytes() * num); ptr = data();}

				//@param : get pointer to underlying buffer
				//@return: pointer to buffer
				char* data() {return buff.data();}
		};

		//@brief: class to read patterns from an ifstream (this actually performs as well as a memory map on the linux systems tested)
		class IfStreamedPatternFile : public StreamedPatternFile {
			std::ifstream ifs ;//underlying file
			const size_t  fByt;//size of underlying file in bytes
			//if this class is modified to include pre/post padding bytes it will be significantly more flexible (and can accommodate the oxford format)

			public:
				//@brief     : open a memory mapped pattern file
				//@param name: name of file to map
				//@note      : all other members are uninitalized
				IfStreamedPatternFile(std::string name) : ifs(name, std::ios::in | std::ios::binary), fByt(fileSize(name)), StreamedPatternFile(ifs) {}

				//@brief   : construct the pointer to the data start from an offset in bytes
				//@param of: offset to data start in bytes
				//@note    : sets number of patterns using file size and pattern size
				void setOffset(const size_t of);
		};

		//@brief: oxford has enough peculiarities that I've written a separate class for now (it may be worth generalizing the streamed pattern file interface to make this unnecessary)
		class OxfordPatternFile : public PatternFile {
			mutable std::ifstream                       ifs;//underlying file
			        std::vector<size_t>                 idx;//order of patterns in file (usually close to 0,1,2,3,... but not always)
			mutable std::vector<size_t>::const_iterator nxt;//next pattern to extract

			public:
				//@brief     : open an oxford pattern file and get order
				//@param name: name of file to map
				OxfordPatternFile(std::string name);

				//@brief    : extract the next batch of patterns into a buffer
				//@param out: location to write
				//@param cnt: maximum number of patterns to extract
				//@return   : vector of the indices of each pattern extracted (e.g. {0,2,1,3} for the first 4 patterns but out of order)
				//@note     : implementations should be thread safe (only for other extract calls)
				std::vector<size_t> extract(char * const out, const size_t cnt) const {return extract(out, cnt, NULL, NULL);}

				//@brief    : extract the next batch of patterns into a buffer
				//@param out: location to write
				//@param cnt: maximum number of patterns to extract
				//@param vx : location to put x coordinates (or NULL to ignore coordinates)
				//@param vy : location to put y coordinates (or NULL to ignore coordinates)
				//@return   : vector of the indices of each pattern extracted (e.g. {0,2,1,3} for the first 4 patterns but out of order)
				//@note     : implementations should be thread safe (only for other extract calls)
				std::vector<size_t> extract(char * const out, const size_t cnt, std::vector<double>* vx, std::vector<double>* vy) const;
		};

	}//ebsd

}//emsphinx

////////////////////////////////////////////////////////////////////////
//                          Implementations                           //
////////////////////////////////////////////////////////////////////////

#include <algorithm>
#include <numeric>
#include <cctype>
#include <sstream>
#include <stack>

namespace emsphinx {
	
	namespace ebsd {

		namespace detail {
			//@brief     : get the extension of a file name
			//@param name: file name to get extension of
			//@return    : extension (all lower case)
			std::string getFileExt(const std::string name) {
				size_t pos = name.find_last_of(".");//find the last '.' in the name
				if(std::string::npos == pos) return "";//handle files with no extension
				std::string ext = name.substr(pos+1);//extract the file extension
				std::transform(ext.begin(), ext.end(), ext.begin(), [](char c){return std::tolower(c);});//convert to lowercase
				return ext;
			}
		}

		//@brief   : set pattern shape
		//@param px: width of pattern in pixels
		//@param py: height of pattern in pixels
		//@param bt: bit depth of pixels
		void PatternFile::setShape(const size_t px, const size_t py, const Bits bt) {
			bits = bt;//pixel type
			w    = px;//width of patterns
			h    = py;//height of patterns
			switch(bt) {
				case Bits::U8 : byt = 1; break;
				case Bits::U16: byt = 2; break;
				case Bits::F32: byt = 4; break;
				case Bits::UNK: 
				default       : throw std::runtime_error("unvalid bit type");
			}
			byt *= w * h;
		}

		//@brief: virtual constructor (this is an abstract intermediate class)
		ContigousPatternFile::~ContigousPatternFile() {}

		//@brief: virtual constructor (this is an abstract intermediate class)
		StreamedPatternFile::~StreamedPatternFile() {}

		//@brief: virtual constructor (this is an abstract intermediate class)
		ChunkedPatternFile::~ChunkedPatternFile() {}

		//@brief     : get pattern info from a file without reading entire thing
		//@param name: file name to get info from
		//@param w   : location to write width  (-1 if the file have width  information)
		//@param h   : location to write height (-1 if the file have height information)
		//@param bit : location to write bitdepth
		//@param num : location to write number (or file size in bytes if width/height are unknown)
		//@param aux : pattern file path for h5 files
		//@return    : true if the file can be read (based on extension), false otherwise
		bool PatternFile::GetFileDims(const std::string name, int& w, int& h, Bits& bit, uint64_t& num, std::string aux) {
			w = h = -1;
			num = 0;
			bit = Bits::UNK;
			const std::string ext = detail::getFileExt(name);
			const uint64_t fileBytes = fileSize(name);

			//handle easy types first (these have a header with all required data)
			if("up1" == ext || "up2" == ext || "ebsp" == ext) {
				std::shared_ptr<PatternFile> pat = Read(name);
				w   = (int)pat->width    ();
				h   = (int)pat->height   ();
				num =      pat->numPat   ();
				bit =      pat->pixelType();
			} else if("data" == ext) {
				bit = Bits::F32;
				num = fileBytes;
			} else if("h5" == ext || "hdf" == ext || "hdf5" == ext) {
				//open the dataset and get information
				hsize_t dims[3];
				H5::DataSet dSet = H5::H5File(name, H5F_ACC_RDONLY).openDataSet(aux);//open the file and get the dataset we're after
				if(3 != dSet.getSpace().getSimpleExtentNdims()) throw std::runtime_error("hdf pattern dataset must be 3D");
				dSet.getSpace().getSimpleExtentDims(dims);//read extent in each dimension
				w   = (int)dims[2];
				h   = (int)dims[1];
				num =      dims[0];

				//determine data type
				H5::DataType type = dSet.getDataType();
				if     (type == H5::DataType(H5::PredType::NATIVE_UINT8 )) bit = ImageSource::Bits::U8 ;
				else if(type == H5::DataType(H5::PredType::NATIVE_UCHAR )) bit = ImageSource::Bits::U8 ;
				else if(type == H5::DataType(H5::PredType::NATIVE_UINT16)) bit = ImageSource::Bits::U16;
				else if(type == H5::DataType(H5::PredType::NATIVE_FLOAT )) bit = ImageSource::Bits::F32;
				if(ImageSource::Bits::UNK == bit) throw std::runtime_error("only uint8, uint16, and float hdf patterns are supported");
			} else {
				return false;
			}
			return true;
		}

		//@brief     : read experimental patterns from a file
		//@param name: name of pattern file to read (appropriate reader will be selected from the file name)
		//@param aux : auxiliary information (for hdf5 datasets this is the path to the dataset)
		//@param px  : width of patterns in pixels (0 to determine automatically)
		//@param py  : height of patterns in pixels (0 to determine automatically)
		//@return    : shared pointer to a pattern file
		std::shared_ptr<PatternFile> PatternFile::Read(const std::string name, const std::string aux, const size_t px, const size_t py) {
			//extract extension
			const std::string ext = detail::getFileExt(name);
			const uint64_t fileBytes = fileSize(name);

			//read based on file extension
			if(0 == ext.compare("up1") || 0 == ext.compare("up2")) {
				//open file and read header
				int32_t header[4];
				{
					//read header
					std::ifstream is(name, std::ios::in | std::ios::binary);
					is.read((char*)header, 4 * sizeof(int32_t));
				}

				//parse header (logic here courtesy of stuart)
				int32_t& vers    = header[0];//up* file version
				int32_t& width   = header[1];//pattern width in pixels
				int32_t& height  = header[2];//pattern width in pixels
				int32_t& dStart  = header[3];//data start position
				if(width < 1 || height < 1 || width > 5000 || height > 5000 || dStart < 0) {
					if(vers > 2) {
						height = width;
						width = vers;
						dStart = 8;
					}
					if(width < 1 || height < 1 || width > 5000 || height > 5000 || dStart < 0) throw std::runtime_error("invalid UP file");
				}
				if((0 != px && px != width) || (0 != py && py != height)) throw std::runtime_error("patterns aren't expected shape");
				const Bits b = 0 == ext.compare("up1") ? ImageSource::Bits::U8 : ImageSource::Bits::U16;

				//construct with either ifstream
				std::shared_ptr<IfStreamedPatternFile> ptr = std::make_shared<IfStreamedPatternFile>(name);//memory map entire file
				ptr->setShape((size_t)width, (size_t)height, b);
				ptr->setOffset(dStart);
				ptr->flp = true;//EDAX files need to be flipped
				return std::shared_ptr<PatternFile>(ptr);
			} else if(0 == ext.compare("h5") || 0 == ext.compare("hdf") || 0 == ext.compare("hdf5")) {
				//open the dataset and get information
				H5::H5File file = H5::H5File(name, H5F_ACC_RDONLY);//open the file
				H5::DataSet dSet = file.openDataSet(aux);//get the dataset we're after
				hsize_t dims[3];
				if(3 != dSet.getSpace().getSimpleExtentNdims()) throw std::runtime_error("hdf pattern dataset must be 3D");
				dSet.getSpace().getSimpleExtentDims(dims);//read extent in each dimension
				if((0 != px && px != dims[2]) || (0 != py && py != dims[1])) throw std::runtime_error("patterns aren't expected shape");

				//determine if patterns need to be flipped using vendor
				bool vendorFlip = false;
				std::string vendor = GetVendor(name);
				if     ("EDAX"     == vendor) vendorFlip = true ;
				else if("Oxford"   == vendor) vendorFlip = false;
				else if("Bruker"   == vendor) vendorFlip = false;
				else if("DREAM.3D" == vendor) vendorFlip = false;
				else if("EMsoft"   == vendor) vendorFlip = true;
				else throw std::runtime_error("unknown EBSD vendor: " + vendor);

				//determine data type
				Bits bits = ImageSource::Bits::UNK;
				H5::DataType type = dSet.getDataType();
				if     (type == H5::DataType(H5::PredType::NATIVE_UINT8 )) bits = ImageSource::Bits::U8 ;
				else if(type == H5::DataType(H5::PredType::NATIVE_UCHAR )) bits = ImageSource::Bits::U8 ;
				else if(type == H5::DataType(H5::PredType::NATIVE_UINT16)) bits = ImageSource::Bits::U16;
				else if(type == H5::DataType(H5::PredType::NATIVE_FLOAT )) bits = ImageSource::Bits::F32;
				if(ImageSource::Bits::UNK == bits) throw std::runtime_error("only uint8, uint16, and float hdf patterns are supported");

				//get offset from file start
				haddr_t offset;
				try {
					offset = dSet.getOffset();//get the offset to the dataset start
				} catch(H5::DataSetIException&) {
					offset = (haddr_t)-1;//HAADR_UNDEF if not contigous (<0)
				}

				//try to raw access the dataset
				H5::DSetCreatPropList props = dSet.getCreatePlist();
				if(0 == props.getExternalCount()                                               &&   //we can't use a memory map if the dataset references external files
				   ((H5D_COMPACT == props.getLayout() || H5D_CONTIGUOUS == props.getLayout())) &&   //we can't memory map datasets unless they are contigous in memory
				   0 != props.getNfilters()                                                    &&   //we can't memory map compressed datasets
				   offset != (haddr_t)-1                                                         ) {//negative value (HADDR_UNDEF) means that hdf5 couldn't compute the data offset
					//we can try memory mapping the file
					try {
						std::shared_ptr<IfStreamedPatternFile> ptr = std::make_shared<IfStreamedPatternFile>(name);//memory map entire file
						ptr->setShape(dims[2], dims[1], bits);
						ptr->setOffset(dSet.getOffset());
						ptr->setNum(dims[0]);
						ptr->flp = vendorFlip;
						return std::shared_ptr<PatternFile>(ptr);
					} catch (...) {
						//fall back to normal hdf5 reading if we failed
					}
				}

				//if we made it this far we either can't get raw access or failed to do so
				//just read entire dataset into memory (it may be good to add a ChunkedPatternFile derived option for large files)
				std::shared_ptr<BufferedPatternFile> ptr = std::make_shared<BufferedPatternFile>();
				ptr->setShape(dims[2], dims[1], bits);
				ptr->setNum(dims[0]);
				ptr->allocate();
				dSet.read((void*)ptr->data(), H5::PredType::NATIVE_UINT8);
				ptr->flp = vendorFlip;
				return std::shared_ptr<PatternFile>(ptr);
			} else if(0 == ext.compare("data")) {
				//raw 32bit float patterns with no header, construct with ifstream
				std::shared_ptr<IfStreamedPatternFile> ptr = std::make_shared<IfStreamedPatternFile>(name);//memory map entire file
				ptr->setShape(px, py, Bits::F32);
				ptr->setOffset(0);
				ptr->flp = false;
				return std::shared_ptr<PatternFile>(ptr);
			} else if(0 == ext.compare("ebsp")) {//oxford raw pattern format
				//the format for this reader was parsed from the dream3d reader in OxfordReader.cpp
				// const std::string ext2 = detail::getFileExt(aux);
				// if(0 != ext2.compare("dat")) throw std::runtime_error(name + " doesn't have associated `dat' file");
				//there is a binary ctf type file in the associated .dat file
				//for each pixel in the orientation map there are 25 bytes in a planar format:
				//  -uint8_t : phase
				//  -3x float: euler angles
				//  -float   : MAD
				//  -uint8_t : BC
				//  -uint8_t : BS
				//  -uint8_t : bands
				//  -uint8_t : unknown
				//  -float   : unknown
				//with first numPat bytes with the phase of each pixel, then 12*numPat bytes with eulers, etc
				std::shared_ptr<OxfordPatternFile> ptr = std::make_shared<OxfordPatternFile>(name);
				return ptr;
			}
			throw std::runtime_error("couldn't find reader for '" + name + "'");
		}

		std::vector<std::string> PatternFile::SearchH5(const std::string name) {
			H5::Exception::dontPrint();//silence caught errors

			//@brief: function to loop over an hdf5 location adding children names to a list
			//@param gid   : location id
			//@param name  : name of child
			//@param opData: operator data (pointer to vector to store names in)
			H5G_iterate_t iterFunc = [](hid_t gid, const char * name, void *opData)->herr_t{
				std::vector<std::string>* pDsets = (std::vector<std::string>*)opData;
				pDsets->push_back(name);
				return 0;
			};

			//open file and allocate space for valid names
			std::vector<std::string> choices;
			H5::H5File file(name.c_str(), H5F_ACC_RDONLY);

			//recursively loop over locations in files checking for validity
			std::stack<std::string> locs;
			locs.push(".");//add root
			while(!locs.empty()) {
				//get top of stack
				std::string loc = locs.top();
				locs.pop();

				//get all children of top
				int idx = 0;
				std::vector<std::string> children;
				file.iterateElems(loc.c_str(), &idx, iterFunc, &children);

				//loop over children classifying
				for(std::string& child : children) {
					std::string name = loc + "/" + child;//build full path
					try {
						file.openGroup(name.c_str());//try to open as a group
						locs.push(name);//if openGroup didn't throw it is a group, add to recursion
					} catch(...) {//this location is a dataset, check if it is valid
						H5::DataSpace space = file.openDataSet(name.c_str()).getSpace();
						if(3 == space.getSimpleExtentNdims()) {//pattern datasets must be 3d
							hsize_t dims[3];
							space.getSimpleExtentDims(dims);
							if(dims[2] > 4) choices.push_back(name);//ebsd H5 files can have RGB(A) images for e.g. coordinate system
						}
					}
				}
			}
			return choices;
		}

		//@brief     : get the vendor string from an hdf file
		//@param name: file to get vendor string from
		//@return    : vendor string or "" if none were found
		std::string PatternFile::GetVendor(std::string name) {
			//find the manufacturer dataset
			H5::H5File file = H5::H5File(name, H5F_ACC_RDONLY);//open the file
			int idx = 0;
			int pad = -1;//1 for "Manufacturer", 2 for " Manufacturer", 0 for nothing
			file.iterateElems(".", &idx, [](int, const char* nm, void* pPad)->int{
				if(std::string("Manufacturer") == nm) {
					*(int*)pPad = 0;//save type
					return 1;//non zero -> stop searching
				} else if(std::string(" Manufacturer") == nm) {//some EDAX files have a stray leading space...
					*(int*)pPad = 1;//save type
					return 1;//non zero -> stop searching
				}
				return 0;//0 -> keep searching
			}, &pad);
			if(-1 == pad) throw std::runtime_error(name + " doesn't have a Manufacturer string");
			std::string manStr = pad == 1 ? " Manufacturer" : "Manufacturer";

			//determine the vendor and use to set pattern flip flag
			try {
				std::string vendor;
				file.openDataSet(manStr).read( vendor, H5::StrType(0, H5T_VARIABLE), H5::DataSpace(H5S_SCALAR));
				return vendor;
			} catch (H5::Exception&) {//H5T_VARIABLE + std::string seems to fail here...
				std::string vendor(128, 0);
				file.openDataSet(manStr).read((void*)vendor.data(), H5::StrType(0, vendor.size()), H5::DataSpace(H5S_SCALAR));//try reading as fixed length
				vendor.erase(std::find(vendor.begin(), vendor.end(), '\0'), vendor.end());//get rid of extra padding nulls
				return vendor;
			}
		}

		//@brief    : read experimental patterns from individual image files
		//@param fmt: expression of pattern files to read (used with printf formatting to create files)
		//@param px : width of scan (how many patterns)
		//@param py : height of scan (how many patterns)
		//@return   : shared pointer to a pattern file
		std::shared_ptr<PatternFile> PatternFile::FromImages(const std::string fmt, const size_t px, const size_t py) {
			//generate a list of filenames
			std::vector<std::string> fileNames;
			std::string name;//location to write filename
			for(size_t j = 0; j < py; j++) {
				for(size_t i = 0; i < px; i++) {
					const size_t bytes = std::snprintf(NULL, 0, fmt.data(), i, j) + 1;//how much space will it take to create the string (+ null character)
					if(bytes > name.size()) name.resize(bytes);
					std::snprintf(&name[0], name.size(), fmt.data(), i, j);
					fileNames.push_back(name);
				}
			}
			if(fileNames.empty()) throw std::runtime_error("no file names created for image files");

			//read based on file extension
			const std::string ext = detail::getFileExt(fmt);//extract extension
			if(0 == ext.compare("bmp")) {
				//read first pattern header
				bmp::Header hdr;
				{
					std::ifstream is(fileNames.front(), std::ios::in | std::ios::binary);
					hdr.read(is);
				}

				//parse shape etc
				std::shared_ptr<BufferedPatternFile> ptr = std::make_shared<BufferedPatternFile>();
				const size_t bitCount = hdr.bitCount;
				switch(bitCount) {
					case 8 : break;
					case 24: break;//potentially a grayscale image as rgb
					default: throw std::runtime_error("unsupported bitmap pattern bit depth");
				}
				ptr->setShape(hdr.width, hdr.height, ImageSource::Bits::U8);
				ptr->setNum(fileNames.size());

				//allocate data and get pointer
				ptr->allocate();
				const size_t imBytes = ptr->imBytes();
				char* buff = ptr->data();

				//loop over images reading
				for(const std::string nm : fileNames) {
					//open file and read header
					std::ifstream is(nm, std::ios::in | std::ios::binary);
					hdr.read(is);

					//make sure the pattern size is the same
					if(hdr.width != ptr->width() || hdr.height != ptr->height()) throw std::runtime_error("pattern shape mismatch");
					if(hdr.bitCount != bitCount) throw std::runtime_error("pattern bitdepth mismatch");

					//read the pattern
					hdr.readImage(is, buff, true);
					buff += imBytes;
				}

				std::ofstream os("test.raw", std::ios::out | std::ios::binary);
				os.write(ptr->data(), ptr->imBytes() * fileNames.size());

				//return upcasted pointer
				return std::shared_ptr<PatternFile>(ptr);
			}

			throw std::runtime_error("couldn't find reader for '" + name + "' (currently only bmp images are supported)");
		}

		//@brief    : extract the next batch of patterns into a buffer
		//@param ptr: location to write
		//@param cnt: maximum number of patterns to extract
		//@return   : index of first and last pattern extract (exlusive so second - first patterns were extracted)
		//@note     : implementations should be thread safe (only for other extract calls)
		std::vector<size_t> ContigousPatternFile::extract(char * const out, const size_t cnt) const {
			std::lock_guard<std::mutex> lock(mut);//only 1 thread can get patterns at once
			const size_t numExt = std::min(numPat() - idx, cnt);//determine how many patterns we can extract
			char const * const ptrNew = ptr + numExt * imBytes();//get end of range to copy
			std::copy(ptr, ptrNew, out);//copy patterns to output
			std::vector<size_t> res(numExt);//build vector to hold extraction index list
			std::iota(res.begin(), res.end(), idx);//fill index list
			idx += numExt;//update index of next patterns
			ptr = ptrNew;//update pointer to next pattern
			if(idx == numPat()) ptr = 0;//we've run out of patterns
			return res;
		}

		//@brief    : extract the next batch of patterns into a buffer
		//@param ptr: location to write
		//@param cnt: maximum number of patterns to extract
		//@return   : vector of the indices of each pattern extracted (e.g. {0,2,1,3} for the first 4 patterns but out of order)
		//@note     : implementations should be thread safe (only for other extract calls)
		std::vector<size_t> StreamedPatternFile::extract(char * const out, const size_t cnt) const {
			std::lock_guard<std::mutex> lock(mut);//only 1 thread can get patterns at once
			const size_t numExt = std::min(numPat() - idx, cnt);//determine how many patterns we can extract
			const size_t readBytes = numExt * imBytes();//get end of range to copy
			is.read(out, readBytes);//read from istream into the output
			std::vector<size_t> res(numExt);//build vector to hold extraction index list
			std::iota(res.begin(), res.end(), idx);//fill index list
			idx += numExt;//update index of next patterns
			return res;
		}

		//@brief   : construct the pointer to the data start from an offset in bytes
		//@param of: offset to data start in bytes
		//@note    : sets number of patterns using file size and pattern size
		void IfStreamedPatternFile::setOffset(const size_t of) {
			num = (fByt - of) / imBytes();
			ifs.seekg(of);
			idx = 0;
		}

		//@brief     : open an oxford pattern file and get order
		//@param name: name of file to map
		OxfordPatternFile::OxfordPatternFile(std::string name) : ifs(name, std::ios::in | std::ios::binary) {//open the input file
			//the oxford format starts an 8 bytes header followed by a sequence of uint64_t's with the offset to the data for each pattern
			//the first pattern is always immediately after the offsets
			//start by reading the header and first offset to determine the pattern count
			uint64_t offset;
			uint8_t header[8];
			ifs.read((char*)header, 8);//read header bytes
			if(//0xFF != header[0] || sometimes this is 0xFE...
			   0xFF != header[1] ||
			   0xFF != header[2] ||
			   0xFF != header[3] ||
			   0xFF != header[4] ||
			   0xFF != header[5] ||
			   0xFF != header[6] ||
			   0xFF != header[7]) throw std::runtime_error("unexpected header for " + name);//check header
			if(!ifs.read((char*)&offset, sizeof(offset))) throw std::runtime_error("failed to read first pattern position from " + name);//read offset to first pattern
			const uint64_t numPat = (offset - 8) / 8;//compute pattern count from offset to first pattern

			//now read all the offsets at once
			ifs.seekg(8);//go back to first index
			std::vector<uint64_t> offsets(numPat);//save space to hold all offsets
			if(!ifs.read((char*)offsets.data(), sizeof(uint64_t) * offsets.size())) throw std::runtime_error("failed to read offsets from " + name);

			//next seek to the first pattern and get some info
			//each pattern block has a 16 byte header, some (maybe 0) padding space, the pattern, and then 18 tail bytes:
			// uint8_t x position
			// double  x position (in microns)
			// uint8_t y position
			// double  y position (in microns)
			ifs.seekg(offsets[0]);
			uint32_t leadIn, width, height, bytes;
			if(!ifs.read((char*)&leadIn, sizeof(uint32_t))) throw std::runtime_error("failed to read lead in of first pattern in " + name);
			if(!ifs.read((char*)&height, sizeof(uint32_t))) throw std::runtime_error("failed to read width of first pattern in " + name);
			if(!ifs.read((char*)&width , sizeof(uint32_t))) throw std::runtime_error("failed to read height of first pattern in " + name);
			if(!ifs.read((char*)&bytes , sizeof(uint32_t))) throw std::runtime_error("failed to read bytes of first pattern in " + name);

			//parse the bit depth
			Bits bits;
			if(bytes == width * height) {
				bits = Bits::U8;
			} else if(bytes == width * height * 2) {
				bits = Bits::U16;
			} else {
				throw std::runtime_error("couldn't determine pixel type for patterns in " + name);
			}

			//next compute the size of a the first pattern
			uint64_t blockBytes = 16     +  //4 * uint32_t for leadIn, width, height, + bytes
			                      leadIn + //padding bytes between the header + data
			                      bytes  + //actual image data
			                      18     ; //tail data

			//now convert from offsets to index of each pattern in order
			idx.assign(numPat, numPat);//save space to hold all offsets and fill with numPat (largest index + 1)
			for(size_t i = 0; i < numPat; i++) {//loop over offsets
				uint64_t off = offsets[i] - offset;//convert from absolute offset to offset to from first pattern
				if(0 != off % blockBytes) throw std::runtime_error("inconsistent block sizes aren't currently supported for .ebsp files");
				off /= blockBytes;//we know of the offset of pattern i (in patterns)
				idx[off] = i;//what we actually need to know is which pattern is at each position
			}

			//make sure we got all the patterns
			for(const size_t& i : idx) {
				if(i == numPat) throw std::runtime_error("not all patterns were found in " + name);
			}

			//if we made it this far we can stream the file, save all the information required
			setShape(width, height, bits);
			setNum(numPat);
			flp = false;

			//move to first pattern
			nxt = idx.cbegin();
			ifs.seekg(offset);
		}

		//@brief    : extract the next batch of patterns into a buffer
		//@param out: location to write
		//@param cnt: maximum number of patterns to extract
		//@param vx : location to put x coordinates (or NULL to ignore coordinates)
		//@param vy : location to put y coordinates (or NULL to ignore coordinates)
		//@return   : vector of the indices of each pattern extracted (e.g. {0,2,1,3} for the first 4 patterns but out of order)
		//@note     : implementations should be thread safe (only for other extract calls)
		std::vector<size_t> OxfordPatternFile::extract(char * const out, const size_t cnt, std::vector<double>* vx, std::vector<double>* vy) const {
			if(nxt == idx.cend()) return std::vector<size_t>();//handle trivial case (nothing available)

			//extract as many patterns as possible (up to cnt)
			std::lock_guard<std::mutex> lock(mut);//only 1 thread can get patterns at once
			std::vector<size_t> ret;
			uint32_t leadIn, width, height, bytes;
			const size_t pByt = imBytes();
			for(size_t i = 0; i < cnt; i++) {//loop over requested values
				//read block header
				ifs.read((char*)&leadIn, sizeof(uint32_t));
				ifs.read((char*)&height, sizeof(uint32_t));
				ifs.read((char*)&width , sizeof(uint32_t));
				ifs.read((char*)&bytes , sizeof(uint32_t));

				//make sure the pattern is the same shape as the first one
				if(width != this->width() || height != this->height() || bytes != pByt) {
					std::stringstream ss;
					ss << "pattern " << std::distance(idx.cbegin(), nxt) << " isn't the same shape as the first pattern";
					throw std::runtime_error(ss.str());
				}

				//skip lead in and read
				ifs.ignore(leadIn);//skip padding space
				ifs.read(out + i * pByt, pByt);//read the pattern
				// ifs.ignore(18);//skip tail data

				double x, y;
				uint8_t ix, iy;
				ifs.read((char*)&ix, sizeof(uint8_t));
				ifs.read((char*)& x, sizeof(double ));
				ifs.read((char*)&iy, sizeof(uint8_t));
				ifs.read((char*)& y, sizeof(double ));
				if(NULL != vx) vx->push_back(x);
				if(NULL != vy) vy->push_back(y);

				//save index of pattern and advance
				ret.push_back(*nxt);
				++nxt;
			}
			return ret;
		}

	}//ebsd

}//emsphinx

#endif//_EBSD_PATTERN_H_
