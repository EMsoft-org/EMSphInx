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

#include <sstream>
#include <fstream>
#include <cmath>

#include "modality/ebsd/pattern.hpp"
#include "util/timer.hpp"
#include "H5Cpp.h"

//@brief    : in place bin a pattern (binned pattern will be converted to floats and summed)
//@param in : raw buffer holding pattern (of pixel type T)
//@param out: location to write output (cannot overlap with input)
//@param bin: binning factor
//@param w  : binned pattern width  (input width  is w * bin)
//@param h  : binned pattern height (input height is h * bin)
template <typename T>
void binFloat(T const* in, float* out, const size_t bin, const size_t w, const size_t h) {
	if(0 == bin){
		throw std::runtime_error("bin factor must be positive");
	} else if(1 == bin) {
		std::transform(in, in + w * h, out, [](const T& i){return (float)i;});//handle trivial case
	} else {
		for(size_t j = 0; j < h; j++) {//loop over output rows
			std::fill(out, out + w, 0.0f);//initialize output with zero
			for(size_t k = 0; k < bin; k++) {//loop over input rows associated with this output row
				for(size_t i = 0; i < w; i++) {//loop over output columns
					out[i] += std::accumulate(in, in + bin, 0.0f);
					in += bin;
				}
			}
			out += w;
		}
	}
}

//@brief    : in place bin a pattern (without changing its type, result with be averaged)
//@param in : raw buffer holding pattern (of pixel type T)
//@param out: location to write output (cannot overlap with input)
//@param bin: binning factor
//@param w  : binned pattern width  (input width  is w * bin)
//@param h  : binned pattern height (input height is h * bin)
template <typename T>
void binAvg(T const* in, T* out, const size_t bin, const size_t w, const size_t h) {
	if(0 == bin){
		throw std::runtime_error("bin factor must be positive");
	} else if(1 == bin) {
		std::copy(in, in + w * h, out);
		// std::transform(in, in + w * h, out, [](const T& i){return (float)i;});//handle trivial case
	} else {
		std::vector<double> row(w);//space large enough to hold a single row
		for(size_t j = 0; j < h; j++) {//loop over output rows
			std::fill(row.begin(), row.end(), 0.0);//initialize output with zero
			for(size_t k = 0; k < bin; k++) {//loop over input rows associated with this output row
				for(size_t i = 0; i < w; i++) {//loop over output columns
					row[i] += std::accumulate(in, in + bin, 0.0);
					in += bin;
				}
			}
			std::for_each(row.begin(), row.end(), [bin](double& v){v /= bin * bin;});//convert from sum to average
			if(std::is_integral<T>::value) std::for_each(row.begin(), row.end(), [](double& v){v = std::round(v);});//round to nearest int if needed
			for(size_t i = 0; i < w; i++) out[i] = (T)row[i];//convert back to T
			out += w;
		}
	}
}

//@brief    : in place flip a pattern
//@param buf: raw buffer holding pattern
//@param w  : pattern width
//@param h  : pattern height
//@param b  : bytes per pixel
void flipPat(char*buf, const size_t w, const size_t h, const size_t b) {
	const size_t rowBytes = w * b;
	char* revBuf = buf + w * (h-1) * b;//pointer to start of last row
	for(size_t j = 0; j < h/2; j++) {
		std::swap_ranges(buf, buf + rowBytes, revBuf);//swap rows
		buf    += rowBytes;//increment  lower row
		revBuf -= rowBytes;//decremenet upper row
	}
}

int main(int argc, char *argv[]) {
	const bool binToFloat = false;//should binned values be converted up to a float or kept as their current type
	const bool flip = true;//should binned values be converted up to a float or kept as their current type

	//sanity check argument count
	if(!(3 == argc || 4 == argc)) {
		std::cout << "usage: " << argv[0] << " inputFile outputFile [binning]\n";
		std::cout << "\tinputFile  - pattern file to read (*.up1, *.up2, *.data, or *.ebsp)\n";
		std::cout << "\toutputFile - output file (*.hdf)\n";
		std::cout << "\tbinning    - [optional] binning size (must evenly divide into pattern size)\n";
		return EXIT_FAILURE;
	}

	//prase arguments
	size_t binning = 1;
	std::string inputFile (argv[1]);
	std::string outputFile(argv[2]);
	if(argc == 4) {
		if(!(std::istringstream(argv[3]) >> binning)) {
			std::cout << "couldn't parse binning from 3rd argument `" << argv[3] << "'\n";
			return EXIT_FAILURE;
		}
	}

	//sanity check binning size
	if(0 == binning) {//dont % 0
		std::cout << "binning must be a positive\n";
		return EXIT_FAILURE;
	}

	//attempt to open pattern file (handles format via extension)
	Timer t;
	std::shared_ptr<emsphinx::ebsd::PatternFile> pats = emsphinx::ebsd::PatternFile::Read(inputFile);
	if(emsphinx::ImageSource::Bits::UNK == pats->pixelType()) {
		std::cout << "unknown pixel type\n";
		return EXIT_FAILURE;
	}

	//now print some info
	std::cout << "found " << pats->numPat() << " patterns in " << inputFile << ":\n";
	std::cout << "\twidth : " << pats->width () << '\n';
	std::cout << "\thegiht: " << pats->height() << '\n';
	std::cout << "\ttype  : ";
	switch(pats->pixelType()) {
		case emsphinx::ImageSource::Bits::U8 : std::cout << "8 bit\n"  ; break;
		case emsphinx::ImageSource::Bits::U16: std::cout << "16 bit\n" ; break;
		case emsphinx::ImageSource::Bits::F32: std::cout << "float\n"  ; break;
		case emsphinx::ImageSource::Bits::UNK: // intentional fall through
		default: throw std::logic_error("unknown pixel type");
	}
	std::cout << "\tbytes : " << pats->imBytes() << "\n";
	double mb = double(pats->imBytes() * pats->numPat()) / (1024*1024);
	std::cout << "total size: ";
	if(mb > 1024) {
		std::cout << mb / 1024 << " GB\n";
	} else {
		std::cout << mb        << " MB\n";
	}

	//make sure our binning is valid
	if(0 != (pats->width() % binning) || 0 != (pats->height() % binning)) {
		std::cout << "binning doesn't evenly divide pattern size\n";
		return EXIT_FAILURE;
	}

	std::cout << t.poll() << "s to parse pattern file\n";

	//build up data type and space
	H5::DataType dType;
	switch(pats->pixelType()) {//determine data set type
		case emsphinx::ImageSource::Bits::U8 : dType = H5::PredType::NATIVE_UINT8 ; break;
		case emsphinx::ImageSource::Bits::U16: dType = H5::PredType::NATIVE_UINT16; break;
		case emsphinx::ImageSource::Bits::F32: dType = H5::PredType::NATIVE_FLOAT ; break;
		case emsphinx::ImageSource::Bits::UNK: // intentional fall through
		default: throw std::logic_error("unknown pixel type");
	}
	if(binToFloat) {
		if(binning != 1) dType = H5::PredType::NATIVE_FLOAT;//if we're binning bin into floats
	}

	hsize_t dims[3] = {pats->numPat(), pats->height() / binning, pats->width() / binning};
	H5::DataSpace dSpace(3, dims);

	//open output file and create dataset for patterns
	H5::H5File hf(outputFile, H5F_ACC_TRUNC);
	H5::DSetCreatPropList props;
	props.setAllocTime(H5D_ALLOC_TIME_EARLY);//this lets us get a raw binary offset later with getOffset()
	H5::DataSet dSet = hf.createDataSet("patterns", dType, dSpace, props);

	//create data spaces and arrays for hyperslab selection
	dims[0] = 1;//buffer holds a single pattern
	hsize_t count [3] = {1, dims[1], dims[2]};//block count
	hsize_t start [3] = {      0,       0, 0};//block sizes
	H5::DataSpace buffSpace(3, dims);//create H5 data space for memory buffer

	//finally loop over patterns writing into file
	std::vector<char> buff(pats->imBytes());
	if(1 == binning) {//no binning ==> raw copy
		for(size_t i = 0; i < pats->numPat(); i++) {
			start[0] = pats->extract(buff.data(), 1)[0];//read the next pattern and get index (shouldn't be empty since we're only extracting with one thread)
			if(flip) flipPat(buff.data(), pats->width(), pats->height(), pats->pixBytes());
			dSpace.selectHyperslab(H5S_SELECT_SET, count, start);//tell hdf5 where to write
			dSet.write((void*) buff.data(), dType, buffSpace, dSpace);//copy the pattern out
		}
	} else {//binning
		std::vector<float> binBuf((size_t)dims[1] * (size_t)dims[2]);//this is big enough for either case

		if(binToFloat) {
			switch(pats->pixelType()) {
				case emsphinx::ImageSource::Bits::U8 :
					for(size_t i = 0; i < pats->numPat(); i++) {
						start[0] = pats->extract(buff.data(), 1)[0];//read the next pattern and get index (shouldn't be empty since we're only extracting with one thread)
						if(flip) flipPat(buff.data(), pats->width(), pats->height(), pats->pixBytes());
						binFloat((uint8_t*) buff.data(), binBuf.data(), binning, (size_t)dims[2], (size_t)dims[1]);
						dSpace.selectHyperslab(H5S_SELECT_SET, count, start);//tell hdf5 where to write
						dSet.write((void*) binBuf.data(), H5::PredType::NATIVE_FLOAT, buffSpace, dSpace);//copy the pattern out
					}
				break;

				case emsphinx::ImageSource::Bits::U16:
					for(size_t i = 0; i < pats->numPat(); i++) {
						start[0] = pats->extract(buff.data(), 1)[0];//read the next pattern and get index (shouldn't be empty since we're only extracting with one thread)
						if(flip) flipPat(buff.data(), pats->width(), pats->height(), pats->pixBytes());
						binFloat((uint16_t*) buff.data(), binBuf.data(), binning, (size_t)dims[2], (size_t)dims[1]);
						dSpace.selectHyperslab(H5S_SELECT_SET, count, start);//tell hdf5 where to write
						dSet.write((void*) binBuf.data(), H5::PredType::NATIVE_FLOAT, buffSpace, dSpace);//copy the pattern out
					}
				break;

				case emsphinx::ImageSource::Bits::F32:
					for(size_t i = 0; i < pats->numPat(); i++) {
						start[0] = pats->extract(buff.data(), 1)[0];//read the next pattern and get index (shouldn't be empty since we're only extracting with one thread)
						if(flip) flipPat(buff.data(), pats->width(), pats->height(), pats->pixBytes());
						binFloat((float*) buff.data(), binBuf.data(), binning, (size_t)dims[2], (size_t)dims[1]);
						dSpace.selectHyperslab(H5S_SELECT_SET, count, start);//tell hdf5 where to write
						dSet.write((void*) binBuf.data(), H5::PredType::NATIVE_FLOAT, buffSpace, dSpace);//copy the pattern out
					}
				break;

				case emsphinx::ImageSource::Bits::UNK: // intentional fall through
				default: throw std::logic_error("unknown pixel type");
			}
		} else {
			switch(pats->pixelType()) {
				case emsphinx::ImageSource::Bits::U8 :
					for(size_t i = 0; i < pats->numPat(); i++) {
						start[0] = pats->extract(buff.data(), 1)[0];//read the next pattern and get index (shouldn't be empty since we're only extracting with one thread)
						if(flip) flipPat(buff.data(), pats->width(), pats->height(), pats->pixBytes());
						binAvg((uint8_t *) buff.data(), (uint8_t *)binBuf.data(), binning, (size_t)dims[2], (size_t)dims[1]);
						dSpace.selectHyperslab(H5S_SELECT_SET, count, start);//tell hdf5 where to write
						dSet.write((void*) binBuf.data(), dType, buffSpace, dSpace);//copy the pattern out
					}
				break;

				case emsphinx::ImageSource::Bits::U16:
					for(size_t i = 0; i < pats->numPat(); i++) {
						start[0] = pats->extract(buff.data(), 1)[0];//read the next pattern and get index (shouldn't be empty since we're only extracting with one thread)
						if(flip) flipPat(buff.data(), pats->width(), pats->height(), pats->pixBytes());
						binAvg((uint16_t*) buff.data(), (uint16_t*)binBuf.data(), binning, (size_t)dims[2], (size_t)dims[1]);
						dSpace.selectHyperslab(H5S_SELECT_SET, count, start);//tell hdf5 where to write
						dSet.write((void*) binBuf.data(), dType, buffSpace, dSpace);//copy the pattern out
					}
				break;

				case emsphinx::ImageSource::Bits::F32:
					for(size_t i = 0; i < pats->numPat(); i++) {
						start[0] = pats->extract(buff.data(), 1)[0];//read the next pattern and get index (shouldn't be empty since we're only extracting with one thread)
						if(flip) flipPat(buff.data(), pats->width(), pats->height(), pats->pixBytes());
						binAvg((float   *) buff.data(), (float   *)binBuf.data(), binning, (size_t)dims[2], (size_t)dims[1]);
						dSpace.selectHyperslab(H5S_SELECT_SET, count, start);//tell hdf5 where to write
						dSet.write((void*) binBuf.data(), dType, buffSpace, dSpace);//copy the pattern out
					}
				break;

				case emsphinx::ImageSource::Bits::UNK: // intentional fall through
				default: throw std::logic_error("unknown pixel type");
			}
		}
	}
	std::cout << t.poll() << "s to repack patterns\n";
	return EXIT_SUCCESS;
}
