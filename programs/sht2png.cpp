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

#include <fstream>
#include <algorithm>
#include <cmath>

#include "idx/master.hpp"
#include "sht_file.hpp"

#define MINIZ_NO_STDIO
#define MINIZ_NO_TIME
#define MINIZ_NO_ZLIB_APIS
#include "miniz/miniz.c"

//@brief     : write a gray or rgb image to a png
//@param im  : image to write
//@param w   : image width
//@param h   : image height
//@param samp: samples per pixel (1 for gray, 3 for rgb)
//@param name: file name to write to (*.png)
void writePng(uint8_t* im, const size_t w, const size_t h, const size_t samp, std::string name) {
	//convert to png in memory
	size_t pngSize = 0;
	const mz_uint compr = MZ_BEST_COMPRESSION;//compression level [0,10]
	const mz_bool flip  = MZ_FALSE;//flip the image?
	void *pPNG_data = tdefl_write_image_to_png_file_in_memory_ex((void*)im, (int)w, (int)h, (int)samp, &pngSize, compr, flip);
	if(!pPNG_data) throw std::runtime_error("failed to create PNG image");

	//write to file
	std::ofstream os(name, std::ios::out | std::ios::binary);
	os.write((char*)pPNG_data, pngSize);
	mz_free(pPNG_data);//cleanup memory allocated by png creation
}

int main(int argc, char *argv[]) {
	try {
		//sanity check argument count
		//this should be adjusted in the future to allow batch conversion since there is a good amount of overhead in building the transformer
		if(4 != argc) {
			std::cout << "usage: " << argv[0] << " inputFile nhFile shFile\n";
			std::cout << "\tinputFile - spherical hamrnoics file to write (*.Sht)\n";
			std::cout << "\tnhFile    - location to write north hemisphere image (*.png)\n";
			std::cout << "\tshFile    - location to write south hemisphere image (*.png)\n";
			return EXIT_FAILURE;
		}

		//read in the master spectrum
		emsphinx::MasterSpectra<double> spec;
		spec.read(argv[1]);

		std::cout << spec.getKv() << ' ' << spec.getSig() << '\n';

		//build spherical harmonic transformer and reconstruct on square legendre grid
		const size_t dim = spec.getBw() + (spec.getBw() % 2 == 0 ? 3 : 2);
		emsphinx::square::DiscreteSHT<double> sht(dim, spec.getBw(), emsphinx::square::Layout::Legendre);
		std::vector<double> sph(dim * dim * 2, 0.0);
		sht.synthesize(spec.data(), sph.data(), sph.data() + dim * dim);

		//convert from doubles to 8 bit
		std::pair<double*, double*> minMax = std::minmax_element(sph.data(), sph.data() + sph.size());
		const double vMin = *minMax.first;
		const double delt = *minMax.second - vMin;
		const double fact = 255.0 / delt;
		std::vector<uint8_t> sph8(sph.size());
		std::transform(sph.begin(), sph.end(), sph8.begin(), [fact, vMin](const double& v) {return (uint8_t)std::round((v - vMin) * fact);});

		//write outputs
		writePng(sph8.data()            , dim, dim, 1, argv[2]);
		writePng(sph8.data() + dim * dim, dim, dim, 1, argv[3]);

		//now read in the raw SHT file and print header information
		sht::File file;
		std::ifstream is(argv[1], std::ios::in | std::ios::binary);
		file.read(is);

		//print out header info
		std::cout << "file version " << (int)file.header.fileVersion()[0] << '.' << (int)file.header.fileVersion()[1] << "\n";
		std::cout << "written with software version " << std::string(file.header.softwareVersion(), file.header.softwareVersion()+8) << "\n";
		std::cout << "modality: ";
		switch(file.header.modality()) {
			case sht::Modality::Unknown: std::cout << "Unknown\n"; break;
			case sht::Modality::EBSD   : std::cout << "EBSD   \n"; break;
			case sht::Modality::ECP    : std::cout << "ECP    \n"; break;
			case sht::Modality::TKD    : std::cout << "TKD    \n"; break;
			case sht::Modality::PED    : std::cout << "PED    \n"; break;
			case sht::Modality::Laue   : std::cout << "Laue   \n"; break;
			default                    : std::cout << "invalid\n"; break;
		}
		std::cout << "beam eng: "  << (int) file.header.beamEnergy() << '\n';
		std::cout << "angle 1 : "  << (int) file.header.primaryAngle() << '\n';
		std::cout << "angle 2 : "  << (int) file.header.secondaryAngle() << '\n';
		std::cout << "res : "  << (int) file.header.reservedParam() << '\n';
		std::cout << "notes   : `" << file.header.notes.substr(0, file.header.noteLen()) << "'\n";
		std::cout << "doi     : `" << file.header.doi.substr(0, file.header.doiLen()) << "'\n";
		std::cout << '\n';

		//print out master pattern info
		std::cout << "master pattern composed from " << (int)file.mpData.numXtal() << " crystals with effective sg# " << (int)file.mpData.sgEff() << '\n';
		std::cout << "rotations are " << (char) file.mpData.rotSense() << " with pijk = " << (int) file.mpData.pijk() << '\n';
		std::cout << "simulation data " << file.mpData.simMetaSize() << " bytes from vendor ";
		switch(file.mpData.vendor()) {
			case sht::Vendor::Unknown: std::cout << "Unknown"; break;
			case sht::Vendor::EMsoft : std::cout << "EMsoft" ; break;
			default                  : std::cout << "invalid"; break;
		}
		std::cout << " for modality ";
		switch(file.header.modality()) {
			case sht::Modality::Unknown: std::cout << "Unknown\n"; break;
			case sht::Modality::EBSD   : std::cout << "EBSD   \n"; break;
			case sht::Modality::ECP    : std::cout << "ECP    \n"; break;
			case sht::Modality::TKD    : std::cout << "TKD    \n"; break;
			case sht::Modality::PED    : std::cout << "PED    \n"; break;
			case sht::Modality::Laue   : std::cout << "Laue   \n"; break;
			default                    : std::cout << "invalid\n"; break;
		}
		for(int8_t i = 0; i < file.mpData.numXtal(); i++) {
			const sht::CrystalData& xtal = file.mpData.xtals[i];
			std::cout << "\tsg " << (int) xtal.sgNum() << " setting " << (int) xtal.sgSet() << '\n';
			std::cout << "\t\taxis / cell choice: " << (int) xtal.sgAxis() << " / " << (int) xtal.sgCell() << "\n";
			std::cout << "\t\tadditional origin shift: " << xtal.oriX() << ", " << xtal.oriY() << ", " << xtal.oriZ() << "\n";
			std::cout << "\t\tabc: " << xtal.lat()[0] << ", " << xtal.lat()[1] << ", " << xtal.lat()[2] << "\n";
			std::cout << "\t\tabg: " << xtal.lat()[3] << ", " << xtal.lat()[4] << ", " << xtal.lat()[5] << "\n";
			std::cout << "\t\trot: " << xtal.rot()[0] << ", " << xtal.rot()[1] << ", " << xtal.rot()[2] << ", " << xtal.rot()[3] << "\n";
			std::cout << "\t\twgt: " << xtal.weight() << '\n';
			std::cout << "\t\tfrm: " << '`' << xtal.form.substr(0, xtal.formulaLen  ()) << '\'' << '\n';
			std::cout << "\t\tnam: " << '`' << xtal.name.substr(0, xtal.matNameLen  ()) << '\'' << '\n';
			std::cout << "\t\tsym: " << '`' << xtal.symb.substr(0, xtal.structSymLen()) << '\'' << '\n';
			std::cout << "\t\tref: " << '`' << xtal.refs.substr(0, xtal.refsLen     ()) << '\'' << '\n';
			std::cout << "\t\tnot: " << '`' << xtal.note.substr(0, xtal.noteLen     ()) << '\'' << '\n';
			std::cout << "\t\t" << (int) xtal.numAtoms() << " atoms:\n";

			for(const sht::AtomData& at : xtal.atoms) {
				std::cout << "\t\t\t" << (int) at.atZ() << ": " << at.x() / 24.0f << ' ' << at.y() / 24.0f << ' ' << at.z() / 24.0f << ' ' << at.occ() << ' ' << at.charge() << ' ' << at.debWal() << '\n';
			}

			if(sht::Vendor::EMsoft == file.mpData.vendor  () &&
			   sht::Modality::EBSD == file.header.modality() &&
			   NULL != file.mpData.simul[i].get()) {
				sht::EMsoftED* pED2 = (sht::EMsoftED*) file.mpData.simul[i].get();

				std::cout << "\temVers   : " << std::string(pED2->emsoftVersion()) << '\n';
				std::cout << "\tsigStart : " << pED2->sigStart () << '\n';
				std::cout << "\tsigEnd   : " << pED2->sigEnd   () << '\n';
				std::cout << "\tsigStep  : " << pED2->sigStep  () << '\n';
				std::cout << "\tomega    : " << pED2->omega    () << '\n';
				std::cout << "\tkeV      : " << pED2->keV      () << '\n';
				std::cout << "\teHistMin : " << pED2->eHistMin () << '\n';
				std::cout << "\teBinSize : " << pED2->eBinSize () << '\n';
				std::cout << "\tdepthMax : " << pED2->depthMax () << '\n';
				std::cout << "\tdepthStep: " << pED2->depthStep() << '\n';
				std::cout << "\tthickness: " << pED2->thickness() << '\n';
				std::cout << "\ttotNumEl : " << pED2->totNumEl () << '\n';
				std::cout << "\tnumSx    : " << pED2->numSx    () << '\n';
				std::cout << "\tc1       : " << pED2->c1       () << '\n';
				std::cout << "\tc2       : " << pED2->c2       () << '\n';
				std::cout << "\tc3       : " << pED2->c3       () << '\n';
				std::cout << "\tsigDbDiff: " << pED2->sigDbDiff() << '\n';
				std::cout << "\tdMin     : " << pED2->dMin     () << '\n';
				std::cout << "\tnumPx    : " << pED2->numPx    () << '\n';
				std::cout << "\tlatGridType: ";
				switch( pED2->latGridType()) {
					case 1 : std::cout << "square lambert\n";break;
					case 2 : std::cout << "square legendre\n";break;
					default: std::cout << "unknown\n";break;
				} 
			}
		}

	} catch(std::exception& e) {
		std::cout << e.what() << '\n';
	}

	return EXIT_SUCCESS;
}
