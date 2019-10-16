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

#include "idx/master.hpp"
#include "sht_file.hpp"
/*

#define MINIZ_NO_STDIO
#define MINIZ_NO_TIME
#define MINIZ_NO_ZLIB_APIS
#include "miniz/miniz.c"
//@brief   : extract the crystal data from an EMsoft master pattern file
//@param fn: name of file to extract from
sht::CrystalData extractXtal(char * fn) {
	//read scalar data
	int v;
	sht::CrystalData xtal;
	H5::H5File file = H5::H5File(fn, H5F_ACC_RDONLY);//read only access
	file.openGroup("CrystalData").openDataSet("SpaceGroupNumber" ).read(&v        , H5::PredType::NATIVE_INT  ); xtal.sgNum   () = (uint8_t) v;
	file.openGroup("CrystalData").openDataSet("SpaceGroupSetting").read(&v        , H5::PredType::NATIVE_INT  ); xtal.sgSet   () = (uint8_t) v;
	file.openGroup("CrystalData").openDataSet("LatticeParameters").read(xtal.lat(), H5::PredType::NATIVE_FLOAT);
	file.openGroup("CrystalData").openDataSet("Natomtypes"       ).read(&v        , H5::PredType::NATIVE_INT  ); xtal.numAtoms() = (uint16_t) v;

	//read atom data
	std::vector<int> atomTypes(v);
	std::vector<float> atomData(v * 5);
	file.openGroup("CrystalData").openDataSet("AtomData" ).read(atomData .data(), H5::PredType::NATIVE_FLOAT);
	file.openGroup("CrystalData").openDataSet("Atomtypes").read(atomTypes.data(), H5::PredType::NATIVE_INT  );


	xtal.sgAxis() = sht::CrystalData::Axis::Default;
	xtal.sgCell() = sht::CrystalData::Cell::Default;
	xtal.oriX() = xtal.oriY() = xtal.oriZ() = 0.0f;
	xtal.rot()[0] = 1.0f;
	xtal.rot()[1] = xtal.rot()[2] = xtal.rot()[3] = 0.0f;
	xtal.weight() = 1.0f;
	xtal.atoms.resize(v);

	for(int i = 0; i < v; i++) {
		xtal.atoms[i].x     () = atomData [i*5 + 0];
		xtal.atoms[i].y     () = atomData [i*5 + 1];
		xtal.atoms[i].z     () = atomData [i*5 + 2];
		xtal.atoms[i].occ   () = atomData [i*5 + 3];
		xtal.atoms[i].charge() = 0.0f;
		xtal.atoms[i].debWal() = atomData [i*5 + 4];
		xtal.atoms[i].resFp () = 0.0f;
		xtal.atoms[i].atZ   () = atomTypes[i      ];
	}
	return xtal;
}
*/

int main(int argc, char *argv[]) {

	//sanity check argument count
	//this should be adjusted in the future to allow batch conversion since there is a good amount of overhead in building the transformer
	if(3 != argc) {
		std::cout << "usage: " << argv[0] << " inputFile outputFile\n";
		std::cout << "\tinputFile  - master pattern to read (*.h5)\n";
		std::cout << "\toutputFile - spherical hamrnoics file to write (*.spx)\n";
		return EXIT_FAILURE;
	}

	try {
		//read in the master pattern and convert to spectrum
		const size_t bw = 384;
		const bool nrm = true;
		emsphinx::MasterSpectra<double> spec(emsphinx::MasterPattern<double>(argv[1]), bw, nrm);

		//read in crystal data
		float lat[6];
		int32_t sgN, sgS, nAt;
		H5::H5File file = H5::H5File(argv[1], H5F_ACC_RDONLY);//read only access
		file.openGroup("CrystalData").openDataSet("SpaceGroupNumber" ).read(&sgN, H5::PredType::NATIVE_INT32);
		file.openGroup("CrystalData").openDataSet("SpaceGroupSetting").read(&sgS, H5::PredType::NATIVE_INT32);
		file.openGroup("CrystalData").openDataSet("LatticeParameters").read( lat, H5::PredType::NATIVE_FLOAT);
		file.openGroup("CrystalData").openDataSet("Natomtypes"       ).read(&nAt, H5::PredType::NATIVE_INT32);
		std::vector<int32_t> aTy(nAt);
		std::vector<float> aCd(nAt * 5);
		file.openGroup("CrystalData").openDataSet("AtomData" ).read(aCd .data(), H5::PredType::NATIVE_FLOAT);
		file.openGroup("CrystalData").openDataSet("Atomtypes").read(aTy.data(), H5::PredType::NATIVE_INT32);

		//read in simulation data
		float fprm[15];
		int32_t iprm[5];

		file.openGroup("NMLparameters").openGroup("MCCLNameList"      ).openDataSet("sig"       ).read(fprm +  0, H5::PredType::NATIVE_FLOAT);
		fprm[1] = NAN;//sig end
		fprm[2] = NAN;//sig step
		file.openGroup("NMLparameters").openGroup("MCCLNameList"      ).openDataSet("omega"     ).read(fprm +  3, H5::PredType::NATIVE_FLOAT);
		file.openGroup("NMLparameters").openGroup("MCCLNameList"      ).openDataSet("EkeV"      ).read(fprm +  4, H5::PredType::NATIVE_FLOAT);
		file.openGroup("NMLparameters").openGroup("MCCLNameList"      ).openDataSet("Ehistmin"  ).read(fprm +  5, H5::PredType::NATIVE_FLOAT);
		file.openGroup("NMLparameters").openGroup("MCCLNameList"      ).openDataSet("Ebinsize"  ).read(fprm +  6, H5::PredType::NATIVE_FLOAT);
		file.openGroup("NMLparameters").openGroup("MCCLNameList"      ).openDataSet("depthmax"  ).read(fprm +  7, H5::PredType::NATIVE_FLOAT);
		file.openGroup("NMLparameters").openGroup("MCCLNameList"      ).openDataSet("depthstep" ).read(fprm +  8, H5::PredType::NATIVE_FLOAT);
		fprm[9] = std::numeric_limits<float>::infinity();//thickness
		file.openGroup("NMLparameters").openGroup("BetheList"         ).openDataSet("c1"        ).read(fprm + 10, H5::PredType::NATIVE_FLOAT);
		file.openGroup("NMLparameters").openGroup("BetheList"         ).openDataSet("c2"        ).read(fprm + 11, H5::PredType::NATIVE_FLOAT);
		file.openGroup("NMLparameters").openGroup("BetheList"         ).openDataSet("c3"        ).read(fprm + 12, H5::PredType::NATIVE_FLOAT);
		file.openGroup("NMLparameters").openGroup("BetheList"         ).openDataSet("sgdbdiff"  ).read(fprm + 13, H5::PredType::NATIVE_FLOAT);
		file.openGroup("NMLparameters").openGroup("EBSDMasterNameList").openDataSet("dmin"      ).read(fprm + 14, H5::PredType::NATIVE_FLOAT);

		file.openGroup("NMLparameters").openGroup("MCCLNameList"      ).openDataSet("totnum_el" ).read(iprm +  0, H5::PredType::NATIVE_INT32);
		file.openGroup("NMLparameters").openGroup("MCCLNameList"      ).openDataSet("multiplier").read(iprm +  1, H5::PredType::NATIVE_INT32);
		file.openGroup("NMLparameters").openGroup("MCCLNameList"      ).openDataSet("numsx"     ).read(iprm +  2, H5::PredType::NATIVE_INT32);
		file.openGroup("NMLparameters").openGroup("EBSDMasterNameList").openDataSet("npx"       ).read(iprm +  3, H5::PredType::NATIVE_INT32);
		iprm[4] = 1;//lattitude grid type

		//build compression flags
		// int8_t flg[2];
		// flg[0] = spec.nFold();
		// flg[1] = (spec.invSym() ? 1 : 0) + (spec.mirror() ? 2 : 0) + (4 * spec.phase().pg.mmType());

		//@brief    : write a file using EMsoft style EBSD data
		//@prief fn : file name to write
		//@prief nt : notes string
		//@param sgN: space group number [1,230]
		//@param sgS: space group setting [1,2]
		//@param nAt: number of atoms
		//@param aTy: atom types (nAt atomic numbers)
		//@param aCd: atom coordinates, (nAt * 5 floats {x, y, z, occupancy, Debye-Waller in nm^2})
		//@param lat: lattice parameters {a, b, a, alpha, beta, gamma} (in nm / degree)
		//@param fprm: floating point parameters (float32 EMsoftED parameters in order)
		//@param iprm: integer parameters {# electrons, electron multiplier, numsx, npx, latgridtype}
		//@param bw : bandwidth
		//@param flg: symmetry flags {zRot, mirInv}
		//@param alm: actual harmonics (uncompressed format)
		std::string notes("notes string");
		std::string doi("doi string");
		sht::File::EMsoftEBSD(argv[2], notes.c_str(), doi.c_str(), sgN, sgS, nAt, aTy.data(), aCd.data(), lat, fprm, iprm, (int32_t)bw, (double*)spec.data());




/*



		sht::File file;

		//build header info
		file.header.modality() = sht::Modality::EBSD;
		file.header.setDoi("doistr");
		file.header.setNotes("file notes");

		//vendor + simullen set by file

		file.header.beamEnergy  () = (float) spec.getKv ();
		file.header.primaryAngle() = (float) spec.getSig();

		//build material data
		file.material.xtals.push_back(extractXtal(argv[1]));
		file.material.numXtal() = 1;
		file.material.sgEff  () = file.material.xtals.front().sgNum();

		//build simulation meta data
				// std::unique_ptr<SimulationData> simulMeta;//can be null if header.simDataSize() == 0

		//extract harmonics
		const int32_t nHarm = sht::HarmonicsData::NumHarm(bw, (int8_t)spec.nFold(), spec.invSym(), spec.mirror());
		file.harmonics.bw     () = bw;
		file.harmonics.zRot   () = (int8_t) spec.nFold();
		file.harmonics.mirInv () = (spec.invSym() ? 1 : 0) | (spec.mirror() ? 2 : 0);
		file.harmonics.doubCnt() = nHarm * 2;
		file.harmonics.alm.resize(nHarm * 2);
		double * pHrm = file.harmonics.alm.data();
		for(size_t m = 0; m < bw; m++) {
			std::complex<double> * const row = spec.data() + bw * m;
			if(spec.nFold() < 2 ? false : 0 != m % spec.nFold()) continue;//systemic zeros
			for(size_t l = m; l < bw; l++) {
				if( (spec.invSym() && 0 != l % 2) || (spec.mirror() && 0 != (l + m) % 2) ) continue;
				*pHrm++ = row[l].real();
				*pHrm++ = row[l].imag();
			}
		}

		//if we made it this far everything needed was parsed, write out the result
		std::ofstream os(argv[2], std::ios::out | std::ios::binary);
		file.write(os);
		*/
		return EXIT_SUCCESS;
	} catch (std::exception& e) {
		std::cout << e.what() << '\n';
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;

}
