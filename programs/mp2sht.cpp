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

		//save header data
		float fprm[25];
		int32_t iprm[11];
		iprm[1] = (int32_t) sht::Modality::EBSD;
		iprm[2] = bw;
		fprm[0] = (float) spec.getKv ();
		fprm[1] = (float) spec.getSig();
		fprm[2] = 0.0f;
		fprm[3] = 0.0f;

		//read in crystal data
		float lat[6];
		H5::H5File file = H5::H5File(argv[1], H5F_ACC_RDONLY);//read only access
		file.openGroup("CrystalData").openDataSet("SpaceGroupNumber" ).read(iprm + 3, H5::PredType::NATIVE_INT32); iprm[0] = iprm[3];//effective space group
		file.openGroup("CrystalData").openDataSet("SpaceGroupSetting").read(iprm + 4, H5::PredType::NATIVE_INT32);
		file.openGroup("CrystalData").openDataSet("LatticeParameters").read(fprm + 4, H5::PredType::NATIVE_FLOAT);
		file.openGroup("CrystalData").openDataSet("Natomtypes"       ).read(iprm + 5, H5::PredType::NATIVE_INT32);
		std::vector<int32_t> aTy(iprm[5]);
		std::vector<float> aCd(iprm[5] * 5);
		file.openGroup("CrystalData").openDataSet("AtomData" ).read(aCd .data(), H5::PredType::NATIVE_FLOAT);
		file.openGroup("CrystalData").openDataSet("Atomtypes").read(aTy.data(), H5::PredType::NATIVE_INT32);

  		//read in simulation data
		file.openGroup("NMLparameters").openGroup("MCCLNameList"      ).openDataSet("sig"       ).read(fprm + 10, H5::PredType::NATIVE_FLOAT);
		fprm[11] = NAN;//sig end
		fprm[12] = NAN;//sig step
		file.openGroup("NMLparameters").openGroup("MCCLNameList"      ).openDataSet("omega"     ).read(fprm + 13, H5::PredType::NATIVE_FLOAT);
		file.openGroup("NMLparameters").openGroup("MCCLNameList"      ).openDataSet("EkeV"      ).read(fprm + 14, H5::PredType::NATIVE_FLOAT);
		file.openGroup("NMLparameters").openGroup("MCCLNameList"      ).openDataSet("Ehistmin"  ).read(fprm + 15, H5::PredType::NATIVE_FLOAT);
		file.openGroup("NMLparameters").openGroup("MCCLNameList"      ).openDataSet("Ebinsize"  ).read(fprm + 16, H5::PredType::NATIVE_FLOAT);
		file.openGroup("NMLparameters").openGroup("MCCLNameList"      ).openDataSet("depthmax"  ).read(fprm + 17, H5::PredType::NATIVE_FLOAT);
		file.openGroup("NMLparameters").openGroup("MCCLNameList"      ).openDataSet("depthstep" ).read(fprm + 18, H5::PredType::NATIVE_FLOAT);
		fprm[19] = std::numeric_limits<float>::infinity();//thickness
		file.openGroup("NMLparameters").openGroup("BetheList"         ).openDataSet("c1"        ).read(fprm + 20, H5::PredType::NATIVE_FLOAT);
		file.openGroup("NMLparameters").openGroup("BetheList"         ).openDataSet("c2"        ).read(fprm + 21, H5::PredType::NATIVE_FLOAT);
		file.openGroup("NMLparameters").openGroup("BetheList"         ).openDataSet("c3"        ).read(fprm + 22, H5::PredType::NATIVE_FLOAT);
		file.openGroup("NMLparameters").openGroup("BetheList"         ).openDataSet("sgdbdiff"  ).read(fprm + 23, H5::PredType::NATIVE_FLOAT);
		file.openGroup("NMLparameters").openGroup("EBSDMasterNameList").openDataSet("dmin"      ).read(fprm + 24, H5::PredType::NATIVE_FLOAT);

		file.openGroup("NMLparameters").openGroup("MCCLNameList"      ).openDataSet("totnum_el" ).read(iprm +  6, H5::PredType::NATIVE_INT32);
		file.openGroup("NMLparameters").openGroup("MCCLNameList"      ).openDataSet("multiplier").read(iprm +  7, H5::PredType::NATIVE_INT32);
		file.openGroup("NMLparameters").openGroup("MCCLNameList"      ).openDataSet("numsx"     ).read(iprm +  8, H5::PredType::NATIVE_INT32);
		file.openGroup("NMLparameters").openGroup("EBSDMasterNameList").openDataSet("npx"       ).read(iprm +  9, H5::PredType::NATIVE_INT32);
		iprm[10] = 1;//lattitude grid type

		//build set of elements in order for formula estimate (should really include multiplicity)
		std::set<int32_t> atomType(aTy.begin(), aTy.end());
		std::string form;
		static const std::vector<std::string> AtSyb = {
			"H" , "He", "Li", "Be", "B" , "C" , "N" , "O" , "F" , "Ne", "Na", "Mg", "Al", "Si", "P" , "S" ,
			"Cl", "Ar", "K" , "Ca", "Sc", "Ti", "V" , "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge",
			"As", "Se", "Br", "Kr", "Rb", "Sr", "Y" , "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
			"In", "Sn", "Sb", "Te", "I" , "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd",
			"Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W" , "Re", "Os", "Ir", "Pt", "Au", "Hg",
			"Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U" , "Np", "Pu", "Am", "Cm",
			"Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn",
			"Nh", "Fl", "Mc", "Lv", "Ts", "Og", 
		};
		for(const int32_t& z : atomType) form += AtSyb[z-1];

		//assemble file
		sht::File shtFile;
		std::string doi = "https://doi.org/10.1016/j.ultramic.2019.112841";
		std::string note = "created with mp2sht";
		char emVers[8] = "5_0_0_0";
		std::string cprm = form + '\0';
		cprm += "\0";//material name e.g. gamma
		cprm += "\0";//structure symbol e.g. L1_2
		cprm += "\0";//references
		cprm += "\0";//note
		shtFile.initFileEMsoft(iprm, fprm, doi.c_str(), note.c_str(), (double*)spec.data());//initialize header + harmonics
		shtFile.addDataEMsoft(iprm + 3, fprm + 4, aTy.data(), aCd.data(), emVers, cprm.c_str());//add crystal + simulation data
		std::ofstream os(argv[2], std::ios::out | std::ios::binary);
		shtFile.write(os);

		return EXIT_SUCCESS;
	} catch (std::exception& e) {
		std::cout << e.what() << '\n';
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;

}
