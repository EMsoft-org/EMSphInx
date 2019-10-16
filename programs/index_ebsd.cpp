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

#include <iomanip>
#include <atomic>

#include "util/threadpool.hpp"
#include "util/timer.hpp"
#include "util/nml.hpp"
#include "util/sysnames.hpp"
#include "modality/ebsd/pattern.hpp"

#include "modality/ebsd/nml.hpp"

#include "modality/ebsd/idx.hpp"

int main(int argc, char *argv[]) {

try {

	typedef double Real;//should we use float, double, or long double for calculations?

	////////////////////////////////////////////////////////////////////////
	//                          Parse Arguments                           //
	////////////////////////////////////////////////////////////////////////

	//check argument count
	if(2 != argc) {
		std::cout << "useage: ";
		std::cout << "\tindex using a nml file : " << argv[0] << " input.nml\n";
		std::cout << "\tgenerate a template nml: " << argv[0] << " -t\n";
		return EXIT_FAILURE;
	}

	//check for template request
	const std::string nmlName(argv[1]);
	if(0 == nmlName.compare("-t")) {
		std::ofstream os(std::string("IndexEBSD") + ".nml");//create programname.nml (in the future argv[0] could be used with std::filesystem to remove full path)
		emsphinx::ebsd::Namelist nml;
		nml.defaults();
		os << nml.to_string();
		return EXIT_SUCCESS;
	}

	//read nml and parse
	emsphinx::ebsd::Namelist nml;
	{
		std::ifstream is(nmlName);
		std::string str((std::istreambuf_iterator<char>(is)), std::istreambuf_iterator<char>());
		std::string warning = nml.from_string(str);
		if(!warning.empty()) {
			std::cout << "\n * * * * * * warning: some namelist parameters weren't used: " << warning << " * * * * * * \n" << std::endl;
		}
	}

	////////////////////////////////////////////////////////////////////////
	//                    Read Inputs / Build Indexers                    //
	////////////////////////////////////////////////////////////////////////

	emsphinx::ebsd::IndexingData<Real> idxData(nml);

	////////////////////////////////////////////////////////////////////////
	//                             Print Info                             //
	////////////////////////////////////////////////////////////////////////

	//lets print some information before indexing gets started
	std::cout << '\n' << std::boolalpha;
	std::cout << "Running program \"" << argv[0] << "\"\n";
	std::cout << "\tCompiled From        : " << __FILE__ << '\n';
	std::cout << "\tGit Branch           : " << emsphinx::GitBranch << '\n';
	std::cout << "\tCommit Hash          : " << emsphinx::GitHash << '\n';
	std::cout << "\tVersion String       : " << emsphinx::Version << '\n';
	std::cout << '\n';
	std::cout << "Geometry\n";
	std::cout << "\tSample Tilt          : " << idxData.phases.front().getSig() << " degrees\n";
	std::cout << "\tScintillator Distance: " << idxData.geom.sDst << " microns\n";
	std::cout << "\tCamera Tilt          : " << idxData.geom.dTlt << " degrees\n";
	std::cout << "\tCamera               : " << nml.patDims[0] << " x " << nml.patDims[1] << " with " << idxData.geom.pX << " micron pixels\n";
	std::cout << "\tPattern Center       : " << idxData.geom.cX << ", " << idxData.geom.cY << " fractional pixels\n";
	std::cout << "\tCircular Mask        : " << idxData.geom.circ << "\n";
	std::cout << "\tVertical Flip        : " << idxData.geom.flip << "\n";
	std::cout << "\n";
	std::cout << "Indexing patterns from \"" << nml.patFile;
	if(!nml.patName.empty()) std::cout << ":/" << nml.patName;
	std::cout << "\"\n";
	std::cout << "\tScan Dimensions      : " << nml.scanDims[0] << " x " << nml.scanDims[1] << " pixels\n";
	std::cout << "\tScan Resolution      : " << nml.scanSteps[0] << " x " << nml.scanSteps[1]  << " microns\n";
	std::cout << "\tPattern bitdepth     : " << idxData.pat->pixBytes() * 8 << '\n';
	std::cout << "\tTotal Patterns       : " << idxData.pat->numPat() << '\n';
	std::cout << "\tAHE grid points      : " << nml.nRegions << '\n';
	std::cout << "\n";
	std::cout << "Against master pattern" << (idxData.phases.size() > 1 ? "s" : "") << ":\n";
	for(size_t i = 0; i < nml.masterFiles.size(); i++) {
		std::cout << "\tFile Name            : " << nml.masterFiles[i] << '\n';
		std::cout << "\tPoint Group          : " << idxData.phases[i].pointGroup().name() << '\n';
		std::cout << "\tZ Rotational Symmetry: " << (int)idxData.phases[i].pointGroup().zRot() << '\n';
		std::cout << "\tEquatorial Mirror    : " << (idxData.phases[i].pointGroup().zMirror() ? "yes" : "no") << '\n';
	}
	std::cout << "\n";
	std::cout << "Indexing with\n";
	std::cout << "\tBandwidth            : " << nml.bw << '\n';
	std::cout << "\tSide Length          : " << fft::fastSize(uint32_t(2 * nml.bw - 1)) << '\n';
	std::cout << "\tROI Mask             : " << ( !nml.roi.hasShape() ? "entire scan" : nml.roi.to_string() )<< '\n';
	std::cout << "\tThread Count         : " << idxData.threadCount << '\n';
	std::cout << "\tBatch Size           : " << nml.batchSize << '\n';
	std::cout.flush();

	// const std::pair<size_t, size_t> scaleSize = idxData.indexers.front()->rescaledSize();
	// std::cout << "\tPatterns Resized To  : " << scaleSize.first << 'x' << scaleSize.second << '\n';
	std::cout << std::endl;

	////////////////////////////////////////////////////////////////////////
	//                            Do Indexing                             //
	////////////////////////////////////////////////////////////////////////

	//build thread pool and get start time
	ThreadPool pool(idxData.threadCount);//pool
	time_t tmStart = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());

	//queue parallel indexing
	Timer t;
	size_t batches = idxData.pat->numPat() / nml.batchSize;//how many batches are needed
	if(batches * nml.batchSize < idxData.pat->numPat()) ++batches;//extra batch for leftovers
	for(size_t i = 0; i < batches; i++) {//loop over batches
		const size_t start = i * nml.batchSize;//first pattern
		const size_t end = std::min(start + nml.batchSize, idxData.pat->numPat());//last pattern
		pool.schedule(std::bind(idxData.workItem, std::placeholders::_1));//queue indexing
	}

	//wait for indexing to complete
	const std::chrono::milliseconds uptFreq(1000);//milliseconds between updates
	while(!pool.waitAll(uptFreq)) {
		//get the time elapsed since we started (without resetting the reference point)
		const double elapsed = t.poll(false);
		const double percent = double(idxData.idxCtr) / idxData.numIdx;
		const double rate = elapsed / percent;//estimate seconds per %
		const double remaining = rate * (1.0 - percent);//estimate seconds left

		//print update
		std::cout << '\r';
		Timer::PrintSeconds(elapsed  , std::cout);
		std::cout << " elapsed, " << std::fixed << std::setprecision(1) << percent * 100          << "% complete, ";
		Timer::PrintSeconds(remaining, std::cout);
		std::cout << " remaining   ";
		std::cout.flush();
	}

	const double total = t.poll();
	time_t tmEnd = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
	std::cout << '\n' << total << "s to index (" << double(idxData.numIdx) / total << " pat/s)\n";

	////////////////////////////////////////////////////////////////////////
	//                            Save Outputs                            //
	////////////////////////////////////////////////////////////////////////

	idxData.save(tmStart, tmEnd, total);

	//done
	return 0;

} catch (std::exception& e) {
	std::cout << e.what() << '\n';
	return EXIT_FAILURE;
}

}
