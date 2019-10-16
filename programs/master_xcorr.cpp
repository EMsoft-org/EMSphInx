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

#include "util/timer.hpp"
#include "idx/master.hpp"
#include "sht/sht_xcorr.hpp"
#include "xtal/diagram.hpp"
#include "constants.hpp"

int main(int argc, char *argv[]) {
	//parse arguments (file names / bandwidth)
	if(4 != argc && 5 != argc) {
		std::cout << "usage: " << argv[0] << "scanFile bandWidth cutoff masterFile1 [masterFile2]\n";
		std::cout << "\tbandWidth  : bandwidth for cross correlation (2*bw-1 should be product of small primes)\n";
		std::cout << "\t             88, 95, 113, 123, 158, 172, 188, 203, 221, 263, and 284 are reasonable values\n";
		std::cout << "\tcutoff     : cutoff for peak consideration [0,1] (relative to maximum cross correlation)\n";
		std::cout << "\tmasterFile1: name of first master pattern file (e.g. Ni.h5)\n";
		std::cout << "\tmasterFile2: name of second master pattern file (e.g. Ni.h5)\n";
		std::cout << "\nnote: only symmetry of first pattern will be used for to improve calculation speed\n";
		return EXIT_FAILURE;
	}
	const uint16_t bw = (uint16_t)std::strtoul(argv[1], NULL, 0);
	const double cutoff = atof(argv[2]);
	if(cutoff < 0.0 || cutoff > 1.0) throw std::runtime_error("cutoff must lie in [0,1]");
	std::string masterFile1(            argv[3]              );
	std::string masterFile2(argc == 5 ? argv[4] : masterFile1);

	//define some values we may want to parse in the future
	std::string outputFile("pseudo_sym.h5");
	const double factor = 0.95;//how far below the target cutoff should we search to allow for off grid peaks

	//sanity check bandwidth (anything is technically ok but it is probably a good idea to restrict the range)
	if(bw < 53 ) throw std::runtime_error("bandwidth too small");
	if(bw > 313) throw std::runtime_error("bandwidth too large");

	//read master pattern(s) and compute their SHT
	Timer t;
	emsphinx::MasterPattern<double> mp1(masterFile1);//read master pattern
	emsphinx::MasterPattern<double> mp2(masterFile2);//read master pattern
	std::cout << t.poll() << "s to read master patterns" << std::endl;
	emsphinx::MasterSpectra<double> p1(mp1, bw, true);//compute harmonics of normed master pattern
	emsphinx::MasterSpectra<double> p2(mp2, bw, true);//compute harmonics of normed master pattern
	std::cout << t.poll() << "s to compute SHT of patterns" << std::endl;
	p1.removeDC();//it should already be almost zero from nrm == true in construction
	p2.removeDC();//it should already be almost zero from nrm == true in construction

	//compute master pattern cross correlation
	double eu[3];//we need somewhere to write peak correlation (hopefully it is 0...)
	emsphinx::sphere::Correlator<double> s2Corr(bw);
	std::cout << t.poll() << "s to build correlator" << std::endl;
	s2Corr.correlate(p1.data(), p2.data(), p1.mirror(), p1.nFold(), eu, false);//correlate
	std::cout << t.poll() << "s to correlate\n";

	//find peak intensity (no rotation for auto-correlation)
	const size_t sl = bw * 2 - 1;//side length of euler space grid
	if(masterFile1 == masterFile2) {//auto correlation
		const size_t idxIdent = (bw-1) * sl * sl + (bw / 2) * sl + bw / 2;//no rotation
		s2Corr.indexEuler(idxIdent, eu);
	}
	const double vMax = s2Corr.refinePeak(p1.data(), p2.data(), p1.mirror(), p1.nFold(), eu);
	std::cout << "maximum intensity: " << vMax << std::endl;

	//now find all local maxima
	struct Maxima {
		double intensity;
		size_t index;
		xtal::Quat<double> qu;
		inline bool operator<(const Maxima& m) const {return intensity > m.intensity;}//reverse sorting order
	};
	Maxima mx;
	double nh[3][3][3];
	std::vector<Maxima> maxima;
	const double vMin = vMax * cutoff * factor;//don't even both checking neighborhoods around pixels that aren't bright enough

	//loop over cross correlation grid searching for local maxima
	fft::vector<double> xc(s2Corr.getXC());
	for(size_t k = 0; k < bw; k++) {//loop over slices
		for(size_t n = 0; n < sl; n++) {//loop over rows
			for(size_t m = 0; m < sl; m++) {//loop over columns
				const size_t idx = k * sl * sl + n * sl + m;//compute vectorized index
				if(xc[idx] >= vMin) {//check if this voxel is bright enough to consider
					s2Corr.extractNeighborhood<1>(idx, nh);//get 3x3x3 neighborhood surrounding point
					if(nh[1][1][1] >= nh[0][0][0] &&//check if this is a local maxima
					   nh[1][1][1] >= nh[0][0][1] &&
					   nh[1][1][1] >= nh[0][0][2] &&
					   nh[1][1][1] >= nh[0][1][0] &&
					   nh[1][1][1] >= nh[0][1][1] &&
					   nh[1][1][1] >= nh[0][1][2] &&
					   nh[1][1][1] >= nh[0][2][0] &&
					   nh[1][1][1] >= nh[0][2][1] &&
					   nh[1][1][1] >= nh[0][2][2] &&
					   nh[1][1][1] >= nh[1][0][0] &&
					   nh[1][1][1] >= nh[1][0][1] &&
					   nh[1][1][1] >= nh[1][0][2] &&
					   nh[1][1][1] >= nh[1][1][0] &&
					   nh[1][1][1] >= nh[1][1][2] &&
					   nh[1][1][1] >= nh[1][2][0] &&
					   nh[1][1][1] >= nh[1][2][1] &&
					   nh[1][1][1] >= nh[1][2][2] &&
					   nh[1][1][1] >= nh[2][0][0] &&
					   nh[1][1][1] >= nh[2][0][1] &&
					   nh[1][1][1] >= nh[2][0][2] &&
					   nh[1][1][1] >= nh[2][1][0] &&
					   nh[1][1][1] >= nh[2][1][1] &&
					   nh[1][1][1] >= nh[2][1][2] &&
					   nh[1][1][1] >= nh[2][2][0] &&
					   nh[1][1][1] >= nh[2][2][1] &&
					   nh[1][1][1] >= nh[2][2][2]
					) {
						//build maxima
						mx.index = idx;
						mx.intensity = nh[1][1][1] / vMax;
						s2Corr.indexEuler(idx, eu);
						xtal::zyz2qu(eu, mx.qu.data());

						//find the closest existing maxima
						size_t iNear = 0;//index of nearest maxima
						double iDist = 360.;//distance to nearest maxima
						for(size_t i = 0; i < maxima.size(); i++) {//loop over maxima
							const double dot = std::min<double>(1, std::fabs( mx.qu.dot(maxima[i].qu) ) );
							const double angle = std::acos(dot) * xtal::Constants<double>::rd2dg;
							if(angle < iDist) {//is this the closest maxima so far?
								iNear = i;
								iDist = angle;
							}
						}

						//make sure this isn't the same as an orientation we've already found
						if(iDist < 2.0) {//are we within 2 degrees (this is extremely arbitrary)
							if(mx.intensity > maxima[iNear].intensity) {
								maxima[iNear] = mx;//keep the better maxima
							}
						} else {
							maxima.push_back(mx);//add a new maxima
						}
					}
				}
			}
		}
	}

	//now sort local maxima by intensity
	std::sort(maxima.begin(), maxima.end());
	std::cout << t.poll() << "s to extract and sort " << maxima.size() << " local maxima\n";
	if(maxima.empty()) {
		std::cout << "no local maxima found in autocorrelation\n";
		return EXIT_FAILURE;
	}

	//work through maxima refining
	for(size_t i = 0; i < maxima.size(); i++) {
		s2Corr.indexEuler(maxima[i].index, eu);
		maxima[i].intensity = s2Corr.refinePeak(p1.data(), p2.data(), p1.mirror(), p1.nFold(), eu) / vMax;
		xtal::zyz2qu(eu, maxima[i].qu.data());
	}
	std::cout << t.poll() << "s to refine\n";

	//now resort and print
	std::sort(maxima.begin(), maxima.end());
	std::cout << std::fixed << std::setprecision(4);
	for(size_t i = 0; i < maxima.size(); i++) {
		if(maxima[i].intensity < cutoff) break;
		std::cout << maxima[i].intensity << ' ' << maxima[i].qu << '\n';
	}

	//normalize entire dataset by brightest pixel
	// std::for_each(xc.begin(), xc.end(), [vMax](double& v){v /= vMax;});

	//create file to write result to
	hsize_t dims[3] = {bw, sl, sl};
	H5::H5File file = H5::H5File(outputFile, H5F_ACC_TRUNC);//overwrite existing if needed
	file.createDataSet("Cross Correlation", H5::PredType::NATIVE_DOUBLE, H5::DataSpace(3, dims)).write(xc.data(), H5::PredType::NATIVE_DOUBLE);

	//now let's make and xdmf file for easy visualization
	std::ofstream os("pseudo_sym.xdmf");
	const double res = 360.0 / sl;//angular resolution in degrees
	os << "<?xml version=\"1.0\"?>\n";
	os << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\"[]>\n";
	os << "<Xdmf xmlns:xi=\"http://www.w3.org/2003/XInclude\" Version=\"2.2\">\n";
	os << " <Domain>\n";
	os << "  <Grid Name=\"MasterPatternMesh\" GridType=\"Uniform\">\n";
	os << "   <Topology TopologyType=\"3DCoRectMesh\" Dimensions=\"" << bw << ' ' << sl << ' ' << sl << "\"></Topology>\n";
	os << "   <Geometry Type=\"ORIGIN_DXDYDZ\">\n";
	os << "     <DataItem Format=\"XML\" Dimensions=\"3\">0 0 0</DataItem>\n";
	os << "     <DataItem Format=\"XML\" Dimensions=\"3\">" << res << ' ' << res << ' ' << res << "</DataItem>\n";
	os << "   </Geometry>\n";
	os << "   <Attribute Name=\"Correlation\" AttributeType=\"Scalar\" Center=\"Cell\">\n";
	os << "    <DataItem Format=\"HDF\" NumberType=\"Float\" Precision=\"8\" Dimensions=\"" << bw << ' ' << sl << ' ' << sl << "\">\n";
	os << "     " << outputFile << ":/Cross Correlation\n";
	os << "    </DataItem>\n";
	os << "   </Attribute>\n";
	os << "  </Grid>\n";
	os << " </Domain>\n";
	os << "</Xdmf>\n";

	//next create diagram
	svg::Color c(0, 0, 0);
	xtal::Diagram diag(mp1, c);
	diag.getHemi().write("true.svg");//save true symmetry

	//get rotational symmetry operators of the point group
	const size_t numOp                   =                             mp1.pointGroup().numRotOps     ();
	xtal::Quat<double> const * const ops = (xtal::Quat<double> const *)mp1.pointGroup().rotOps<double>();

	//loop over detected psuedosymmetry
	const double degCut = 0.999;//cosine(same rotation) ~2.5 degrees
	const double rotCut = 0.05;//cutoff for something to be considered an n fold axis
	for(size_t i = 0; i < maxima.size(); i++) {
		if(maxima[i].intensity < cutoff) break;

		//make sure this isn't a symmetry operator
		bool match = false;
		for(size_t j = 0; j < numOp; j++) {
			const double dot = std::fabs(ops[j].dot(maxima[i].qu));
			if(dot > degCut) {//~2.5 degrees
				match = true;
				break;
			}
		}

		//this isn't a symmetry operator
		if(!match) {
			//make sure we haven't already added this operator
			for(size_t j = 0; j < i; j++) {
				const double dot = std::fabs(maxima[i].qu.dot(maxima[j].qu));
				if(dot > degCut) {//~2.5 degrees
					match = true;
					break;
				}
			}

			if(!match) {//if we made it this far this is a new pseudo symmetry operator
				c.rgb[0] = c.rgb[1] = maxima[i].intensity;//scale color of operator with intensity of pseudo symmetry
				const double order = emsphinx::Constants<double>::pi / (std::acos(maxima[i].qu.w));//get order of rotation axis
				const double nFld = std::round(order);//what is the closest integer order
				if(std::fabs(nFld - nFld) < rotCut) {
					//this is a crystallographic pseudosymmetry operator, use the right symbol
					diag.setColor(c.rgb[0], c.rgb[1], c.rgb[2]);
					diag.addRotor(maxima[i].qu.x, maxima[i].qu.y, maxima[i].qu.z, (int)nFld, 0.67);
				} else {
					//just use a dot for other rotations
					diag.addDot(maxima[i].qu.data()+1, c, 0.67);
				}
			}
		}
	}

	//save updated diagram
	diag.getHemi().write("pseudo.svg");
	return 0;
	//build a weighted histogram of misorientation angles

	//rescale xc to [0,1]
	/*
	auto minMax = std::minmax_element(s2Corr.xc.cbegin(), s2Corr.xc.cend());
	const double vMin = *(minMax.first);
	const double range = *(minMax.second) - *(minMax.first);
	for(double& v : s2Corr.xc) {
		v = (v - vMin) / range;
	}

	*/

	// //build volume elements
	// // const double dv0 = xtal::Constants<double>::pi2 * xtal::Constants<double>::pi2 * xtal::Constants<double>::pi2 / (sl * sl * sl); 
	// const double dv0 = xtal::Constants<double>::pi / (sl * sl * sl);// normalized by 8 pi^2
	// std::vector<double> dv(bw+1);
	// for(size_t k = 0; k <= bw; k++) {
	// 	double eu[3];
	// 	s2Corr.indexEuler(k * sl * sl, eu);
	// 	dv[k] = dv0 * std::fabs(std::sin(eu[1]));//sin(beta)
	// }

	//normalize cross correlation
	/*
	double mean = 0, den = 0;
	for(size_t k = 0; k <= bw; k++) {
		std::vector<double>::const_iterator iter = s2Corr.xc.cbegin() + k * sl * sl;
		mean += std::accumulate(iter, iter + sl * sl, 0.0) * dv[k];
		den += dv[k] * sl * sl;
	}
	mean /= den;
	for(double& v : s2Corr.xc) v -= mean;
	double stdev = 0;
	for(size_t k = 0; k <= bw; k++) {
		std::vector<double>::const_iterator iter = s2Corr.xc.cbegin() + k * sl * sl;
		stdev += std::inner_product(iter, iter + sl * sl, iter, 0.0) * dv[k];
	}
	stdev = std::sqrt(stdev / den);
	for(double& v : s2Corr.xc) v /= stdev;
	*/
	/*
	const double mean = std::accumulate(s2Corr.xc.cbegin(), s2Corr.xc.cend(), 0.0) / s2Corr.xc.size();
	for(double& v : s2Corr.xc) v -= mean;
	const double stdev = std::sqrt(std::inner_product(s2Corr.xc.cbegin(), s2Corr.xc.cend(), s2Corr.xc.cbegin(), 0.0) / s2Corr.xc.size());
	for(double& v : s2Corr.xc) v /= stdev;


	std::ofstream ofs("xcorr.raw", std::ios::out | std::ios::binary);
	ofs.write((char*)s2Corr.xc.data(), s2Corr.xc.size() * sizeof(double));

	//construct misorientation angle bins
	const size_t numBins = sl;
	double maxDiso = 180;
	std::vector<double> bins(numBins);//bin limits
	std::vector<double> cnts(numBins, 0);//counts
	std::vector<double> wgts(numBins, 0);//weights
	for(size_t i = 0; i < numBins; i++) bins[i] = maxDiso * (i+1) / numBins;

	//loop over euler grid
	xtal::Quat<double> qi = xtal::Quat<double>::Identity();
	for(size_t k = 0; k <= bw; k++) {
		for(size_t m = 0; m < sl; m++) {
			for(size_t n = 0; n < sl; n++) {
				//get euler angle at this grid point
				const size_t idx = k * sl * sl + m * sl + n;
				double eu[3];
				s2Corr.indexEuler(idx, eu);

				//convert to quaternion and get misorientation
				xtal::Quat<double> qu;
				xtal::zyz2qu(eu, qu.data());
				p.pointGroup().disoQu(qu.data(), qi.data(), qu.data());
				
				//find misorientation bin
				const double w = std::acos(qu.w) * 2 * xtal::Constants<double>::rd2dg;//misorientatino angle
				const size_t i = std::distance(bins.cbegin(), std::upper_bound(bins.cbegin(), bins.cend(), w));

				//compute volume element and accumulate
				if(i < bins.size()) {
					wgts[i] += dv[k];
					cnts[i] += dv[k] * std::exp(s2Corr.xc[idx]);
					// cnts[i] += dv[k] * s2Corr.xc[idx];
				}
			}
		}
	}
	const double sumWgt = std::accumulate(cnts.cbegin(), cnts.cend(), 0.0);

	//write histogram to file
	std::ofstream csv("hist.csv");
	csv << "w,cnt,wgt\n";
	for(size_t i = 0; i < bins.size(); i++) {
		csv << bins[i] << ',' << cnts[i]/sumWgt << ',' << wgts[i] << '\n';
	}
	*/
	return 0;
}
