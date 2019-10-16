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

#ifndef _EBSD_IDX_H_
#define _EBSD_IDX_H_

#define MINIZ_NO_STDIO
#define MINIZ_NO_TIME
#define MINIZ_NO_ZLIB_APIS
#include "miniz/miniz.c"

#include <fstream>

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

#include <atomic>
#include <vector>
#include <functional>

#include "util/threadpool.hpp"
#include "util/image.hpp"
#include "xtal/orientation_map.hpp"
#include "modality/ebsd/pattern.hpp"
#include "modality/ebsd/nml.hpp"
#include "modality/ebsd/imprc.hpp"
#include "modality/ebsd/detector.hpp"
#include "idx/indexer.hpp"

namespace emsphinx {

	namespace ebsd {
		
		template <typename Real>
		struct IndexingData {
		private:
			//@brief     : initialize indexing
			//@param nml : namelist to initailze with
			//@param pIpf: location to write ipf colorsnamelist to initailze with
			void initialize(Namelist& nml, char * pIpf);

		public:

		    IndexingData(IndexingData const&) = delete;
		    IndexingData& operator=(IndexingData const&) = delete;

		    IndexingData(Namelist& nml, char * pIpf = NULL) {initialize(nml, pIpf);}

			std::vector< MasterSpectra<Real> >                        phases     ;//master patterns to index against
			std::shared_ptr<PatternFile>                              pat        ;//raw patterns to index
			std::vector< xtal::OrientationMap<Real> >                 om         ;//orientation maps to index into (1 per phase + 1 combined for multi phase)
			std::vector<char>                                         idxMask    ;//mask of point to index
			Geometry<Real>                                            geom       ;//detector geometry
			std::vector< std::unique_ptr< emsphinx::Indexer<Real> > > indexers   ;//actual indexers
			std::vector< std::vector< char > >                        patBufs    ;//work space for patterns
			std::function<void(size_t)>                               workItem   ;//indexing work item
			Namelist                                                  namelist   ;

			size_t                                                    numIdx     ;//number of point to index
			std::atomic<uint64_t>                                     idxCtr     ;//indexed pattern counter
			size_t                                                    threadCount;//number of worker threads
			std::vector< char >                                       ipfBuf     ;//ipf map

			//@brief  : save results 
			//@param s: start time
			//@param e: end time
			//@param t: total time in s
			void save(time_t& s, time_t& e, double t);
		};

		//@brief    : templated work item for an ebsd indexing run
		//@param om : orientation map(s) to index into [1 map for each result to save]
		//@param pat: patterns to index
		//@param idx: indexers
		//@param buf: location to extract patterns for each thread
		//@param cnt: batch size in # of patterns
		//@param msk: indexing bitmask (0x01 to index, 0x02 to refine, 0x03 for both)
		//@param ctr: completed pattern counter
		//@param ipf: location to write ipf map
		template <typename TPix, typename Real>
		std::function<void(size_t)> ebsdWorkItem(
			std::vector<xtal::OrientationMap<Real> >                 & om ,//map is shared (each pixel is accessed by only 1 thread)
			std::shared_ptr<PatternFile>                               pat,//patterns are shared (thread safe)
			std::vector< std::unique_ptr< emsphinx::Indexer<Real> > >& idx,//1 indexer per thread
			std::vector<                  std::vector< char >  >     & buf,//1 buffer per thread
			const size_t                                               cnt,//1 buffer per thread
			char const * const                                         msk,//indexing bitmask
			std::atomic<uint64_t>&                                     ctr,//counter for indexed pattern count
			char       * const                                         ipf //ipf map
			);


		//@brief: helper to compute image quality in parallel
		template <typename Real>
		class ThreadedIqCalc {
			ThreadPool                                                      pool ;//worker threads
			std::vector< std::shared_ptr< image::ImageQualityCalc<Real> > > calcs;//image quality calculators
			std::vector<                  std::vector<char>               > buffs;//pattern read buffers
			std::vector<                  std::vector<Real>               > dBufs;//pattern conversion buffers
			std::atomic<int64_t>                                            cmplt;//count of completed patterns

			public:
				//@brief     : construct a threaded IQ calculator
				//@param nThd: number of worker threads
				ThreadedIqCalc(const size_t nThd = std::max<size_t>(ThreadPool::Concurrency() / 2, 1)) : pool(nThd), calcs(nThd), buffs(nThd), dBufs(nThd) {cmplt.store(-1);}

				//@brief: destructor cancels work immediatley (you must keep the object in scope long enough for completion)
				~ThreadedIqCalc() {pool.clear();}

				//@brief    : set the pattern file and start computing (non-blocking)
				//@param pat: patterns to compute for
				//@param iq : vector to start writing image quality
				//@param upd: callback function for updates (called once when calculations have started and once when all calculations are queued)
				void start(std::shared_ptr<PatternFile> pat, std::vector<Real>& iq, std::function<void()> upd = NULL);

				//@brief        : wait for the thread pool to finish
				//@param timeout: how long to wait
				//@return       : number of completed patterns (-1 for unstarted)
				template <class Rep, class Period = std::ratio<1> >
				int64_t wait(std::chrono::duration<Rep, Period> timeout) const {pool.waitAll(timeout); return cmplt.load();}
		};

		//@brief     : initialize indexing
		//@param nml : namelist to initailze with
		//@param pIpf: location to write ipf colorsnamelist to initailze with
		template <typename Real>
		void IndexingData<Real>::initialize(Namelist& nml, char * pIpf) {
			namelist = nml;
			nml.writeFileHeader();//save meta data up front

			//read master patterns
			for(const std::string& file : nml.masterFiles) {
				phases.push_back(MasterSpectra<Real>(file));//compute SHT of master patterns once
				phases.back().resize(nml.bw);
			}
			for(size_t i = 1; i < phases.size(); i++) {
				if(phases.front().getSig() != phases[i].getSig() || phases.front().getKv() != phases[i].getKv()) throw std::runtime_error("all master patterns must have the same tilt and kV");
			}

			//add psuedo-symmetry
			size_t numMatches = 1;
			if(!nml.pSymFile.empty()) {
				if(1 != phases.size()) throw std::runtime_error("psuedo-symmetry files currently only supported for single phase indexing");
				phases[0].addPseudoSym(nml.pSymFile);
			}

			//read patterns and sanity check number / size
			pat = PatternFile::Read(nml.patFile, nml.patName, nml.patDims[0], nml.patDims[1]);
			if(pat->numPat() < nml.scanDims[0] * nml.scanDims[1]) throw std::runtime_error("not enough patterns for scan");
			else if(pat->numPat() > nml.scanDims[0] * nml.scanDims[1]) {
				// std::cout << "\n * * * * * * warning: " << pat->numPat() - nml.scanDims[0] * nml.scanDims[1] << " more pattern(s) than scan points * * * * * * \n\n";
			}
			if(pat->width() != nml.patDims[0] || pat->height() != nml.patDims[1]) throw std::runtime_error("detector/pattern size mismatch");

			//build mask and convert from 0/1 to 0/3 if needed (for old style indexing bitmask)
			idxMask = nml.roi.buildMask(nml.scanDims[0], nml.scanDims[1]);
			numIdx = std::count_if(idxMask.cbegin(), idxMask.cend(), [](const uint8_t& v){return v != 0;});
			if(0 == numIdx) {
				std::fill(idxMask.begin(), idxMask.end(), 1);
				numIdx = idxMask.size();
			}
			if(nml.refine) {
				for(char& v : idxMask) {
					if(v == 1) v = 3;
				}
			}

			//read tilt from master pattern monte carlo data and build geometry
			geom.sampleTilt(phases.front().getSig());
			geom.cameraTilt(nml.thetac);
			geom.cameraSize(nml.patDims[0], nml.patDims[1], nml.delta);
			if("EMsoft" == nml.ven) {
				geom.patternCenter      (nml.pctr[0], nml.pctr[1], nml.pctr[2]);
			} else if("EDAX" == nml.ven || "tsl" == nml.ven) {
				geom.patternCenterTSL   (nml.pctr[0], nml.pctr[1], nml.pctr[2]);
			} else if("Oxford" == nml.ven) {
				geom.patternCenterOxford(nml.pctr[0], nml.pctr[1], nml.pctr[2]);
			} else if("Bruker" == nml.ven || "tsl" == nml.ven) {
				geom.patternCenterBruker(nml.pctr[0], nml.pctr[1], nml.pctr[2]);
			} else throw std::runtime_error("unknown vendor for pattern center `" + nml.ven + "'");
			geom.maskPattern(nml.circRad == 0);
			geom.flipPattern(pat->flipY());

			//determine threading / batch values
			threadCount = nml.nThread == 0 ? ThreadPool::Concurrency() : nml.nThread;
			if(0 == nml.batchSize) nml.batchSize = (int32_t)Indexer<Real>::BatchEstimate(nml.bw, threadCount, numIdx);

			//create orientation map from input scan dimensions and add phases from master patterns
			xtal::OrientationMap<Real> newOm(nml.scanDims[0], nml.scanDims[1], nml.scanSteps[0], nml.scanSteps[1]);
			for(const MasterSpectra<Real>& mp : phases) newOm.phsList.push_back(mp.phase()); 
			newOm.calib.sTlt = phases.front().getSig();
			const size_t extraScans = 1 == phases.size() ? phases.front().pseudoSym().size() : 1;//do we need extra scans for psuedosymmetry or multi phase
			om.assign(phases.size() + extraScans, newOm);

			//build image processor
			std::unique_ptr< PatternProcessor<Real> > prc(new PatternProcessor<Real>());
			prc->setSize(geom.w, geom.h, nml.circRad, nml.gausBckg, nml.nRegions);

			//build back projector
			const size_t gridDim = nml.bw + (nml.bw % 2 == 0 ? 3 : 2);
			std::array<Real, 4> quNp = geom.northPoleQuat();
			std::unique_ptr< BackProjector<Real> > prj(new BackProjector<Real>(geom, gridDim, std::sqrt(Real(2)), quNp.data()));

			//build cross correlators
			std::vector< std::unique_ptr< sphere::PhaseCorrelator<Real> > > corrs;
			if(nml.normed) {
				//build a spherical harmonic transformer once
				square::DiscreteSHT<Real> sht(gridDim, nml.bw, square::Layout::Legendre);
				std::vector<Real> sph(gridDim * gridDim * 2, 0);//spherical grid

				//compute the SHT of the window function once
				std::vector<Real> win(pat->numPix(), Real(1));//image of solid 1
				prj->unproject(win.data(), sph.data());//unproject window function onto sphere
				std::vector< std::complex<Real> > mlm(nml.bw * nml.bw);
				sht.analyze(sph.data(), mlm.data());

				//now loop over the phases computing SHT(function^2) and building normalized correlators
				std::vector< std::complex<Real> > flm2(nml.bw * nml.bw);
				for(const MasterSpectra<Real>& p : phases) {
					sht.synthesize(p.data(), sph.data());//reconstruct function on legendre grid
					std::transform(sph.cbegin(), sph.cend(), sph.cbegin(), sph.begin(), std::multiplies<Real>());//square function

					sht.analyze(sph.data(), flm2.data());//compute SHT(function^2)
					std::unique_ptr< sphere::NormalizedCorrelator<Real> > pCorr(new sphere::NormalizedCorrelator<Real>((int)nml.bw, p.data(), flm2.data(), p.mirror(), p.nFold(), mlm.data() ) );//build correlator
					corrs.push_back( std::move( pCorr ) );
				}
			} else {
				//build unnormalized spherical cross correlators
				for(const MasterSpectra<Real>& p : phases) {
					std::shared_ptr< std::vector< std::complex<Real> > > flm = std::make_shared< std::vector< std::complex<Real> > >(p.data(), p.data() + nml.bw * nml.bw);//copy harmonics into shared pointer
					std::unique_ptr< sphere::UnNormalizedCorrelator<Real> > pCorr(new sphere::UnNormalizedCorrelator<Real>((int)nml.bw, flm, p.mirror(), p.nFold() ) );//build correlator
					corrs.push_back( std::move( pCorr ) );
				}
			}

			//build a single indexer for each thread
			std::unique_ptr<emsphinx::Indexer<Real> > idx(new emsphinx::Indexer<Real>((int)nml.bw, prc->clone(), prj->clone(), corrs));//make a single indexer (good for 1 thread)
			indexers.push_back(std::move(idx));//move original indexer onto 
			for(size_t i = 1; i < threadCount; i++) indexers.push_back( std::move(indexers.front()->clone()) );//duplicate n-1 times

			//build thread work item
			idxCtr.store(0);
			patBufs.assign(threadCount, std::vector<char>(nml.batchSize * pat->imBytes()));
			if(NULL == pIpf) {
				ipfBuf.resize(nml.scanDims[0] * nml.scanDims[1] * 3);
				pIpf = ipfBuf.data();
			}
			switch(pat->pixelType()) {//templated lambda functions don't exist yet so we'll need to select the appropriate type here
				case PatternFile::Bits::U8 : workItem = ebsdWorkItem<uint8_t , Real>(om, pat, indexers, patBufs, nml.batchSize, idxMask.data(), std::ref(idxCtr), pIpf); break;
				case PatternFile::Bits::U16: workItem = ebsdWorkItem<uint16_t, Real>(om, pat, indexers, patBufs, nml.batchSize, idxMask.data(), std::ref(idxCtr), pIpf); break;
				case PatternFile::Bits::F32: workItem = ebsdWorkItem<float   , Real>(om, pat, indexers, patBufs, nml.batchSize, idxMask.data(), std::ref(idxCtr), pIpf); break;
				case PatternFile::Bits::UNK: throw std::runtime_error("unsupported pattern pixel type");
			}
		}

		//@brief  : save results 
		//@param s: start time
		//@param e: end time
		//@param t: total time in s
		template <typename Real>
		void IndexingData<Real>::save(time_t& s, time_t& e, double t) {
			//reopen HDF output
			H5::H5File file(namelist.opath + namelist.dataFile, H5F_ACC_RDWR);//we made the file earlier
			{//write indexing start / stop times
				std::stringstream ss;
				ss << std::put_time(std::localtime(&s), "%a %b %d %Y %T");
				std::string startStr = ss.str();
				ss.str("");
				ss << std::put_time(std::localtime(&e), "%a %b %d %Y %T");
				std::string endStr = ss.str();
				ss.str("");
				double rate = double(numIdx) / t;
				file.openGroup("EMheader").createDataSet("StartTime", H5::StrType(0, H5T_VARIABLE), H5::DataSpace(H5S_SCALAR)).write(startStr, H5::StrType(0, H5T_VARIABLE));
				file.openGroup("EMheader").createDataSet("StopTime" , H5::StrType(0, H5T_VARIABLE), H5::DataSpace(H5S_SCALAR)).write(endStr  , H5::StrType(0, H5T_VARIABLE));
				file.openGroup("EMheader").createDataSet("PatPerS"  , H5::PredType::NATIVE_DOUBLE , H5::DataSpace(H5S_SCALAR)).write(&rate   , H5::PredType::NATIVE_DOUBLE );
			}

			//choose IPF reference direction and hsl function
			Real n[3] = {0, 0, 1};
			std::function<void(Real const*const, Real*const)> h2r = xtal::sph2rgb<Real>;

			//save out vendor file and images if needed
			std::vector<uint8_t> ipfMap = om.front().ipfColor(n, h2r);
			std::vector<uint8_t> xcMap  = image::to8Bit(om.front().metric.data(), om.front().metric.size());
			std::vector<uint8_t> iqMap  = image::to8Bit(om.front().imQual.data(), om.front().imQual.size());
			if(!namelist.vendorFile.empty()) om.front().write(namelist.opath + namelist.vendorFile);
			if(!namelist.ipfName   .empty()) writePng(ipfMap.data(), om.front().width, om.front().height, 3, namelist.opath + namelist.ipfName );
			if(!namelist.qualName  .empty()) writePng(xcMap .data(), om.front().width, om.front().height, 1, namelist.opath + namelist.qualName);

			//write write all scans to hdf
			size_t padNum = (size_t) ceil(std::log10(om.size()+1));//how many digits does the longest number take
			for(size_t i = 0; i < om.size(); i++) {
				//make maps
				ipfMap = om.front().ipfColor(n, h2r);
				xcMap  = image::to8Bit(om.front().metric.data(), om.front().metric.size());
				iqMap  = image::to8Bit(om.front().imQual.data(), om.front().imQual.size());

				//build string prefix "0i_"
				std::stringstream ss;
				ss << std::setfill('0') << std::setw(padNum) << i + 1;

				//write scan data
				om[i].name = "Match " + ss.str();
				hsize_t dims[3] = {om[i].height, om[i].width, 3};
				H5::Group grp = file.createGroup("Scan " + ss.str());
				om[i].writeH5(grp);
				// H5IMmake_image_8bit(grp.getId(), "IPF Map", om[i].width, om[i].height, ipfMap.data());//this requires HL liberary
				grp.createDataSet("IPF Map", H5::PredType::NATIVE_UINT8, H5::DataSpace(3, dims)).write(ipfMap.data(), H5::PredType::NATIVE_UINT8);
				grp.createDataSet("XC Map" , H5::PredType::NATIVE_UINT8, H5::DataSpace(2, dims)).write(xcMap .data(), H5::PredType::NATIVE_UINT8);
				grp.createDataSet("IQ Map" , H5::PredType::NATIVE_UINT8, H5::DataSpace(2, dims)).write(iqMap .data(), H5::PredType::NATIVE_UINT8);
			}

		}

		//@brief    : templated work item for an ebsd indexing run
		//@param om : orientation map(s) to index into [1 map for each result to save]
		//@param pat: patterns to index
		//@param idx: indexers
		//@param buf: location to extract patterns for each thread
		//@param cnt: batch size in # of patterns
		//@param msk: indexing bitmask (0x01 to index, 0x02 to refine, 0x03 for both)
		//@param ctr: completed pattern counter
		//@param ipf: location to write ipf map
		template <typename TPix, typename Real>
		std::function<void(size_t)> ebsdWorkItem(
			std::vector<xtal::OrientationMap<Real> >                 & om ,//map is shared (each pixel is accessed by only 1 thread)
			std::shared_ptr<PatternFile>                               pat,//patterns are shared (thread safe)
			std::vector< std::unique_ptr< emsphinx::Indexer<Real> > >& idx,//1 indexer per thread
			std::vector<                  std::vector< char >  >     & buf,//1 buffer per thread
			const size_t                                               cnt,//1 buffer per thread
			char const * const                                         msk,//indexing bitmask
			std::atomic<uint64_t>&                                     ctr,//counter for indexed pattern count
			char       * const                                         ipf //ipf map
			)  {
			//@brief      : index a batch of pattern
			//@param count: maximum number of patterns to index
			//@param id   : id of current thread
			return [&om,pat,&idx,&buf,cnt,msk,&ctr,ipf](const size_t id){//work function
				//choose IPF reference direction and hsl function
				Real n[3] = {0, 0, 1}, nx[3], rgb[3];
				std::function<void(Real const*const, Real*const)> h2r = xtal::sph2rgb<Real>;

				//extract patterns and get pointer to first pattern
				const std::vector<size_t> indices = pat->extract(buf[id].data(), cnt);
				const size_t bytes = pat->imBytes();
				char const * ptr = buf[id].data();

				//loop over patterns indexing
				std::vector< emsphinx::Result<Real> > res(om.size());//we need a result for each map
				for(const size_t i : indices) {//loop over indices to index (may not be contiguous or in order)
					if(0x00 != msk[i]) {//do we need to index and/or refine this point?
						if(0x01 & msk[i]) {//reindexing requested
							const bool ref = 0x02 & msk[i];//should we refine?
							try {
								idx[id]->indexImage((TPix*)ptr, res.data(), res.size(), ref);//index pattern
								for(size_t j = 0; j < res.size(); j++) {
									std::copy(res[j].qu, res[j].qu+4, om[j].qu[i].data());//save orientation
									om[j].metric[i] =               res[j].corr ;//save cross correlation
									om[j].imQual[i] =               res[j].iq   ;//image quality
									om[j].phase [i] = (uint_fast8_t)res[j].phase;//save phase
								}
								std::fill(ipf + 3 * i, ipf + 3 * i + 3, 0);
								xtal::quat::rotateVector(res[0].qu, n, nx);
								if (-1 == res[0].phase) {
									std::fill(rgb, rgb + 3, Real(0));
								} else {
									om[0].phsList[res[0].phase].pg.ipfColor(nx, rgb, h2r);
								}
								std::transform(rgb, rgb+3, ipf + 3 * i, [](const Real& v){return (char) std::round(v * 255);});
							} catch (std::exception&) {// e) {
								// std::cout << i << ": " << e.what() << '\n';//don't let one exception stop everything
								for(size_t j = 0; j < res.size(); j++) {
									om[j].qu[i].w = 1;
									om[j].qu[i].x = om[j].qu[i].y = om[j].qu[i].z = 0;
									om[j].metric[i] = 0;
									om[j].imQual[i] = 0;
									om[j].phase [i] = (uint_fast8_t)-1;
								}
								std::fill(ipf + 3 * i, ipf + 3 * i + 3, 0);
							}
						} else if(0x02 & msk[i]) {//only refinement requested
							try {
								//this currently only uses a single result but it could be modified
								res[0].phase = om[0].phase[i];
								std::copy(om[0].qu[i].data(), om[0].qu[i].data()+4, res[0].qu);//copy original orientation
								idx[id]->refineImage((TPix*)ptr, res[0]);//refine orientation using pattern
								std::copy(res[0].qu, res[0].qu+4, om[0].qu[i].data());//save orientation
								om[0].metric[i] = res[0].corr ;//save cross correlation
								om[0].imQual[i] = res[0].iq   ;//image quality
							} catch (std::exception&) {// e) {
								// std::cout << i << ": " << e.what() << '\n';//don't let one exception stop everything
							}
						}
						ctr++;//increment indexed pattern count
					}
					ptr += bytes;//move to next pattern
				}
			};
		}

		//@brief    : set the pattern file and start computing (non-blocking)
		//@param pat: pattern file to set
		//@param iq : vector to start writing image quality
		//@param upd: callback function for updates (called once when calculations have started and once when all calculations are queued)
		template <typename Real>
		void ThreadedIqCalc<Real>::start(std::shared_ptr<PatternFile> pat, std::vector<Real>& iq, std::function<void()> upd) {
			//build calculators
			const size_t batSz = 10;//arbitrarily set batch size for now
			calcs.front() = std::make_shared< image::ImageQualityCalc<Real> >(pat->width(), pat->height());//build a single image calculator
			for(size_t i = 1; i < calcs.size(); i++) calcs[i] = std::make_shared< image::ImageQualityCalc<Real> >(*calcs.front());//copy the rest
			for(size_t i = 0; i < calcs.size(); i++) {//allocate buffers
				buffs[i].resize(batSz * pat->imBytes ());
				dBufs[i].resize(        pat->numPix  ());
			}

			//compute number of batches
			const size_t totNm = pat->numPat();
			size_t batches = totNm / batSz;//how many batches are needed
			if(batches * batSz < totNm) ++batches;//extra batch for leftovers

			//build a work item to process a single batch
			std::function<void(size_t)> work = [batSz,pat,&iq,this](const size_t id) {
				//extract a batch of patterns and get pointer to first pattern
				const std::vector<size_t> indices = pat->extract(this->buffs[id].data(), batSz);
				const size_t bytes = pat->imBytes();
				char const * ptr = this->buffs[id].data();

				//loop over patterns computing image quality
				for(const size_t i : indices) {//loop over indices to index (may not be contiguous or in order)
					switch(pat->pixelType()) {//convert pattern to double
						case PatternFile::Bits::UNK: std::fill(this->dBufs[id].begin(), this->dBufs[id].end(), 0); break;
						case PatternFile::Bits::U8 : std::copy((uint8_t  const*)ptr, (uint8_t  const*)(ptr+bytes), this->dBufs[id].begin()); break;
						case PatternFile::Bits::U16: std::copy((uint16_t const*)ptr, (uint16_t const*)(ptr+bytes), this->dBufs[id].begin()); break;
						case PatternFile::Bits::F32: std::copy((float    const*)ptr, (float    const*)(ptr+bytes), this->dBufs[id].begin()); break;
					}
					iq[i] = this->calcs[id]->compute(this->dBufs[id].data());//compute image quality
					(this->cmplt)++;//increment count
					ptr += bytes;//move to next pattern
				}
			};

			//queue work
			cmplt.store(0);//we haven't finished any patterns
			iq.resize(totNm);//make sure there is space to write outputs
			if(NULL != upd) upd();
			for(size_t i = 0; i < batches; i++) {//loop over batches
				const size_t start = i * batSz;//first pattern
				const size_t end = std::min(start + batSz, pat->numPat());//last pattern
				pool.schedule(std::bind(work, std::placeholders::_1));//queue computatino
			}
			if(NULL != upd) upd();
		}

	}//ebsd

}//emsphinx



#endif//_EBSD_IDX_H_

