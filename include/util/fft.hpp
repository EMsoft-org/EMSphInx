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

#ifndef _FFTW_WRAP_
#define _FFTW_WRAP_

#include <complex>
#include <type_traits>
#include <vector>

#include "fftw3.h"

//preprocessor macros to use fftw float, double, and/or long
//these can be hard coded here if needed but should be define by cmake
// #define EM_USE_F
// #define EM_USE_D
// #define EM_USE_L

namespace fft {
	namespace flag {
		enum class Plan {
			Estimate   = FFTW_ESTIMATE   ,//only estimate fastest algorithm (fastest plan creation but slowest execution)
			Measure    = FFTW_MEASURE    ,//measure a few paths to determine fastest algorithm [~seconds to create]
			Patient    = FFTW_PATIENT    ,//measure many paths to determine fastest algorithm [several times slower than measure]
			Exhaustive = FFTW_EXHAUSTIVE ,//measure all paths to determine fastest algorithm (slowest plan creation but fastest execution) [significantly longer]
			// WisdomOnly = FFTW_WISDOM_ONLY,//only create a plan if wisdom is available
		};
		
		enum class Input {
			Destroy   = FFTW_DESTROY_INPUT ,//out of place transforms can destroy input
			Preserve  = FFTW_PRESERVE_INPUT,//out of place transform must preserve input
			// Unaligned = FFTW_UNALIGNED     ,//prevent special alignment restrictions (disables SIMD acceleration)
		};
	}

	namespace detail {
		//templated helper to hold an fftw plan
		template <typename Real> struct Plan {static_assert(std::is_same<Real, float>::value || std::is_same<Real, double>::value || std::is_same<Real, long double>::value, "Real must be float, double");};

		//helper to manage fftw wisdom
		template <typename Real>
		struct Wisdom {
			//@brief: import existing wisdom from a file
			void read();

			//@brief: export accumulated wisdom to a file
			void write();

			Wisdom() {read();}//import existing wisdom on creation
			~Wisdom() {write();}//export existing wisdom on cleanup

			//@brief    : custom signal handler for file reading
			//@param sig: signal number (SIGILL)
			//@note     : converts to exception in case the wisdom file contains SIMD -> illegal instruction
			static void IllegalInstructionHandler(int sig);
		};

		//single global instance for types of interest to automatically load/save wisdom
	#ifdef EM_USE_F
		Wisdom<     float > fWisdom;
	#endif
	#ifdef EM_USE_D
		Wisdom<     double> dWisdom;
	#endif
	#ifdef EM_USE_L
		Wisdom<long double> lWisdom;
	#endif
	}

	//@brief: find the closest FFT size that is fast
	//@param x: minimum FFT
	//@return : smallest fast FFT size that is greater or equal to x
	uint32_t fastSize(const uint32_t x);

	//helper to allow stl containers with fftw_malloc
	template <class T> struct allocator {
		typedef T value_type;//typedef for allocator_traits
		                          allocator(                   ) noexcept {                                       }//constructor for allocator_traits
		template <class U>        allocator(const allocator<U>&) noexcept {                                       }//template constructor for allocator_traits
		                   T*     allocate (      size_t n     )          {return (T*) fftw_malloc(n * sizeof(T));}//allocate with fftw_malloc instead of new
		                   void deallocate (T* p, size_t n     )          {            fftw_free  ((void*) p    );}//free with fftw_free instead of delete
	};
	template<class T> using vector = std::vector<T, allocator<T> >;//use fft::vector<T> instead of std::vector<T> to use fftw allocation

	//templated helper to wrap a pair of forward/reverse transformations for real data
	template <typename Real> struct RealFFT {
		detail::Plan<Real> pFor, pInv;

		//@brief        : construct the forward and reverse FFT plans
		//@param n      : length of FFT
		//@param pFlag  : fft planning flag
		//@param destroy: can the input be destroyed during execution
		RealFFT(const size_t n, const flag::Plan pFlag, const bool destroy = false);

		//@brief        : compute an FFT from input data
		//@param signal : data to compute FFT of
		//@param spectra: location to write FFT of signal
		void forward(Real* signal, std::complex<Real>* spectra) const;

		//@brief        : compute an inverse FFT from input FFT
		//@param spectra: spectra to reconstruct signal from
		//@param signal : location to write reconstructed data
		void inverse(std::complex<Real>* spectra, Real* signal) const;
	};

	//templated helper to wrap an inverse transformation of 3D real data
	template <typename Real> struct RealFFT3D {
		detail::Plan<Real> pInv;

		//@brief      : construct the forward and reverse FFT plans
		//@param n    : side length of 3D FFT
		//@param pFlag: fft planning flag
		RealFFT3D(const size_t n, const flag::Plan pFlag);//mutli dimensional c2r transforms cannot preserve input

		//@brief        : compute an inverse FFT from input FFT
		//@param spectra: spectra to reconstruct signal from
		//@param signal : location to write reconstructed data
		void inverse(std::complex<Real>* spectra, Real* signal) const;
	};

	//templated helper to wrap a inverse transformation of 3D real data (seperated into components for each dimension)
	template <typename Real> struct SepRealFFT3D {
		const int vN, vH;//full and halfcomplex sizes
		detail::Plan<Real> pX, pY, pZ;//planes for transforms along x, y, and z axis

		//@brief      : construct the FFT plans
		//@param n    : side length of 3D FFT
		//@param pFlag: fft planning flag
		SepRealFFT3D(const size_t n, const flag::Plan pFlag);//mutli dimensional c2r transforms cannot preserve input

		//@brief        : compute an inverse FFT from input FFT
		//@param spectra: spectra to reconstruct signal from
		//@param signal : location to write reconstructed data
		//@param dx     : frequency of YZ planes with non zero elements
		void inverse(std::complex<Real>* spectra, Real* signal, const size_t dx = 1) const;
	};

	template <typename Real> struct DCT2D {
		detail::Plan<Real> plan;

		//@brief        : construct DCT plan
		//@param width  : 2d array width (fast index)
		//@param height : 2d array height (slow index)
		//@param pFlag  : fft planning flag
		//@param reverse: true/false for reverse/forward cosine transform
		//@param destroy: can the input be destroyed during execution
		DCT2D(const size_t width, const size_t height, const flag::Plan pFlag, const bool reverse, const bool destroy = false);

		//@brief       : execute the planned cosine transform on a new array
		//@param input : input signal/spectra
		//@param output: location to write output signal/spectra
		void execute(Real* input, Real* output) const;
	};

	template <typename Real> struct DCT {
		detail::Plan<Real> plan;

		//@brief        : construct DCT plan
		//@param length : 1d array size
		//@param pFlag  : fft planning flag
		//@param reverse: true/false for reverse/forward cosine transform
		//@param destroy: can the input be destroyed during execution
		DCT(const size_t length, const flag::Plan pFlag, const bool reverse, const bool destroy = false);

		//@brief       : execute the planned cosine transform on a new array
		//@param input : input signal/spectra
		//@param output: location to write output signal/spectra
		void execute(Real* input, Real* output) const;
	};
}

////////////////////////////////////////////////////////////////////////
//                          Implementations                           //
////////////////////////////////////////////////////////////////////////

#include <vector>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <mutex>
#include <set>
#include <algorithm>
#include <csignal>

#include "sysnames.hpp"

namespace fft {
	namespace detail {
		//float
	#ifdef EM_USE_F
		template <> struct Plan <float> {
			fftwf_plan p  ;//fftw plane
			std::mutex mut;//mutex for cleanup when plan is shared across threads

			Plan() : p(NULL) {}//make sure we don't try to clean up random memory on destruction

			//destructor (needs to be thread safe)
			~Plan() {
				std::unique_lock<std::mutex> lock(mut);//make sure another thread isn't currently calling the destructor
				if(NULL != p) fftwf_destroy_plan(p);//don't destroy plan if it was already cleaned up of never allocoated
			}
		};
	#endif

		//double
	#ifdef EM_USE_D
		template <> struct Plan <double> {
			fftw_plan  p  ;//fftw plane 
			std::mutex mut;//mutex for cleanup when plan is shared across threads

			Plan() : p(NULL) {}//make sure we don't try to clean up random memory on destruction

			//destructor (needs to be thread safe)
			~Plan() {
				std::unique_lock<std::mutex> lock(mut);//make sure another thread isn't currently calling the destructor
				if(NULL != p) fftw_destroy_plan(p);}//don't destroy plan if it was already cleaned up of never allocoated

		};
	#endif

		//long double
	#ifdef EM_USE_L
		template <> struct Plan <long double> {
			fftwl_plan p  ;//fftw plane
			std::mutex mut;//mutex for cleanup when plan is shared across threads

			Plan() : p(NULL) {}//make sure we don't try to clean up random memory on destruction

			//destructor (needs to be thread safe)
			~Plan() {
				std::unique_lock<std::mutex> lock(mut);//make sure another thread isn't currently calling the destructor
				if(NULL != p) fftwl_destroy_plan(p);//don't destroy plan if it was already cleaned up of never allocoated
			}
		};
	#endif

			//@brief    : custom signal handler for file reading
			//@param sig: signal number (SIGILL)
			//@note     : converts to exception in case the wisdom file contains SIMD -> illegal instruction
	#ifdef EM_USE_F
		template <> void Wisdom<     float >::IllegalInstructionHandler(int sig) {
			if(SIGILL == sig) {
				const std::string fileName = getSharedDataDir() + "fftwf.wisdom";
				throw std::runtime_error("illegal instruction while reading wisdom from " + fileName + ", try deleting wisdom file");
			}
			exit(sig);
		}
	#endif
	#ifdef EM_USE_D
		template <> void Wisdom<     double>::IllegalInstructionHandler(int sig) {
			if(SIGILL == sig) {
				const std::string fileName = getSharedDataDir() + "fftw.wisdom";
				throw std::runtime_error("illegal instruction while reading wisdom from " + fileName + ", try deleting wisdom file");
			}
			exit(sig);
		}
	#endif
	#ifdef EM_USE_L
		template <> void Wisdom<long double>::IllegalInstructionHandler(int sig) {
			if(SIGILL == sig) {
				const std::string fileName = getSharedDataDir() + "fftwl.wisdom";
				throw std::runtime_error("illegal instruction while reading wisdom from " + fileName + ", try deleting wisdom file");
			}
			exit(sig);
		}
	#endif

		//@brief: helper function to check if a file exists
		//@param name: name of file to look for
		//@return: true/false if file does/doesn't exist
		bool fileExists(std::string name) {
			std::ifstream is(name);
			return is.good();
		}

		//@brief: import existing wisdom from a file
	#ifdef EM_USE_F
		template <> void Wisdom<     float >::read() {
			const std::string fileName = getSharedDataDir() + "fftwf.wisdom";
			if(fileExists(fileName)) { //only try if the file exists
				signal(SIGILL, &Wisdom<     float >::IllegalInstructionHandler);//switch to custom illegal instruction handler
				const bool imported = fftwf_import_wisdom_from_filename(fileName.c_str());//try to read wisdom
				signal(SIGILL, SIG_DFL);//switch back to default illegal instruction handler
				if(!imported) throw std::runtime_error("failed to read wisdom from " + fileName + ", try deleting wisdom file");
			}
		}
	#endif
	#ifdef EM_USE_D
		template <> void Wisdom<     double>::read() {
			const std::string fileName = getSharedDataDir() + "fftw.wisdom";
			if(fileExists(fileName)) { //only try if the file exists
				signal(SIGILL, &Wisdom<     double>::IllegalInstructionHandler);//switch to custom illegal instruction handler
				const bool imported = fftw_import_wisdom_from_filename (fileName.c_str());//try to read wisdom
				signal(SIGILL, SIG_DFL);//switch back to default illegal instruction handler
				if(!imported) throw std::runtime_error("failed to read wisdom from " + fileName + ", try deleting wisdom file");
			}
		}
	#endif
	#ifdef EM_USE_L
		template <> void Wisdom<long double>::read() {
			const std::string fileName = getSharedDataDir() + "fftwl.wisdom";
			if(fileExists(fileName)) { //only try if the file exists
				signal(SIGILL, &Wisdom<long double>::IllegalInstructionHandler);//switch to custom illegal instruction handler
				const bool imported = fftwl_import_wisdom_from_filename(fileName.c_str());//try to read wisdom
				signal(SIGILL, SIG_DFL);//switch back to default illegal instruction handler
				if(!imported) throw std::runtime_error("failed to read wisdom from " + fileName + ", try deleting wisdom file");
			}
		}
	#endif

		//@brief: export accumulated wisdom to a file
	#ifdef EM_USE_F
		template <> void Wisdom<     float >::write() {
			const std::string fileName = getSharedDataDir() + "fftwf.wisdom";
			if(!fftwf_export_wisdom_to_filename(fileName.c_str()))
				std::cerr << "failed to write wisdom to " << fileName << '\n';//destructor can't throw, at least print a warning
		}
	#endif
	#ifdef EM_USE_D
		template <> void Wisdom<     double>::write() {
			const std::string fileName = getSharedDataDir() + "fftw.wisdom";
			if(!fftw_export_wisdom_to_filename (fileName.c_str()))
				std::cerr << "failed to write wisdom to " << fileName << '\n';//destructor can't throw, at least print a warning
		}
	#endif
	#ifdef EM_USE_L
		template <> void Wisdom<long double>::write() {
			const std::string fileName = getSharedDataDir() + "fftwl.wisdom";
			if(!fftwl_export_wisdom_to_filename(fileName.c_str()))
				std::cerr << "failed to write wisdom to " << fileName << '\n';//destructor can't throw, at least print a warning
		}
	#endif
	}

	//@brief: find the closest FFT size that is fast
	//@param x: minimum FFT
	//@return : smallest fast FFT size that is greater or equal to x
	uint32_t fastSize(const uint32_t x) {
		//handle special/easy cases
		if(x <= 16) return std::max<uint32_t>(1, x);//fftw explicitly implements all ffts from 1 -> 16

		//start by computing the next power of 2 up from x
		//https://graphics.stanford.edu/~seander/bithacks.html#RoundUpPowerOf2
		uint32_t v2x = x;//start with first power of 2 >= x
		v2x--;
		v2x |= v2x >>  1;
		v2x |= v2x >>  2;
		v2x |= v2x >>  4;
		v2x |= v2x >>  8;
		v2x |= v2x >> 16;
		v2x++;//v2x now holds first poewr of 2 >= x

		//now compute the log_2 of v2x
		//https://graphics.stanford.edu/~seander/bithacks.html#IntegerLog
		static const uint32_t b[] = {0xAAAAAAAA, 0xCCCCCCCC, 0xF0F0F0F0,  0xFF00FF00, 0xFFFF0000};
		uint32_t ln2x = (v2x & b[0]) != 0;
		for(int i = 4; i > 0; i--) ln2x |= ((v2x & b[i]) != 0) << i;

		//next compute log_2(log_2(v2x)) [since we'll be squaring in the last step]
		uint32_t maxIter = (ln2x & b[0]) != 0;
		for(int i = 4; i > 0; i--) maxIter |= ((ln2x & b[i]) != 0) << i;

		//now compute all combinations of 2^i * 3^j * 5^k... for fast (small) primes
		//small primes for fftw are: 2, 3, 5, 7, 11, and 13
		//i, j, k ... only need to be checked for i < r;
		std::set<uint32_t> s;
		s.insert( 2);
		s.insert( 3);
		s.insert( 5);
		s.insert( 7);
		s.insert(11);
		s.insert(13);
		uint32_t vMin = v2x;//our initial guess for smallest fast size is the next power of 2 up
		for(uint32_t iter = 0; iter < maxIter; iter++) {//loop over required iterations
			std::set<uint32_t> sNew;//set to hold new elements to be added
			for(const uint32_t& i : s) {//loop over current elements once
				for(const uint32_t& j : s) {//loop over current elements twice
					const uint32_t v = i * j;//compute product of elements
					if(v < vMin) {//is this element small enough to care about (smaller than our current best)
						if(v < x) {//values less than x are prefactors
							sNew.insert(v);
						} else {//otherwise (>= x) we've found a new best size
							vMin = v;
						}
					}
				}
			}
			s.insert(sNew.cbegin(), sNew.cend());//add our new pre factors
		}
		return vMin;
	}

	//@brief        : construct the forward and reverse FFT plans
	//@param n      : length of FFT
	//@param pFlag  : fft planning flag
	//@param destroy: can the input be destroyed during execution
#ifdef EM_USE_F
	template <> RealFFT<     float >::RealFFT(const size_t n, const flag::Plan pFlag, const bool destroy) {
		std::vector<     float > work1(n);
		std::vector< std::complex<     float > > work2(n);
		pFor.p = fftwf_plan_dft_r2c_1d((int)n, work1.data(), (fftwf_complex*)work2.data(), (destroy ? FFTW_DESTROY_INPUT : FFTW_PRESERVE_INPUT) | (int)pFlag);
		pInv.p = fftwf_plan_dft_c2r_1d((int)n, (fftwf_complex*)work2.data(), work1.data(), (destroy ? FFTW_DESTROY_INPUT : FFTW_PRESERVE_INPUT) | (int)pFlag);
	}
#endif
#ifdef EM_USE_D
	template <> RealFFT<     double>::RealFFT(const size_t n, const flag::Plan pFlag, const bool destroy) {
		std::vector<     double> work1(n);
		std::vector< std::complex<     double> > work2(n);
		pFor.p = fftw_plan_dft_r2c_1d ((int)n, work1.data(), (fftw_complex *)work2.data(), (destroy ? FFTW_DESTROY_INPUT : FFTW_PRESERVE_INPUT) | (int)pFlag);
		pInv.p = fftw_plan_dft_c2r_1d ((int)n, (fftw_complex *)work2.data(), work1.data(), (destroy ? FFTW_DESTROY_INPUT : FFTW_PRESERVE_INPUT) | (int)pFlag);
	}
#endif
#ifdef EM_USE_L
	template <> RealFFT<long double>::RealFFT(const size_t n, const flag::Plan pFlag, const bool destroy) {
		std::vector<long double> work1(n);
		std::vector< std::complex<long double> > work2(n);
		pFor.p = fftwl_plan_dft_r2c_1d((int)n, work1.data(), (fftwl_complex*)work2.data(), (destroy ? FFTW_DESTROY_INPUT : FFTW_PRESERVE_INPUT) | (int)pFlag);
		pInv.p = fftwl_plan_dft_c2r_1d((int)n, (fftwl_complex*)work2.data(), work1.data(), (destroy ? FFTW_DESTROY_INPUT : FFTW_PRESERVE_INPUT) | (int)pFlag);
	}
#endif

	//@brief        : compute an FFT from input data
	//@param signal : data to compute FFT of
	//@param spectra: location to write FFT of signal
#ifdef EM_USE_F
	template <> void RealFFT<     float >::forward(     float * signal, std::complex<     float >* spectra) const {fftwf_execute_dft_r2c(pFor.p, signal, (fftwf_complex*)spectra);}
#endif
#ifdef EM_USE_D
	template <> void RealFFT<     double>::forward(     double* signal, std::complex<     double>* spectra) const {fftw_execute_dft_r2c (pFor.p, signal, (fftw_complex *)spectra);}
#endif
#ifdef EM_USE_L
	template <> void RealFFT<long double>::forward(long double* signal, std::complex<long double>* spectra) const {fftwl_execute_dft_r2c(pFor.p, signal, (fftwl_complex*)spectra);}
#endif

	//@brief        : compute an inverse FFT from input FFT
	//@param spectra: spectra to reconstruct signal from
	//@param signal : location to write reconstructed data
#ifdef EM_USE_F
	template <> void RealFFT<     float >::inverse(std::complex<     float >* spectra,      float * signal) const {fftwf_execute_dft_c2r(pInv.p, (fftwf_complex*)spectra, signal);}
#endif
#ifdef EM_USE_D
	template <> void RealFFT<     double>::inverse(std::complex<     double>* spectra,      double* signal) const {fftw_execute_dft_c2r (pInv.p, (fftw_complex *)spectra, signal);}
#endif
#ifdef EM_USE_L
	template <> void RealFFT<long double>::inverse(std::complex<long double>* spectra, long double* signal) const {fftwl_execute_dft_c2r(pInv.p, (fftwl_complex*)spectra, signal);}
#endif

	//@brief      : construct the forward and reverse FFT plans
	//@param n    : side length of 3D FFT
	//@param pFlag: fft planning flag
#ifdef EM_USE_F
	template <> RealFFT3D<     float >::RealFFT3D(const size_t n, const flag::Plan pFlag) {
		std::vector<     float > work1(n*n*n);
		std::vector< std::complex<     float > > work2(n*n*n);
		pInv.p = fftwf_plan_dft_c2r_3d((int)n, (int)n, (int)n, (fftwf_complex*)work2.data(), work1.data(), (int)pFlag);
	}
#endif
#ifdef EM_USE_D
	template <> RealFFT3D<     double>::RealFFT3D(const size_t n, const flag::Plan pFlag) {
		std::vector<     double> work1(n*n*n);
		std::vector< std::complex<     double> > work2(n*n*n);
		pInv.p = fftw_plan_dft_c2r_3d ((int)n, (int)n, (int)n, (fftw_complex *)work2.data(), work1.data(), (int)pFlag);
	}
#endif
#ifdef EM_USE_L
	template <> RealFFT3D<long double>::RealFFT3D(const size_t n, const flag::Plan pFlag) {
		std::vector<long double> work1(n*n*n);
		std::vector< std::complex<long double> > work2(n*n*n);
		pInv.p = fftwl_plan_dft_c2r_3d((int)n, (int)n, (int)n, (fftwl_complex*)work2.data(), work1.data(), (int)pFlag);
	}
#endif

	//@brief        : compute an inverse FFT from input FFT
	//@param spectra: spectra to reconstruct signal from
	//@param signal : location to write reconstructed data
#ifdef EM_USE_F
	template <> void RealFFT3D<     float >::inverse(std::complex<     float >* spectra,      float * signal) const {fftwf_execute_dft_c2r(pInv.p, (fftwf_complex*)spectra, signal);}
#endif
#ifdef EM_USE_D
	template <> void RealFFT3D<     double>::inverse(std::complex<     double>* spectra,      double* signal) const {fftw_execute_dft_c2r (pInv.p, (fftw_complex *)spectra, signal);}
#endif
#ifdef EM_USE_L
	template <> void RealFFT3D<long double>::inverse(std::complex<long double>* spectra, long double* signal) const {fftwl_execute_dft_c2r(pInv.p, (fftwl_complex*)spectra, signal);}
#endif

	//@brief      : construct the FFT plans
	//@param n    : side length of 3D FFT
	//@param pFlag: fft planning flag
#ifdef EM_USE_F
	template <> SepRealFFT3D<     float >::SepRealFFT3D(const size_t n, const flag::Plan pFlag) : vN(int(n)), vH(int(n)/2 + 1) {
		std::vector<                   float   > work1(vN*vN*vN);
		std::vector< std::complex<     float > > work2(vN*vN*vH);
		int rank           = 1                                   ;//individual transforms are all 1D
		int nn[1]          = {vN}                                ;//individual transforms are all of length n
		int howmany        = vN                                  ;//how many transformations will be performed (one down each z for each y at a single x)
		fftwf_complex* in  = (fftwf_complex*) work2.data()       ;//input data
		int* inembed       = NULL                                ;//dimensions of super array that input  is row major subarray of (null for not a subarray)
		int istride        = vN * vH                             ;//stride between sequential elements (z spacing)
		int idist          = vH                                  ;//kth fft input at in + k * idist (y spacing)
		fftwf_complex* out = (fftwf_complex*) work1.data()       ;//output data
		int* oenembed      = NULL                                ;//dimensions of super array that output is row major subarray of (null for not a subarray)
		int ostride        = 1                                   ;//output stride
		int odist          = vN                                  ;//kth fft outputs to out + k * odist
		int sign           = FFTW_BACKWARD                       ;//inverse transform
		unsigned flags     = FFTW_DESTROY_INPUT | (unsigned)pFlag;//planning flags
		pZ.p = fftwf_plan_many_dft    (rank, nn, howmany, in , inembed, istride, idist, out         , oenembed, ostride, odist  , sign, flags);//1st: transform down z for all y at a single x (into output array)
		pX.p = fftwf_plan_many_dft_c2r(rank, nn, howmany, in , inembed, 1      , idist, work1.data(), oenembed, ostride, odist  ,       flags);//3rd: transform down x for all y at a single z (into output array)
		pY.p = fftwf_plan_many_dft    (rank, nn, vH     , out, inembed, vN     , 1    , in          , oenembed, vH     , vN * vH, sign, flags);//2nd: transform down y for all z at a single x (into original 3d input from output)
	}
#endif
#ifdef EM_USE_D
	template <> SepRealFFT3D<     double>::SepRealFFT3D(const size_t n, const flag::Plan pFlag) : vN(int(n)), vH(int(n)/2 + 1) {
		std::vector<                   double  > work1(vN*vN*vN);
		std::vector< std::complex<     double> > work2(vN*vN*vH);
		int rank           = 1                                   ;//individual transforms are all 1D
		int nn[1]          = {vN}                                ;//individual transforms are all of length n
		int howmany        = vN                                  ;//how many transformations will be performed (one down each z for each y at a single x)
		fftw_complex * in  = (fftw_complex *) work2.data()       ;//input data
		int* inembed       = NULL                                ;//dimensions of super array that input  is row major subarray of (null for not a subarray)
		int istride        = vN * vH                             ;//stride between sequential elements (z spacing)
		int idist          = vH                                  ;//kth fft input at in + k * idist (y spacing)
		fftw_complex * out = (fftw_complex *) work1.data()       ;//output data
		int* oenembed      = NULL                                ;//dimensions of super array that output is row major subarray of (null for not a subarray)
		int ostride        = 1                                   ;//output stride
		int odist          = vN                                  ;//kth fft outputs to out + k * odist
		int sign           = FFTW_BACKWARD                       ;//inverse transform
		unsigned flags     = FFTW_DESTROY_INPUT | (unsigned)pFlag;//planning flags
		pZ.p = fftw_plan_many_dft     (rank, nn, howmany, in , inembed, istride, idist, out         , oenembed, ostride, odist  , sign, flags);//1st: transform down z for all y at a single x (into output array)
		pX.p = fftw_plan_many_dft_c2r (rank, nn, howmany, in , inembed, 1      , idist, work1.data(), oenembed, ostride, odist  ,       flags);//3rd: transform down x for all y at a single z (into output array)
		pY.p = fftw_plan_many_dft     (rank, nn, vH     , out, inembed, vN     , 1    , in          , oenembed, vH     , vN * vH, sign, flags);//2nd: transform down y for all z at a single x (into original 3d input from output)
	}
#endif
#ifdef EM_USE_L
	template <> SepRealFFT3D<long double>::SepRealFFT3D(const size_t n, const flag::Plan pFlag) : vN(int(n)), vH(int(n)/2 + 1) {
		std::vector<              long double  > work1(vN*vN*vN);
		std::vector< std::complex<long double> > work2(vN*vN*vH);
		int rank           = 1                                   ;//individual transforms are all 1D
		int nn[1]          = {vN}                                ;//individual transforms are all of length n
		int howmany        = vN                                  ;//how many transformations will be performed (one down each z for each y at a single x)
		fftwl_complex* in  = (fftwl_complex*) work2.data()       ;//input data
		int* inembed       = NULL                                ;//dimensions of super array that input  is row major subarray of (null for not a subarray)
		int istride        = vN * vH                             ;//stride between sequential elements (z spacing)
		int idist          = vH                                  ;//kth fft input at in + k * idist (y spacing)
		fftwl_complex* out = (fftwl_complex*) work1.data()       ;//output data
		int* oenembed      = NULL                                ;//dimensions of super array that output is row major subarray of (null for not a subarray)
		int ostride        = 1                                   ;//output stride
		int odist          = vN                                  ;//kth fft outputs to out + k * odist
		int sign           = FFTW_BACKWARD                       ;//inverse transform
		unsigned flags     = FFTW_DESTROY_INPUT | (unsigned)pFlag;//planning flags
		pZ.p = fftwl_plan_many_dft    (rank, nn, howmany, in , inembed, istride, idist, out         , oenembed, ostride, odist  , sign, flags);//1st: transform down z for all y at a single x (into output array)
		pX.p = fftwl_plan_many_dft_c2r(rank, nn, howmany, in , inembed, 1      , idist, work1.data(), oenembed, ostride, odist  ,       flags);//3rd: transform down x for all y at a single z (into output array)
		pY.p = fftwl_plan_many_dft    (rank, nn, vH     , out, inembed, vN     , 1    , in          , oenembed, vH     , vN * vH, sign, flags);//2nd: transform down y for all z at a single x (into original 3d input from output)
	}
#endif

	//@brief        : compute an inverse FFT from input FFT
	//@param spectra: spectra to reconstruct signal from
	//@param signal : location to write reconstructed data
	//@param dx     : frequency of YZ planes with non zero elements
#ifdef EM_USE_F
	template <> void SepRealFFT3D<     float >::inverse(std::complex<     float >* spectra,      float * signal, const size_t dx) const {
		const int delta = (int)delta;
		for(int i = 0; i < vH; i+= delta) {//loop over yz planes doing 2d transforms
			fftwf_execute_dft(pZ.p, (fftwf_complex*)spectra + i, (fftwf_complex*)signal     );
			fftwf_execute_dft(pY.p, (fftwf_complex*)signal     , (fftwf_complex*)spectra + i);
		}
		for(int i = 0; i < vH; i++) fftwf_execute_dft_c2r(pX.p, (fftwf_complex*)spectra + vN * vH * i, signal + vN * vN * i);//loop up xy planes doing batches of 1d c2r transforms, does 1 extra plane for extraction of 3x3x3 neighborhood at upper glide
	}
#endif
#ifdef EM_USE_D
	template <> void SepRealFFT3D<     double>::inverse(std::complex<     double>* spectra,      double* signal, const size_t dx) const {
		const int delta = (int)dx;
		for(int i = 0; i < vH; i+= delta) {//loop over yz planes doing 2d transforms
			fftw_execute_dft (pZ.p, (fftw_complex *)spectra + i, (fftw_complex *)signal     );
			fftw_execute_dft (pY.p, (fftw_complex *)signal     , (fftw_complex *)spectra + i);
		}
		for(int i = 0; i < vH; i++) fftw_execute_dft_c2r (pX.p, (fftw_complex *)spectra + vN * vH * i, signal + vN * vN * i);//loop up xy planes doing batches of 1d c2r transforms, does 1 extra plane for extraction of 3x3x3 neighborhood at upper glide
	}
#endif
#ifdef EM_USE_L
	template <> void SepRealFFT3D<long double>::inverse(std::complex<long double>* spectra, long double* signal, const size_t dx) const {
		const int delta = (int)dx;
		for(int i = 0; i < vH; i+= delta) {//loop over yz planes doing 2d transforms
			fftwf_execute_dft(pZ.p, (fftwl_complex*)spectra + i, (fftwl_complex*)signal     );
			fftwf_execute_dft(pY.p, (fftwl_complex*)signal     , (fftwl_complex*)spectra + i);
		}
		for(int i = 0; i < vH; i++) fftwl_execute_dft_c2r(pX.p, (fftwl_complex*)spectra + vN * vH * i, signal + vN * vN * i);//loop up xy planes doing batches of 1d c2r transforms, does 1 extra plane for extraction of 3x3x3 neighborhood at upper glide
	}
#endif

	//@brief        : construct DCT plan
	//@param width  : 2d array width (fast index)
	//@param height : 2d array height (slow index)
	//@param pFlag  : fft planning flag
	//@param reverse: true/false for reverse/forward cosine transform
	//@param destroy: can the input be destroyed during execution
#ifdef EM_USE_F
	template <> DCT2D<     float >::DCT2D(const size_t width, const size_t height, const flag::Plan pFlag, const bool reverse, const bool destroy) {
		std::vector<     float > work1(width*height), work2(width*height);
		plan.p = fftwf_plan_r2r_2d((int)height , (int)width , work1.data(), work2.data(), reverse ? FFTW_REDFT01 : FFTW_REDFT10, reverse ? FFTW_REDFT01 : FFTW_REDFT10, (destroy ? FFTW_DESTROY_INPUT : FFTW_PRESERVE_INPUT) | (int)pFlag);
	}
#endif
#ifdef EM_USE_D
	template <> DCT2D<     double>::DCT2D(const size_t width, const size_t height, const flag::Plan pFlag, const bool reverse, const bool destroy) {
		std::vector<     double> work1(width*height), work2(width*height);
		plan.p = fftw_plan_r2r_2d ((int)height , (int)width , work1.data(), work2.data(), reverse ? FFTW_REDFT01 : FFTW_REDFT10, reverse ? FFTW_REDFT01 : FFTW_REDFT10, (destroy ? FFTW_DESTROY_INPUT : FFTW_PRESERVE_INPUT) | (int)pFlag);
	}
#endif
#ifdef EM_USE_L
	template <> DCT2D<long double>::DCT2D(const size_t width, const size_t height, const flag::Plan pFlag, const bool reverse, const bool destroy) {
		std::vector<long double> work1(width*height), work2(width*height);
		plan.p = fftwl_plan_r2r_2d((int)height , (int)width , work1.data(), work2.data(), reverse ? FFTW_REDFT01 : FFTW_REDFT10, reverse ? FFTW_REDFT01 : FFTW_REDFT10, (destroy ? FFTW_DESTROY_INPUT : FFTW_PRESERVE_INPUT) | (int)pFlag);
	}
#endif

	//@brief       : execute the planned cosine transform on a new array
	//@param input : input signal/spectra
	//@param output: location to write output signal/spectra
#ifdef EM_USE_F
	template <> void DCT2D<     float >::execute(     float * input,      float * output) const {fftwf_execute_r2r(plan.p, input, output);}
#endif
#ifdef EM_USE_D
	template <> void DCT2D<     double>::execute(     double* input,      double* output) const {fftw_execute_r2r (plan.p, input, output);}
#endif
#ifdef EM_USE_L
	template <> void DCT2D<long double>::execute(long double* input, long double* output) const {fftwl_execute_r2r(plan.p, input, output);}
#endif

	//@brief        : construct DCT plan
	//@param length : 1d array size
	//@param pFlag  : fft planning flag
	//@param reverse: true/false for reverse/forward cosine transform
	//@param destroy: can the input be destroyed during execution
#ifdef EM_USE_F
	template <> DCT<     float >::DCT(const size_t length, const flag::Plan pFlag, const bool reverse, const bool destroy) {
		std::vector<     float > work1(length), work2(length);
		plan.p = fftwf_plan_r2r_1d((int)length, work1.data(), work2.data(), reverse ? FFTW_REDFT01 : FFTW_REDFT10, (destroy ? FFTW_DESTROY_INPUT : FFTW_PRESERVE_INPUT) | (int)pFlag);
	}
#endif
#ifdef EM_USE_D
	template <> DCT<     double>::DCT(const size_t length, const flag::Plan pFlag, const bool reverse, const bool destroy) {
		std::vector<     double> work1(length), work2(length);
		plan.p = fftw_plan_r2r_1d ((int)length, work1.data(), work2.data(), reverse ? FFTW_REDFT01 : FFTW_REDFT10, (destroy ? FFTW_DESTROY_INPUT : FFTW_PRESERVE_INPUT) | (int)pFlag);
	}
#endif
#ifdef EM_USE_L
	template <> DCT<long double>::DCT(const size_t length, const flag::Plan pFlag, const bool reverse, const bool destroy) {
		std::vector<long double> work1(length), work2(length);
		plan.p = fftwl_plan_r2r_1d((int)length, work1.data(), work2.data(), reverse ? FFTW_REDFT01 : FFTW_REDFT10, (destroy ? FFTW_DESTROY_INPUT : FFTW_PRESERVE_INPUT) | (int)pFlag);
	}
#endif

	//@brief       : execute the planned cosine transform on a new array
	//@param input : input signal/spectra
	//@param output: location to write output signal/spectra
#ifdef EM_USE_F
	template <> void DCT<     float >::execute(     float * input,      float * output) const {fftwf_execute_r2r(plan.p, input, output);}
#endif
#ifdef EM_USE_D
	template <> void DCT<     double>::execute(     double* input,      double* output) const {fftw_execute_r2r (plan.p, input, output);}
#endif
#ifdef EM_USE_L
	template <> void DCT<long double>::execute(long double* input, long double* output) const {fftwl_execute_r2r(plan.p, input, output);}
#endif
}

#endif//_FFTW_WRAP_
