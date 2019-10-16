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

#ifndef _BMP_H_
#define _BMP_H_

#include <vector>
#include <iostream>

namespace bmp {
	//bitmap info header
	struct Header {
		//@brief   : read info from an istream
		//@param is: istream to read from
		//@return  : is
		std::istream& read(std::istream& is);

		//@brief    : read info from a raw pointer
		//@param ptr: raw buffer to read from
		//@return   : true/false if read was successful/failed
		void read(char const * const ptr);

		//@brief    : determine the number of bytes needed to store the image
		//@param gry: true to get the number of bytes to write as a grayscale image instead
		uint32_t bytes(const bool gry = false) const;

		//@brief    : read the image into a raw buffer
		//@param is : istream to read from
		//@param buf: location to write image
		//@param gry: true to convert to grayscale (red channel for color, cast for 1/4 bit)
		//@return   : is
		std::istream& readImage(std::istream& is, char * const buf, const bool gry = false) const;

		//@brief    : read the image into a raw buffer
		//@param ptr: raw buffer to read from
		//@param buf: location to write image
		//@param gry: true to convert to grayscale (red channel for color, cast for 1/4 bit)
		void readImage(char const * const ptr, char * const buf, const bool gry = false) const;

		//Header
		char     signature[2];//magic bytes
		uint32_t fileSize    ;//size of the file in bytes
		uint16_t res1, res2  ;//reserved space
		uint32_t offset      ;//offset to image data

		//Info Header

		// modified BITMAPCOREHEADER (16 bytes) (actual BITMAPCOREHEADER is 12 bytes)
		uint32_t  size       ;// size of this structure in bytes
		uint32_t  width      ;// bitmap width in pixels  // only int16 in bitmapcoreheader only
		uint32_t  height     ;// bitmap height in pixels // only int16 in bitmapcoreheader only
		uint16_t  planes     ;// must be 1
		uint16_t  bitCount   ;// bits per pixel (1, 4, 8, or 24)

		// BITMAPINFOHEADER (40 bytes) adds these fields and extends bitCount to include (0, 16, and 32)
		uint32_t  compression;// type of compression (0:none,1:8bit rle,2:4bit rle,3:indexed,4:jpg,5:png)
		uint32_t  sizeImage  ;// size of image in bytes
		uint32_t  xRes       ;// x resolution in pixels per meter
		uint32_t  yRes       ;// y resolution in pixels per meter
		uint32_t  colorsUsed ;// number of lut values used
		uint32_t  colorsImprt;// number of lut values needed to display the image

		// BITMAPV4HEADER (108 bytes) adds these fields
		uint32_t  redMask    ;// red   bits in rgb image
		uint32_t  greenMask  ;// green bits in rgb image
		uint32_t  blueMask   ;// blue  bits in rgb image
		uint32_t  alphaMask  ;// alpha bits in rgba image
		uint32_t  colorSpace ;// color space (flag for if cie endpoints are given)
		uint32_t  redX       ;// x coordinate of red in cie
		uint32_t  redY       ;// y '                      '
		uint32_t  redZ       ;// z '                      '
		uint32_t  greenX     ;// '
		uint32_t  greenY     ;//                green
		uint32_t  greenZ     ;//                          '
		uint32_t  blueX      ;// '
		uint32_t  blueY      ;//                 blue
		uint32_t  blueZ      ;//                          '
		uint32_t  gammaRed   ;// red   gamma curve value
		uint32_t  gammaGreen ;// green gamma curve value
		uint32_t  gammaBlue  ;// blue  gamma curve value

		// BITMAPV5HEADER (124 bytes) adds these fields and support for additional colorSpace types
		uint32_t  intent     ;// rendering intent
		uint32_t  profileData;// offset in bytes from beginning of header to start of profile data
		uint32_t  profileSize;// size in bytes of profile data
		uint32_t  reserved   ;// should be 0

		private:
			//@brief   : read BITMAPCOREHEADER from an istream
			//@param is: istream to read from
			//@return  : is
			std::istream& readCore(std::istream& is);

			//@brief   : read BITMAPINFOHEADER from an istream
			//@param is: istream to read from
			//@return  : is
			std::istream& readInfo(std::istream& is);

			//@brief   : read BITMAPV4HEADER from an istream
			//@param is: istream to read from
			//@return  : is
			std::istream& readV4(std::istream& is);

			//@brief   : read BITMAPV5HEADER from an istream
			//@param is: istream to read from
			//@return  : is
			std::istream& readV5(std::istream& is);

			//@brief: byteswap to little endian from big if needed
			void headerToLittle();

			//@brief: byteswap to little endian from big if needed
			void infoToLittle();
	};
}

struct Bitmap {
	//@brief         : read a bitmap from a file
	//@param fileName: file to read from
	//@param gry     : true to get the number of bytes to write as a grayscale image instead
	void read(std::string fileName);

	bmp::Header       header;
	std::vector<char> buff  ;//data in row major order (each row padded to nearest byte)
};

////////////////////////////////////////////////////////////////////////
//                          Implementations                           //
////////////////////////////////////////////////////////////////////////

#include <fstream>

namespace bmp {
	namespace detail {
		//@brief : check if the system is big or little endian
		//@return: true/false for big/little endian system
		bool bigEndian() {
			static const union {
				uint16_t i   ;
				char     c[2];
			} u = {0x0001};
			static const bool big = 0x00 == u.c[0];
			return big;
		}

		//@brief    : byteswap fixed size integers
		//@param uxx: xx bit unsigned integer to byteswap
		//@return   : byteswapped integer
		uint16_t swapU16(const uint16_t& u16) {return ((u16<< 8)&0xFF00) | ((u16>> 8)&0x00FF);}
		uint32_t swapU32(const uint32_t& u32) {return ((u32<<24)&0xFF000000) | ((u32<< 8)&0x00FF0000) | ((u32>> 8)&0x0000FF00) | ((u32>>24)&0x000000FF);}
	} 

	////////////////////////////////////////////////////////////////////////
	//                               Header                               //
	////////////////////////////////////////////////////////////////////////

	//@brief   : read info from an istream
	//@param is: istream to read from
	//@return  : is
	std::istream& Header::read(std::istream& is) {
		//read header
		is.read(signature       , 2);
		is.read((char*)&fileSize, 4);
		is.read((char*)&res1    , 2);
		is.read((char*)&res2    , 2);
		is.read((char*)&offset  , 4);
		headerToLittle();
		if('B' != signature[0] || 'M' != signature[1]) throw std::runtime_error("not a valid bitmap file");

		//read size
		is.read((char*)&size, 4);
		if(detail::bigEndian()) size = detail::swapU32(size);

		//select reader based on header size
		switch(size) {
			case  12: readCore(is); break;//BITMAPCOREHEADER
			case  40: readInfo(is); break;//BITMAPINFOHEADER
			case 108: readV4  (is); break;//BITMAPV4HEADER
			case 124: readV5  (is); break;//BITMAPV5HEADER
			default : throw std::runtime_error("unsupported header size");
		}
		infoToLittle();
		return is;
	}

	//@brief    : read info from a raw pointer
	//@param ptr: raw buffer to read from
	//@return   : true/false if read was successful/failed
	void Header::read(char const * const ptr) {
		//read header
		signature[0] = ptr[0];
		signature[1] = ptr[1];
		std::copy(ptr +  2, ptr +  6, (char*)&fileSize);
		std::copy(ptr +  6, ptr +  8, (char*)&res1    );
		std::copy(ptr +  8, ptr + 10, (char*)&res2    );
		std::copy(ptr + 10, ptr + 14, (char*)&offset  );
		headerToLittle();
		if(!('B' == signature[0] && 'M' == signature[1])) throw std::runtime_error("not a valid bitmap file");

		//read size
		char const * const pInfo = ptr + 14;
		std::copy(pInfo, pInfo + 4, (char*)&size);
		if(detail::bigEndian()) size = detail::swapU32(size);

		//read fields based on header size
		switch(size) {
			case 124://BITMAPV5HEADER
				// BITMAPV5HEADER (124 bytes) adds these fields and support for additional colorSpace types
				std::copy(pInfo + 108, pInfo + 112, (char*)&intent     );
				std::copy(pInfo + 112, pInfo + 116, (char*)&profileData);
				std::copy(pInfo + 116, pInfo + 120, (char*)&profileSize);
				std::copy(pInfo + 120, pInfo + 124, (char*)&reserved   );
				//intentional fall through

			case 108://BITMAPV4HEADER
				// BITMAPV4HEADER (108 bytes) adds these fields
				std::copy(pInfo +  40, pInfo +  44, (char*)&redMask    );
				std::copy(pInfo +  44, pInfo +  48, (char*)&greenMask  );
				std::copy(pInfo +  48, pInfo +  52, (char*)&blueMask   );
				std::copy(pInfo +  52, pInfo +  56, (char*)&alphaMask  );
				std::copy(pInfo +  56, pInfo +  60, (char*)&colorSpace );
				std::copy(pInfo +  60, pInfo +  64, (char*)&redX       );
				std::copy(pInfo +  64, pInfo +  68, (char*)&redY       );
				std::copy(pInfo +  68, pInfo +  72, (char*)&redZ       );
				std::copy(pInfo +  72, pInfo +  76, (char*)&greenX     );
				std::copy(pInfo +  76, pInfo +  80, (char*)&greenY     );
				std::copy(pInfo +  80, pInfo +  84, (char*)&greenZ     );
				std::copy(pInfo +  84, pInfo +  88, (char*)&blueX      );
				std::copy(pInfo +  88, pInfo +  92, (char*)&blueY      );
				std::copy(pInfo +  92, pInfo +  96, (char*)&blueZ      );
				std::copy(pInfo +  96, pInfo + 100, (char*)&gammaRed   );
				std::copy(pInfo + 100, pInfo + 104, (char*)&gammaGreen );
				std::copy(pInfo + 104, pInfo + 108, (char*)&gammaBlue  );
				//intentional fall through

			case  40://BITMAPINFOHEADER
				// BITMAPINFOHEADER (40 bytes) adds these fields and extends bitCount to include (0, 16, and 32)
				std::copy(pInfo +  16, pInfo +  20, (char*)&compression);
				std::copy(pInfo +  20, pInfo +  24, (char*)&sizeImage  );
				std::copy(pInfo +  24, pInfo +  28, (char*)&xRes       );
				std::copy(pInfo +  28, pInfo +  32, (char*)&yRes       );
				std::copy(pInfo +  32, pInfo +  36, (char*)&colorsUsed );
				std::copy(pInfo +  36, pInfo +  40, (char*)&colorsImprt);

				// modified BITMAPCOREHEADER (16 bytes) (actual BITMAPCOREHEADER is 12 bytes)
				std::copy(pInfo +   4, pInfo +   8, (char*)&width      );
				std::copy(pInfo +   8, pInfo +  12, (char*)&height     );
				std::copy(pInfo +  12, pInfo +  14, (char*)&planes     );
				std::copy(pInfo +  14, pInfo +  16, (char*)&bitCount   );
				break;

			case  12: {//BITMAPCOREHEADER
				uint16_t tmp;
				std::copy(pInfo +   4, pInfo +   6, (char*)&tmp        ); width  = tmp;
				std::copy(pInfo +   6, pInfo +   8, (char*)&tmp        ); height = tmp;
				std::copy(pInfo +   8, pInfo +  10, (char*)&planes     );
				std::copy(pInfo +  10, pInfo +  12, (char*)&bitCount   );
			} break;

			default: throw std::runtime_error("unsupported header size");
		}
		infoToLittle();
	}
	
	//@brief    : determine the number of bytes needed to store the image
	//@param gry: true to get the number of bytes to write as a grayscale image instead
	uint32_t Header::bytes(const bool gry) const {
		if(gry) {
			switch(bitCount) {
				case  0: return 0;

				case 16: return width * height * 2;//16 bit

				case  1: //bitmask -> 8bit
				case  4: //4 bit -> 8 bit
				case  8: //already 8 bit
				case 24: //rgb -> 8 bit
				case 32: return width * height;//rgba -> 8 bit

				default: throw std::runtime_error("unsupported bits per pixel");
			}
		} else {
			switch(bitCount) {
				case  0: return 0;
				case  1: return (width * height + 7)/8;
				case  4: return (width * height + 1)/2;
				case  8: return width * height    ;
				case 16: return width * height * 2;
				case 24: return width * height * 3;
				case 32: return width * height * 4;
				default: throw std::runtime_error("unsupported bits per pixel");
			}
		}
	}

	//@brief    : read the image into a raw buffer
	//@param is : istream to read from
	//@param buf: location to write image
	//@param gry: true to convert to grayscale (red channel for color, cast for 1/4 bit)
	//@return   : is
	std::istream& Header::readImage(std::istream& is, char * const buf, const bool gry) const {
		//make sure we can read this bitmap
		if(0 != compression) throw std::runtime_error("unsupported compression");
		//compute row padding
		uint32_t rowBytes = width;
		switch(bitCount) {
			case  0: rowBytes = 0;               return is;
			case  1: throw std::runtime_error("unsupported bit depth for reading");//these are going to be annoying to implement properly and probably aren't needed
			case  4: throw std::runtime_error("unsupported bit depth for reading");//these are going to be annoying to implement properly and probably aren't needed
			case  8:                              break;
			case 16: rowBytes *= 2;               break;
			case 24: rowBytes *= 3;               break;
			case 32: rowBytes *= 4;               break;
			default: throw std::runtime_error("unsupported bits per pixel");
		}
		const size_t padBytes = (4 - (rowBytes % 4)) % 4;//rows in the file are padded to multiples of 32 bits
		const size_t bmpRowBytes = rowBytes + padBytes;

		//read data (if we made it this far we have 8, 16, 24, or 32 bit pixels)
		is.seekg(offset);
		uint32_t pad;
		char       * pOut = buf;
		if(gry && (bitCount > 16)) {
			//read only the red channel
			const size_t stride = bitCount / 8;
			for(uint32_t j = 0; j < height; j++) {
				for(uint32_t i = 0; i < width; i++) {
					is.read(pOut + i, 1);
					is.read((char*)&pad, stride - 1);
				}
				is.read((char*)&pad, padBytes);
				pOut +=    width;
			}
		} else {
			//read raw data 
			for(uint32_t j = 0; j < height; j++) {
				is.read(pOut, rowBytes);
				is.read((char*)&pad, padBytes);
				pOut +=    rowBytes;
			}
		}
		return is;
	}

	//@brief    : read the image into a raw buffer
	//@param ptr: raw buffer to read from
	//@param buf: location to write image
	//@param gry: true to convert to grayscale (red channel for color, cast for 1/4 bit)
	void Header::readImage(char const * const ptr, char * const buf, const bool gry) const {
		//make sure we can read this bitmap
		if(0 != compression) throw std::runtime_error("unsupported compression");

		//compute row padding
		uint32_t rowBytes = width;
		switch(bitCount) {
			case  0: rowBytes = 0;               return;
			case  1: throw std::runtime_error("unsupported bit depth for reading");//these are going to be annoying to implement properly and probably aren't needed
			case  4: throw std::runtime_error("unsupported bit depth for reading");//these are going to be annoying to implement properly and probably aren't needed
			case  8:                              break;
			case 16: rowBytes *= 2;               break;
			case 24: rowBytes *= 3;               break;
			case 32: rowBytes *= 4;               break;
			default: throw std::runtime_error("unsupported bits per pixel");
		}
		const size_t padBytes = (4 - (rowBytes % 4)) % 4;//rows in the file are padded to multiples of 32 bits
		const size_t bmpRowBytes = rowBytes + padBytes;

		//read data (if we made it this far we have 8, 16, 24, or 32 bit pixels)
		char const * pIn  = ptr + offset;
		char       * pOut = buf;
		if(gry && (bitCount > 16)) {
			//read only the red channel
			const size_t stride = bitCount / 8;
			for(uint32_t j = 0; j < height; j++) {
				for(uint32_t i = 0; i < width; i++) pOut[i] = pIn[i * stride];
				pIn  += bmpRowBytes;
				pOut +=    width   ;
			}
		} else {
			//read raw data 
			for(uint32_t j = 0; j < height; j++) {
				std::copy(pIn, pIn + rowBytes, pOut);
				pIn  += bmpRowBytes;
				pOut +=    rowBytes;
			}
		}
	}

	//@brief   : read BITMAPCOREHEADER from an istream
	//@param is: istream to read from
	//@return  : is
	std::istream& Header::readCore(std::istream& is) {
		uint16_t tmp = 0;
		is.read((char*)&tmp     , 2); width  = tmp;
		is.read((char*)&tmp     , 2); height = tmp;
		is.read((char*)&planes  , 2);
		is.read((char*)&bitCount, 2);
		if(planes != 1) throw std::runtime_error("bitmap planes must be 1");
		if(!(bitCount ==  1 ||
		     bitCount ==  4 ||
		     bitCount ==  8 ||
		     bitCount == 24 )) throw std::runtime_error("bitmap bitcount must be 1, 4, 8, or 24");
		return is;
	}

	//@brief   : read BITMAPINFOHEADER from an istream
	//@param is: istream to read from
	//@return  : is
	std::istream& Header::readInfo(std::istream& is) {
		//read modified core header
		is.read((char*)&width   , 4);
		is.read((char*)&height  , 4);
		is.read((char*)&planes  , 2);
		is.read((char*)&bitCount, 2);
		if(planes != 1) throw std::runtime_error("bitmap planes must be 1");
		if(!(bitCount ==  0 ||
		     bitCount ==  1 ||
		     bitCount ==  4 ||
		     bitCount ==  8 ||
		     bitCount == 16 ||
		     bitCount == 24 ||
		     bitCount == 32 )) throw std::runtime_error("bitmap bitcount must be 1, 4, 8. 16, 24, or 32");

		//read extra fields
		is.read((char*)&compression, 4);
		is.read((char*)&sizeImage  , 4);
		is.read((char*)&xRes       , 4);
		is.read((char*)&yRes       , 4);
		is.read((char*)&colorsUsed , 4);
		is.read((char*)&colorsImprt, 4);
		return is;
	}

	//@brief   : read BITMAPV4HEADER from an istream
	//@param is: istream to read from
	//@return  : is
	std::istream& Header::readV4(std::istream& is) {
		//read BITMAPINFOHEADER
		readInfo(is);

		//read extra fields
		is.read((char*)&redMask   , 4);
		is.read((char*)&greenMask , 4);
		is.read((char*)&blueMask  , 4);
		is.read((char*)&alphaMask , 4);
		is.read((char*)&colorSpace, 4);
		is.read((char*)&redX      , 4);
		is.read((char*)&redY      , 4);
		is.read((char*)&redZ      , 4);
		is.read((char*)&greenX    , 4);
		is.read((char*)&greenY    , 4);
		is.read((char*)&greenZ    , 4);
		is.read((char*)&blueX     , 4);
		is.read((char*)&blueY     , 4);
		is.read((char*)&blueZ     , 4);
		is.read((char*)&gammaRed  , 4);
		is.read((char*)&gammaGreen, 4);
		is.read((char*)&gammaBlue , 4);
		return is;
	}

	//@brief   : read BITMAPV5HEADER from an istream
	//@param is: istream to read from
	//@return  : is
	std::istream& Header::readV5(std::istream& is) {
		//read BITMAPV4HEADER
		readInfo(is);

		//read extra fields
		is.read((char*)&intent     , 4);
		is.read((char*)&profileData, 4);
		is.read((char*)&profileSize, 4);
		is.read((char*)&reserved   , 4);
		return is;
	}

	//@brief: byteswap to little endian from big if needed
	void Header::headerToLittle() {
		if(detail::bigEndian()) {
			fileSize = detail::swapU32(fileSize);
			res1     = detail::swapU16(res1    );
			res2     = detail::swapU16(res2    );
			offset   = detail::swapU32(offset  );
		}
	}

	//@brief: byteswap to little endian from big if needed
	void Header::infoToLittle() {
		if(detail::bigEndian()) {
			switch(size) {
				case 124://BITMAPV5HEADER
					// BITMAPV5HEADER (124 bytes) adds these fields and support for additional colorSpace types
					intent      = detail::swapU32(intent     );
					profileData = detail::swapU32(profileData);
					profileSize = detail::swapU32(profileSize);
					reserved    = detail::swapU32(reserved   );
					//intentional fall through

				case 108://BITMAPV4HEADER
					// BITMAPV4HEADER (108 bytes) adds these fields
					redMask     = detail::swapU32(redMask    );
					greenMask   = detail::swapU32(greenMask  );
					blueMask    = detail::swapU32(blueMask   );
					alphaMask   = detail::swapU32(alphaMask  );
					colorSpace  = detail::swapU32(colorSpace );
					redX        = detail::swapU32(redX       );
					redY        = detail::swapU32(redY       );
					redZ        = detail::swapU32(redZ       );
					greenX      = detail::swapU32(greenX     );
					greenY      = detail::swapU32(greenY     );
					greenZ      = detail::swapU32(greenZ     );
					blueX       = detail::swapU32(blueX      );
					blueY       = detail::swapU32(blueY      );
					blueZ       = detail::swapU32(blueZ      );
					gammaRed    = detail::swapU32(gammaRed   );
					gammaGreen  = detail::swapU32(gammaGreen );
					gammaBlue   = detail::swapU32(gammaBlue  );
					//intentional fall through

				case  40://BITMAPINFOHEADER
					// BITMAPINFOHEADER (40 bytes) adds these fields and extends bitCount to include (0, 16, and 32)
					compression = detail::swapU32(compression);
					sizeImage   = detail::swapU32(sizeImage  );
					xRes        = detail::swapU32(xRes       );
					yRes        = detail::swapU32(yRes       );
					colorsUsed  = detail::swapU32(colorsUsed );
					colorsImprt = detail::swapU32(colorsImprt);

					// modified BITMAPCOREHEADER (16 bytes) (actual BITMAPCOREHEADER is 12 bytes)
					width       = detail::swapU32(width      );
					height      = detail::swapU32(height     );
					planes      = detail::swapU16(planes     );
					bitCount    = detail::swapU16(bitCount   );
					break;

				case  12://BITMAPCOREHEADER
					width       = detail::swapU16(width      );
					height      = detail::swapU16(height     );
					planes      = detail::swapU16(planes     );
					bitCount    = detail::swapU16(bitCount   );
					break;

				default: throw std::runtime_error("unsupported header size");
			}
		}
	}
}

#endif//_BMP_H_
