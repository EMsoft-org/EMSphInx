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

#include <ostream>
#include <stdexcept>
#include <string>

namespace base64 {

	//@brief      : base64 encode some data
	//@param ptr  : data to base64 encode
	//@param count: number of bytes to encode
	//@param os   : location to write encoded data
	//@return     : bytes written to os
	size_t encode(char const * ptr, size_t count, std::ostream& os) {
		static const char b64[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
		for(size_t i = 0; i < count; i+=3) {
			//get next 3 bytes
			char b[3] = {
				                ptr[i  ]          ,
				(i+1 < count) ? ptr[i+1] : char(0),//dont go out of bounds
				(i+2 < count) ? ptr[i+2] : char(0),//dont go out of bounds
			};

			//turn into 4 base64 indices
			int c[4] = {
				                       ((b[0] & 0xFC) >> 2),
				((b[0] & 0x03) << 4) | ((b[1] & 0xF0) >> 4),
				((b[1] & 0x0F) << 2) | ((b[2] & 0xC0) >> 6),
				((b[2] & 0x3F) << 0)
			};

			//write to output
			os.put(b64[c[0]]);
			os.put(b64[c[1]]);
			if(i+3 < count || i+3 == count) {
				os.put(b64[c[2]]);
				os.put(b64[c[3]]);
			} else if(i+1 < count) {
				os.put(b64[c[2]]);
				os.put('=');
			} else {
				os.put('=');
				os.put('=');
			}
		}
		return 4 * ((count + 2) / 3);//number of output bytes
	}

	//@brief      : base64 decode some data
	//@param ptr  : data to base64 decode
	//@param count: number of bytes to decode
	//@param os   : location to write encoded data
	//@return     : bytes written to os
	size_t decode(char const * ptr, size_t count, std::ostream& os) {
		if(0 != count % 4) throw std::runtime_error("base64 encoded data must be multiple of 4 bytes long (pad with '=' if too short)");
		const size_t maxBytes = (count / 4) * 3;//this is the number of encoded bytes if there was no padding required
		for(size_t i = 0; i < count; i+=4) {//loop over 4 character chungs (24 bits encoded as 32 bits)
			//extract the next block of characters
			const char block[4] = {ptr[i], ptr[i+1], ptr[i+2], ptr[i+3]};

			//convert characters to indices
			char c[4];
			for(size_t j = 0; j < 4; j++) {
				//this could probably use a lookup table but that might be wrong (if the compiler isn't ascii, admittedly unlikely)
				switch(block[j]) {
					case 'A': c[j] = 0x00; break;
					case 'B': c[j] = 0x01; break;
					case 'C': c[j] = 0x02; break;
					case 'D': c[j] = 0x03; break;
					case 'E': c[j] = 0x04; break;
					case 'F': c[j] = 0x05; break;
					case 'G': c[j] = 0x06; break;
					case 'H': c[j] = 0x07; break;
					case 'I': c[j] = 0x08; break;
					case 'J': c[j] = 0x09; break;
					case 'K': c[j] = 0x0A; break;
					case 'L': c[j] = 0x0B; break;
					case 'M': c[j] = 0x0C; break;
					case 'N': c[j] = 0x0D; break;
					case 'O': c[j] = 0x0E; break;
					case 'P': c[j] = 0x0F; break;
					case 'Q': c[j] = 0x10; break;
					case 'R': c[j] = 0x11; break;
					case 'S': c[j] = 0x12; break;
					case 'T': c[j] = 0x13; break;
					case 'U': c[j] = 0x14; break;
					case 'V': c[j] = 0x15; break;
					case 'W': c[j] = 0x16; break;
					case 'X': c[j] = 0x17; break;
					case 'Y': c[j] = 0x18; break;
					case 'Z': c[j] = 0x19; break;
					case 'a': c[j] = 0x1A; break;
					case 'b': c[j] = 0x1B; break;
					case 'c': c[j] = 0x1C; break;
					case 'd': c[j] = 0x1D; break;
					case 'e': c[j] = 0x1E; break;
					case 'f': c[j] = 0x1F; break;
					case 'g': c[j] = 0x20; break;
					case 'h': c[j] = 0x21; break;
					case 'i': c[j] = 0x22; break;
					case 'j': c[j] = 0x23; break;
					case 'k': c[j] = 0x24; break;
					case 'l': c[j] = 0x25; break;
					case 'm': c[j] = 0x26; break;
					case 'n': c[j] = 0x27; break;
					case 'o': c[j] = 0x28; break;
					case 'p': c[j] = 0x29; break;
					case 'q': c[j] = 0x2A; break;
					case 'r': c[j] = 0x2B; break;
					case 's': c[j] = 0x2C; break;
					case 't': c[j] = 0x2D; break;
					case 'u': c[j] = 0x2E; break;
					case 'v': c[j] = 0x2F; break;
					case 'w': c[j] = 0x30; break;
					case 'x': c[j] = 0x31; break;
					case 'y': c[j] = 0x32; break;
					case 'z': c[j] = 0x33; break;
					case '0': c[j] = 0x34; break;
					case '1': c[j] = 0x35; break;
					case '2': c[j] = 0x36; break;
					case '3': c[j] = 0x37; break;
					case '4': c[j] = 0x38; break;
					case '5': c[j] = 0x39; break;
					case '6': c[j] = 0x3A; break;
					case '7': c[j] = 0x3B; break;
					case '8': c[j] = 0x3C; break;
					case '9': c[j] = 0x3D; break;
					case '+': c[j] = 0x3E; break;
					case '/': c[j] = 0x3F; break;
					case '=': {//termination
						if(i+4 < count || j < 2) throw std::runtime_error("unexpected base64 terminator");//make sure this didn't show up too early
						if(2 == j && block[3] != '=') throw std::runtime_error("character after first `='");//make sure we have either ..== or ...= (not ..=.)
						c[j] = 0x00;//fill with 0 so we don't disrupt other bytes
					} break;
					default : throw std::runtime_error(std::string("unexpected base64 character `") + ptr[i] + "'");
				}
			}

			//decode and 3 bytes
			const int b[3] = {
				((c[0] << 2) & 0xFC ) | ( (c[1] >> 4) & 0x03 ),
				((c[1] << 4) & 0xF0 ) | ( (c[2] >> 2) & 0x0F ),
				((c[2] << 6) & 0xC0 ) | ( (c[3]     ) & 0x3F )
			};

			//save extracted bytes
			os.put(b[0]);//we always save the first byte
			if('=' == block[3]) {//don't save all bytes
				if('=' == block[2]) {//only 1 byte was encoded
					return maxBytes - 2;
				} else {//only 2 bytes were incoded
					os.put(b[1]);
					return maxBytes - 1;
				}
			} else {//this isn't the end or there were not padding bytes
				os.put(b[1]);
				os.put(b[2]);
			}
		}
		return maxBytes;//for no padding on last block
	}
}

