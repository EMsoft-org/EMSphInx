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


#ifndef _emsoft_h_
#define _emsoft_h_

#include <vector>
#include <type_traits>
#include <iostream>

#include "xtal/rotations.hpp"

namespace emsoft {

	template <typename Real>
	struct AngleFile {

		//@brief : get the format of the stored rotations
		//@return: format as one of Rotation enumerate
		xtal::Rotation getType() const {return rot;}

		//@brief : get the number of rotations
		//@return: number of stored rotations
		size_t getNum() const {return ang.size() / xtal::rotLen(rot);}

		//@brief  : get the ith rotation
		//@param i: rotation to get (no bounds check)
		//@return : pointer to ith rotation
		Real const * operator[](const size_t i) const {return &ang[i*xtal::rotLen(rot)];}

		//@brief : get the stored rotations
		//@return: pointer to rotations
		Real const * getAngles() const {return ang.data();}

		//@brief   : read angle angle file
		//@param is: input stream to read from
		void read(std::istream& is);

		//@brief   : write an angle file
		//@param os: output stream to write to
		void write(std::ostream& os) const;

		//@brief     : read angle angle file
		//@param name: name of file to read from
		void read(const char * name);
		void read(std::string& name) {read(name.c_str());}

		//@brief     : write an angle file
		//@param name: name of file to write to
		void write(const char * name) const;
		void write(std::string& name) const {write(name.c_str());}

		//@brief: create an empty angle file
		AngleFile() {}

		//@brief     : build an AngleFile from an input
		//@param name: name of file to parse
		AngleFile(const char* name) {read(name);}
		AngleFile(std::string name) {read(name);}

		private:
			static_assert(std::is_floating_point<Real>::value, "AngleFile must be templated on a floating point type");
			std::vector<Real> ang;//actual orientations
			xtal::Rotation    rot;//orientation representation type
	};
}

////////////////////////////////////////////////////////////////////////////////
//                           Implementation Details                           //
////////////////////////////////////////////////////////////////////////////////

#include <stdexcept>
#include <iterator>
#include <locale>

namespace emsoft {
	namespace detail {
		//a custom local that treats commas as white space
		//modified from https://stackoverflow.com/questions/1894886/parsing-a-comma-delimited-stdstring
		struct csv_reader: std::ctype<char> {
			csv_reader(): std::ctype<char>(get_table()) {}

			static std::ctype_base::mask const* get_table() {
				static std::ctype_base::mask const * const pMask = std::ctype<char>::classic_table();//get the default type table
				static std::vector<std::ctype_base::mask> rc(pMask, pMask + std::ctype<char>::table_size);//copy the table
				rc[','] = std::ctype_base::space;//treat comma as white space
				return rc.data();
			}
		};
	}

	//@brief   : read angle angle file
	//@param is: input stream to read from
	template <typename Real>
	void AngleFile<Real>::read(std::istream& is) {
		//try to get name and number
		size_t num = 0;
		if(!(is >> rot >> num)) throw std::runtime_error("failed to read name and number of orientations from angle file");
		if(xtal::Rotation::Unknown == rot) throw std::runtime_error("first line of angle file must be one of 'eu', 'om', 'ax', 'ro', 'qu', 'ho', or 'cu'");
		num *= xtal::rotLen(rot);

		//read all available numbers into the angle vector
		ang.reserve(num * xtal::rotLen(rot));//alloate space up front
		is.imbue(std::locale(std::locale(), new detail::csv_reader()));//treat commas as white space
		std::copy(std::istream_iterator<Real>(is), std::istream_iterator<Real>(), std::back_inserter(ang));//this one liner only w

		//make sure the size is what we expected
		if(ang.size() < num) throw std::runtime_error("not enough orientions in angle file");
		if(ang.size() > num) throw std::runtime_error("too many orientions in angle file");

		//convert from degrees to radians if needed
		if(xtal::Rotation::Euler == rot) std::for_each(ang.begin(), ang.end(), [](Real& v){v *= xtal::Constants<Real>::dg2rd;});
	}

	//@brief   : write an angle file
	//@param os: output stream to write to
	template <typename Real>
	void AngleFile<Real>::write(std::ostream& os) const {
		//write type/number
		os << rot << '\n' << getNum() << '\n';

		//write orientations
		const size_t len = xtal::rotLen(rot);
		const size_t num = getNum();
		Real const * p = ang.data();
		if(xtal::Rotation::Euler == rot) {//euler angles (need to convert from radians to degree)
			for(size_t i = 0; i < num; i++) {
				os << p[0] * xtal::Constants<Real>::rd2dg;
				for(size_t j = 1; j < len; j++) os << ',' << p[j] * xtal::Constants<Real>::rd2dg;
				os << '\n';
				p += len;
			}
		} else {//normal representatino
			for(size_t i = 0; i < num; i++) {
				os << p[0];
				for(size_t j = 1; j < len; j++) os << ',' << p[j];
				os << '\n';
				p += len;
			}
		}
	}

	//@brief     : read angle angle file
	//@param name: name of file to read from
	template <typename Real>
	void AngleFile<Real>::read(const char * name) {
		std::ifstream is(name);
		read(is);
	}

	//@brief     : write an angle file
	//@param name: name of file to write to
	template <typename Real>
	void AngleFile<Real>::write(const char * name) const {
		std::ofstream os(name);
		write(os);
	}
}

#endif//_emsoft_h_
