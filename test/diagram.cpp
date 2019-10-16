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

#include "xtal/diagram.hpp"

int main(int argc, char *argv[]) {
	if(2 == argc) {//make a diagram from a master pattern file
		std::string fileName(argv[1]);
		emsphinx::MasterPattern<double> mp(fileName);//read master pattern
		svg::Color c(0, 0, 0);
		xtal::Diagram diag(mp, c);
		diag.getHemi(true ).write("north.svg");
		diag.getHemi(false).write("south.svg");
	} else {//make plain diagrams
		//enumerate point group names to construct
		std::vector<std::string> groups = {
			    "1",
			   "-1",
			  "121",//multiple settings
			  "112",//multiple settings
			  "1m1",//multiple settings
			  "11m",//multiple settings
			"12/m1",//multiple settings
			"112/m",//multiple settings
			  "222",
			  "mm2",
			  "mmm",
			    "4",
			   "-4",
			  "4/m",
			  "422",
			  "4mm",
			 "-42m",// multiple settings
			 "-4m2",// multiple settings
			"4/mmm",
			    "3",
			   "-3",
			  "321",// multiple settings
			  "312",// multiple settings
			  "3m1",// multiple settings
			  "31m",// multiple settings
			 "-3m1",// multiple settings
			 "-31m",// multiple settings
			    "6",
			   "-6",
			  "6/m",
			  "622",
			  "6mm",
			 "-6m2",// multiple settings
			 "-62m",// multiple settings
			"6/mmm",
			   "23",
			   "m3",
			  "432",
			 "-43m",
			  "m3m"
		};

		//loop over point groups building diagram
		for(std::string& pg : groups) {
			std::string name = pg;
			replace(name.begin(), name.end(), '/', '_');//dont allow slashes in file name
			xtal::Diagram dg(pg);
			dg.addLabel(pg);
			dg.getHemi().write(name + ".svg");
		}
	}

	return 0;
}
