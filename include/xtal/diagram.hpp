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

#ifndef _xtal_diagram_h_
#define _xtal_diagram_h_

#include "util/svg.hpp"
#include "sht/square_sht.hpp"//lambert::sphereToSquare
#include "idx/master.hpp"//MasterPattern
#include "util/image.hpp"//interpolation
#include "symmetry.hpp"

namespace xtal {
	struct Diagram {
		//supported projection types (currently only circular projections)
		enum class Type {
			Stereo ,//stereographic (conformal)
			Lambert,//lambert azimuthal (area preserving)
		};

		//@brief   : construct a crystal diagram
		//@param t : type of projection
		//@param d : size of projection
		Diagram(const Type t = Type::Stereo, const size_t d = 100);

		//@brief   : construct a crystal diagram for a specific point group
		//@param pg: point group to construct diagram for
		//@param t : type of projection
		Diagram(PointGroup pg, const Type t = Type::Stereo) : Diagram(t) {add(pg);}

		//@brief   : construct a crystal diagram for a master pattern
		//@param mp: master pattern to construct diagram for
		//@param c : color for symmetry elements
		//@param t : type of projection
		template <typename Real> Diagram(emsphinx::MasterPattern<Real>& mp, svg::Color c = svg::Color(0.375, 0, 0), const Type t = Type::Stereo);

		//@brief   : construct a crystal diagram for an arbitrary color function
		//@param pg: point group to construct diagram for
		//@param cf: coloring function that converts from a unit crystallographic direction to an rgb color (rgb in [0,1])
		//@param c : color for symmetry elements
		//@param t : type of projection
		template <typename Real> Diagram(PointGroup pg, std::function<void(Real const*const,Real *const)> cf, svg::Color c = svg::Color(0, 0, 0), const Type t = Type::Stereo);

		//@brief   : construct an IPF colored crystal diagram
		//@param pg: point group to construct diagram for
		//@param t : type of projection
		static Diagram IpfLegend(PointGroup pg, const Type t = Type::Stereo);

		//@brief  : set the color to use for drawing elements
		//@param r: red   [0,1]
		//@param g: green [0,1]
		//@param b: blue  [0,1]
		void setColor(double r, double g, double b) {clr.setRGB(r, g, b);}

		//@brief    : add a rotational symmetry operator
		//@param x  : x coordinate of rotation axis
		//@param y  : y coordinate of rotation axis
		//@param z  : z coordinate of rotation axis
		//@param fld: order of rotational symmetry axis, negative for inversion symmetry
		//@param scl: scale factor for symbol size ([0.5, 1.5] is reasonable)
		//@note     : currently only 2, 3, 4, 6, -1, -2, -3, -4, and -6 are supported
		void addRotor(double x, double y, double z, int fld, double scl = 1);

		//@brief  : add a mirror plane
		//@param x: x coordinate of rotation axis
		//@param y: y coordinate of rotation axis
		//@param z: z coordinate of rotation axis
		void addMirror(double x, double y, double z);

		//@brief: add an inversion symmetry marker
		//@param scl: scale factor for symbol size ([0.5, 1.5] is reasonable)
		void addInversion(double scl = 1);

		//@brief    : add the symmetry operators associated with a point group
		//@param scl: scale factor for symbol size ([0.5, 1.5] is reasonable)
		//@param pg : point group to add symmetry operators of
		void add(PointGroup pg, const double scl = 1);

		//@brief    : add a dot to the projections
		//@paran n  : direction to add spot for (+/-), will be normalized
		//@param c  : color to use
		//@param scl: scale factor for symbol size ([0.5, 1.5] is reasonable)
		void addDot(double const * const n, svg::Color c, const double scl = 1);

		//@brief    : add a label to the projections
		//@param lbl: label to add
		//@param scl: scale factor for symbol size ([0.5, 1.5] is reasonable)
		void addLabel(std::string lbl, const double scl = 1);

		//@brief    : get the current svg for a hemisphere
		//@param nth: true/false for north/south hemisphere
		//@return   : svg for selected hemisphere
		svg::SVG getHemi(const bool nth = true) const {return nth ? nh : sh;}

		private:
			const size_t dim ;//canvas size (could be a floating point type if needed...)
			Type         type;//diagram projection type
			svg::SVG     nh  ;//diagram of north hemisphere as svg
			svg::SVG     sh  ;//diagram of north hemisphere as svg
			svg::Color   clr ;//color to use for markers / lines

			//@brief    : project a unit direction to the canvas
			//@param x  : x coordinate of rotation axis
			//@param y  : y coordinate of rotation axis
			//@param z  : z coordinate of rotation axis
			//@param X  : X projected coordinate in canvas
			//@param Y  : Y projected coordinate in canvas
			void project(double x, double y, double z, double& X, double& Y) const;

			//@brief    : unproject a vector in the unit circle to the sphere
			//@param X  : X projected coordinate
			//@param Y  : Y projected coordinate
			//@param x  : x coordinate of unit direction on sphere
			//@param y  : y coordinate of unit direction on sphere
			//@param z  : z coordinate of unit direction on sphere
			//@return   : true/false if XY was a valid direction (in the north hemisphere)
			bool unproject(double X, double Y, double& x, double& y, double& z) const;

			//@brief: convert a Catmull-Rom spline to a cubic Bezier curve
			//@param p: Catmull-Rom control points (4)
			//@param t: Catmull-Rom knots (4)
			//@param c: location to write cubic spline control points (4)
			//@note   : c[0] = p[1], c[3] = p[2]
			static void CatmullRom2Bezier(double const * const p, double const * const t, double * const c);
	};
}

#include <cmath>

#include "constants.hpp"

namespace xtal {
	//@brief   : construct a crystal diagram
	//@param t : type of projection
	//@param d : size of projection
	Diagram::Diagram(const Type t, const size_t d) : dim(d), type(t), nh((double)dim, (double)dim), sh((double)dim, (double)dim), clr(0, 0, 0) {
		const double dDim = (double)dim;
		svg::Circle eq(dDim/2, dDim/2, dDim/2);
		eq.setStroke(svg::Stroke(dDim / 500));
		nh.add(eq);
		sh.add(eq);
	}

	//@brief   : construct a crystal diagram for a master pattern
	//@param mp: master pattern to construct diagram for
	//@param c : color for symmetry elements
	//@param t : type of projection
	template <typename Real>
	Diagram::Diagram(emsphinx::MasterPattern<Real>& mp, svg::Color c, const Type t) : Diagram(t, mp.dim) {
		//create image for north/south hemispheres
		svg::Image imNh(dim, dim, svg::Image::GrayAlpha);
		svg::Image imSh(dim, dim, svg::Image::GrayAlpha);

		//determine pixel range so we can rescale to [0,255]
		auto minMaxNh = std::minmax_element(mp.nh.begin(), mp.nh.end());
		auto minMaxSh = std::minmax_element(mp.sh.begin(), mp.sh.end());
		const Real vMin = std::min(*minMaxNh.first , *minMaxSh.first );
		const Real vMax = std::max(*minMaxNh.second, *minMaxSh.second);
		const Real vDel = vMax - vMin;
		std::function<uint8_t(const Real&)> to8Bit = [&](const Real& v){return (uint8_t)std::round((v - vMin) * Real(255) / vDel);};

		//interpolate from square lambert to stereographic
		double n[3], ix, iy;
		image::BiPix<Real> pix;
		for(size_t j = 0; j < dim; j++) {
			const double Y = -((Real(j) / (dim - 1)) * 2 - 1);//[-1,1], negate for image convention
			for(size_t i = 0; i < dim; i++) {
				const double X = (Real(i) / (dim - 1)) * 2 - 1;//[-1,1]
				if(unproject(X, Y, n[0], n[1], n[2])) {//unstereographic project
					emsphinx::square::lambert::sphereToSquare(n[0], n[1], n[2], ix, iy);//square lambert project
					pix.bilinearCoeff(ix, iy, dim, dim);//compute bilinear interpolation in square lambert image
					const size_t idxIm  = j * dim + i;//get index of output image
					imNh.buff[2*idxIm + 0] = to8Bit(pix.interpolate(mp.nh.data()));
					imNh.buff[2*idxIm + 1] = 0xFF;//alpha
					imSh.buff[2*idxIm + 0] = to8Bit(pix.interpolate(mp.sh.data()));
					imSh.buff[2*idxIm + 1] = 0xFF;//alpha
				}
			}
		}

		//add images
		nh.add(imNh);
		sh.add(imSh);

		//add symmetry
		clr = c;
		add(mp.pointGroup(), 0.67);
	}

	//@brief   : construct a crystal diagram for an arbitrary color function
	//@param pg: point group to construct diagram for
	//@param cf: coloring function that converts from a unit crystallographic direction to an rgb color (rgb in [0,1])
	//@param c : color for symmetry elements
	//@param t : type of projection
	template <typename Real>
	Diagram::Diagram(PointGroup pg, std::function<void(Real const*const,Real *const)> cf, svg::Color c, const Type t) : Diagram(t, 512) {
		//create image for north/south hemispheres
		svg::Image imNh(dim, dim, svg::Image::RGBA);
		svg::Image imSh(dim, dim, svg::Image::RGBA);
		//loop over pixels of images computing color
		double n[3];
		for(size_t j = 0; j < dim; j++) {
			const double Y = -((Real(j) / (dim - 1)) * 2 - 1);//[-1,1], negate for image convention
			for(size_t i = 0; i < dim; i++) {
				const double X = (Real(i) / (dim - 1)) * 2 - 1;//[-1,1]
				if(unproject(X, Y, n[0], n[1], n[2])) {//unstereographic project
					//compute colors
					Real rgbNh[3], rgbSh[3];
					cf(n, rgbNh);//compute color for north hemisphere
					n[2] = -n[2];
					cf(n, rgbSh);//compute color for south hemisphere

					//convert to 8 bit
					for(size_t k = 0; k < 3; k++) {
						rgbNh[k] = std::round(rgbNh[k] * 255);
						rgbSh[k] = std::round(rgbSh[k] * 255);
					}

					//save result
					const size_t idxIm  = (j * dim + i) * 4;//get index of output image
					for(size_t k = 0; k < 3; k++) imNh.buff[idxIm + k] = (uint8_t) rgbNh[k];
					imNh.buff[idxIm + 3] = 0xFF;//alpha
					for(size_t k = 0; k < 3; k++) imSh.buff[idxIm + k] = (uint8_t) rgbSh[k];
					imSh.buff[idxIm + 3] = 0xFF;//alpha
				}
			}
		}

		//add images
		nh.add(imNh);
		sh.add(imSh);

		//add symmetry
		clr = c;
		add(pg, 0.67);
	}

	//@brief   : construct an IPF colored crystal diagram
	//@param pg: point group to construct diagram for
	//@param t : type of projection
	Diagram Diagram::IpfLegend(PointGroup pg, const Type t) {
		std::function<void(double const*const,double *const)> cf = std::bind(&PointGroup::ipfColor<double>, &pg, std::placeholders::_1, std::placeholders::_2, sph2rgb<double>);
		return Diagram(pg, cf, svg::Color(0, 0, 0), t);
	}

	//@brief    : add a rotational symmetry operator
	//@param x  : x coordinate of rotation axis
	//@param y  : y coordinate of rotation axis
	//@param z  : z coordinate of rotation axis
	//@param fld: order of rotational symmetry axis, negative for inversion symmetry
	//@param scl: scale factor for symbol size ([0.5, 1.5] is reasonable)
	//@note     : currently only 2, 3, 4, 6, -1, -2, -3, -4, and -6 are supported
	void Diagram::addRotor(double x, double y, double z, int fld, double scl) {
		//compute scaled constant
		const int aFld = std::abs(fld);
		const double rad = double(dim) * scl / 20;

		//project to hemisphere and compute rotation
		double nX, nY, sX, sY;
		if(std::signbit(z)) {
			x = -x;
			y = -y;
			z = -z;
		}
		project( x,  y,  z, nX, nY);
		project(-x, -y, -z, sX, sY);
		const bool rot2Fold = false;//should lenses be perpendicular (true) or parallel (false) to the equator
		double nTheta = 90 + std::atan2(nY - double(dim) / 2, nX - double(dim) / 2) * 180.0 / Constants<double>::pi;
		double sTheta = 90 + std::atan2(sY - double(dim) / 2, sX - double(dim) / 2) * 180.0 / Constants<double>::pi;
		if(std::max(std::fabs(x), std::fabs(y)) < 1e-6) {//don't rotate if we're extremely close to the pole
			if((aFld == 2 || fld == -4) && !rot2Fold) {
				nTheta = sTheta = 90;//except for lens which I'll always keep vertical at the poles
			} else {
				nTheta = sTheta = 0;
			}
		}

		//build base shapes
		double rx = rad;
		double ry = rad / 2.5;
		if(rot2Fold) std::swap(rx, ry);
		svg::Circle  cirK = svg::Circle(rad / 3);                //1   fold symbol
		svg::Ellipse ell  = svg::Ellipse(rx, ry);                //2   fold symbol
		svg::Polygon ply  = svg::Polygon::Regular(aFld   , rad); //n   fold symbol
		svg::Polygon ply2 = svg::Polygon::Regular(aFld / 2, rad);//n/2 fold symbol
		svg::Circle  cirW = svg::Circle(rad / 5);                //inversion symbol

		//lens + polygon black fill only
		cirK.setFill(clr).removeStroke();
		ell .setFill(clr).removeStroke();
		ply .setFill(clr).removeStroke();
		ply2.setFill(clr).removeStroke();

		//make dot white fill (inversion marker)
		cirW.setFill(1, 1, 1).removeStroke();

		//assemble marker
		svg::Group mrk;
		if(fld < 0) {
			//inversion symmetry
			if(aFld < 3) {
				if(2 == aFld) {
					mrk.add(ell);//use lens for 2 fold
					mrk.add(cirW);//add white circle to center
				} else {
					throw std::runtime_error("rotoinversion axis must be at least -2");
				}
			} else if(0 == aFld % 2) {
				//even means there is an aFld / 2 pure rotation, use special symbol
				ply.setFill(1, 1, 1);
				ply.setStrokeColor(clr);
				mrk.add(ply );
				4 == aFld ? mrk.add(ell) : mrk.add(ply2);
			} else {
				//odd means both, just add dot
				mrk.add(ply );//use lens for 2 fold
				mrk.add(cirW);//add white circle to center
			}
		} else {
			//no inversion symmetry
			if(fld < 3) {
				if(2 == fld) {
					mrk.add(ell);//use lens for 2 fold
				} else {
					throw std::runtime_error("cannot add 0 or 1 fold rotational symmetry");
				}
			} else {
				//use polygon for 3+ fold
				mrk.add(ply);
			}
		}

		//place marker in correct location handling equator
		if(z == 0) {
			nh.add(svg::Group(mrk).translate(nX, nY).rotate(nTheta));
			nh.add(svg::Group(mrk).translate(sX, sY).rotate(sTheta));
			sh.add(svg::Group(mrk).translate(nX, nY).rotate(nTheta));
			sh.add(svg::Group(mrk).translate(sX, sY).rotate(sTheta));
		} else {
			nh.add(svg::Group(mrk).translate(nX, nY).rotate(nTheta));
			sh.add(svg::Group(mrk).translate(sX, sY).rotate(sTheta));
		}
	}

	//@brief  : add a mirror plane
	//@param x: x coordinate of rotation axis
	//@param y: y coordinate of rotation axis
	//@param z: z coordinate of rotation axis
	void Diagram::addMirror(double x, double y, double z) {
		//normalize
		const double mag = std::sqrt(x * x + y * y + z * z);
		x /= mag;
		y /= mag;
		z /= mag;

		//check type of arc
		const double dDim = (double)dim;
		const double aZ = std::fabs(z);
		if(std::fabs(aZ - 1.0) < 1e-6) {//plane normal is z axis
			//add line at equator
			svg::Circle eq(dDim / 2, dDim / 2, dDim / 2);
			eq.setStroke(dDim / 100, clr);
			nh.add(eq);
			sh.add(eq);
		} else if(aZ < 1e-6) {//plane normal lies in equator
			double X1, Y1, X2, Y2;
			project(-y,  x, 0, X1, Y1);
			project( y, -x, 0, X2, Y2);
			svg::Line line(X1, Y1, X2, Y2);
			line.setStroke(dDim / 100, clr);
			nh.add(line);
			sh.add(line);
		} else {//general arc (circle in stereographic projection)
			//compute intersection with equator
			const double h = std::hypot(x, y);
			double n[2] = {y/h, -x/h};//intersection with equator at +/- n

			//calculate some points on the great circle
			std::vector<double> ctrl;//x nh, y nh, x sh, y sh
			const size_t extPts = 5;//how many extra control points should we insert
			const size_t numPts = extPts + 2;//internal control points + end points
			const size_t intPts = numPts - 1;//additional non path points (cubic control points shared by neighbors)
			const size_t totPts = numPts + intPts;
			for(size_t i = 0; i < totPts; i++) {
				//build rotation matrix
				const double theta = double(i) * Constants<double>::pi / (totPts - 1);//fractional angle around half of great circle
				const double c  = std::cos(theta);
				const double cc = 1.0 - c;
				const double s  = std::sin(theta);
				double mat[3][2] = {//we don't need the last column (z of the vector being rotated is 0)
					c + x * x * cc        ,     x * y * cc + z * s,
					    y * x * cc - z * s, c + y * y * cc        ,
					    z * x * cc + y * s,     z * y * cc - x * s,
				};

				//compute rotated vector
				double rn[3] = {
					mat[0][0] * n[0] + mat[0][1] * n[1],
					mat[1][0] * n[0] + mat[1][1] * n[1],
					mat[2][0] * n[0] + mat[2][1] * n[1],
				};

				//project
				double X1, Y1, X2, Y2;
				project( rn[0],  rn[1],  rn[2], X1, Y1);
				project(-rn[0], -rn[1], -rn[2], X2, Y2);

				//add to control points
				ctrl.push_back(X1);
				ctrl.push_back(Y1);
				ctrl.push_back(X2);
				ctrl.push_back(Y2);
			}

			//build arcs
			svg::Path nhArc(ctrl[0], ctrl[1]), shArc(ctrl[2], ctrl[3]);
			double t[4] = {0,1,2,3};
			double p[4][4] = {
				ctrl[0], ctrl[0], ctrl[0], ctrl[4],
				ctrl[1], ctrl[1], ctrl[1], ctrl[5],
				ctrl[2], ctrl[2], ctrl[2], ctrl[6],
				ctrl[3], ctrl[3], ctrl[3], ctrl[7],
			};
			double c[4][4];
			for(size_t i = 1; i < totPts; i++) {
				for(size_t j = 0; j < 4; j++) {//loop over the 4 points we need to convert splines for (X1, Y1, X2, Y2)
					std::rotate(p[j], p[j] + 1, p[j] + 4);//shift all points back one
					p[j][3] = ctrl[4 * (i + 1 == totPts ? i : i+1) + j];//grab new point
					CatmullRom2Bezier(p[j], t, c[j]);//compute cubic bezier curve that goes through our control points
				}
				nhArc.curveTo(c[0][1], c[1][1], c[0][2], c[1][2], c[0][3], c[1][3], true);//accumulate arc
				shArc.curveTo(c[2][1], c[3][1], c[2][2], c[3][2], c[2][3], c[3][3], true);//accumulate arc
			}
			nhArc.setStroke(dDim / 100, clr);
			shArc.setStroke(dDim / 100, clr);
			nh.add(nhArc);
			sh.add(shArc);
		}
	}

	//@brief: add an inversion symmetry marker
	//@param scl: scale factor for symbol size ([0.5, 1.5] is reasonable)
	void Diagram::addInversion(double scl) {
		//project to hemisphere
		double nX, nY, sX, sY;
		project(0, 0, 1, nX, nY);
		project(0, 0,-1, sX, sY);

		//build base shapes
		const double rad = double(dim) * scl / 20;
		svg::Circle  cirK = svg::Circle(rad / 3);//outer circle
		svg::Circle  cirW = svg::Circle(rad / 5);//inner circle
		cirK.setFill(clr).removeStroke();//solid dot
		cirW.setFill(1, 1, 1).removeStroke();//white dot

		//assemble marker and place
		svg::Group mrk;
		mrk.add(cirK);//use lens for 2 fold
		mrk.add(cirW);//add white circle to center
		nh.add(svg::Group(mrk).translate(nX, nY));
		sh.add(svg::Group(mrk).translate(sX, sY));
	}

	//@brief    : add the symmetry operators associated with a point group
	//@param pg : point group to add symmetry operators of
	//@param scl: scale factor for symbol size ([0.5, 1.5] is reasonable)
	void Diagram::add(PointGroup pg, const double scl) {
		//add mirror planes
		const uint_fast8_t numMir = pg.numMirror();
		double const * mir = pg.mirrors<double>();
		for(uint_fast8_t i = 0; i < numMir; i++) addMirror(mir[3*i+0], mir[3*i+1], mir[3*i+2]);

		//add rotors
		const uint_fast8_t numAx = pg.numRotAxis();
		double const * ax = pg.rotAxis<double>();
		for(uint_fast8_t i = 0; i < numAx; i++) addRotor(ax[4*i+1], ax[4*i+2], ax[4*i+3], (int)ax[4*i+0], scl);

		//add inversion center
		if(pg.inversion() && 3 != pg.zRot()) addInversion();//add an inversion dot if there isn't already one from the -3 axis
	}

	//@brief    : add a dot to the projections
	//@paran n  : direction to add spot for (+/-), will be normalized
	//@param c  : color to use
	//@param scl: scale factor for symbol size ([0.5, 1.5] is reasonable)
	void Diagram::addDot(double const * const n, svg::Color c, const double scl) {
		//project to hemisphere
		double nX, nY, sX, sY;
		double x = n[0], y = n[1], z = n[2];
		if(std::signbit(z)) {
			x = -x;
			y = -y;
			z = -z;
		}
		project( x,  y,  z, nX, nY);
		project(-x, -y, -z, sX, sY);

		//build dot
		svg::Circle dot = svg::Circle(double(dim) * scl / 60);
		dot.setFill(c).removeStroke();

		//place dot in correct location handling equator
		if(z == 0) {
			nh.add(dot.translate(nX, nY));
			nh.add(dot.translate(sX, sY));
			sh.add(dot.translate(nX, nY));
			sh.add(dot.translate(sX, sY));
		} else {
			nh.add(dot.translate(nX, nY));
			sh.add(dot.translate(sX, sY));
		}
	}

	//@brief    : add a label to the projections
	//@param lbl: label to add
	//@param scl: scale factor for symbol size ([0.5, 1.5] is reasonable)
	void Diagram::addLabel(std::string lbl, const double scl) {
		svg::Text tx;
		tx.text = lbl;
		tx.size = scl * double(dim) / 12;
		nh.add(tx);
		sh.add(tx);
	}

	//@brief    : project a unit direction to the canvas
	//@param x  : x coordinate of rotation axis
	//@param y  : y coordinate of rotation axis
	//@param z  : z coordinate of rotation axis
	//@param X  : X projected coordinate in canvas
	//@param Y  : Y projected coordinate in canvas
	void Diagram::project(double x, double y, double z, double& X, double& Y) const {
		//normalize
		const double mag = std::sqrt(x * x + y * y + z * z);
		x /= mag;
		y /= mag;
		z /= mag;

		//stereographic project to [-1,1]
		X = x / (std::fabs(z) + 1);
		Y = y / (std::fabs(z) + 1);

		//rescale to canvas dimensinos
		X = (1 + X) * double(dim) / 2;
		Y = (1 - Y) * double(dim) / 2;//correct for image coordinate system
	}

	//@brief    : unproject a vector in the unit circle to the sphere
	//@param X  : X projected coordinate
	//@param Y  : Y projected coordinate
	//@param x  : x coordinate of unit direction on sphere
	//@param y  : y coordinate of unit direction on sphere
	//@param z  : z coordinate of unit direction on sphere
	//@return   : true/false if XY was a valid direction (in the north hemisphere)
	bool Diagram::unproject(double X, double Y, double& x, double& y, double& z) const {
		const double h2 = X * X + Y * Y;
		if(h2 > 1.0) return false;
		const double den = h2 + 1.0;
		x = X * 2 / den;
		y = Y * 2 / den;
		z = (1.0 - h2) / den;
		return true;
	}

	//@brief: convert a Catmull-Rom spline to a cubic Bezier curve
	//@param p: Catmull-Rom control points (4)
	//@param t: Catmull-Rom knots (4)
	//@param c: location to write cubic spline control points (4)
	//@note   : c[0] = p[1], c[3] = p[2]
	void Diagram::CatmullRom2Bezier(double const * const p, double const * const t, double * const c) {
		//https://stackoverflow.com/questions/30748316/catmull-rom-interpolation-on-svg-paths
		const double t21 = t[2] - t[1];
		const double t20 = t[2] - t[0];
		const double t31 = t[3] - t[1];
		const double t10 = t[1] - t[0];
		const double t32 = t[3] - t[2];
		const double c1 = (t21 * t21) / (t20 * t10);
		const double c2 =  t10        /  t20       ;
		const double d1 =  t32        /  t31       ;
		const double d2 = (t21 * t21) / (t31 * t32);
		const double m1 = c1 * (p[1]-p[0]) + c2 * (p[2]-p[1]);
		const double m2 = d1 * (p[2]-p[1]) + d2 * (p[3]-p[2]);
		c[0] = p[1]         ;
		c[1] = p[1] + m1 / 3;
		c[2] = p[2] - m2 / 3;
		c[3] = p[2]         ;
	}
}

#endif//_xtal_diagram_h_
