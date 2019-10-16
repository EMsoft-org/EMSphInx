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

#ifndef _SVG_HPP_
#define _SVG_HPP_

#include <vector>
#include <sstream>

namespace svg {
	struct Element;//an slement of an svg
	struct Group  ;//element groups

	class SVG {
		public:
			//@brief  : create an SVG
			//@param w: width  of SVG
			//@param h: height of SVG
			SVG(const double w, const double h);

			//@brief  : copy constructor
			//@param s: svg to copy
			SVG(const SVG& s);
			//@brief   : add an element to the svg
			//@param el: element to add
			void add(const Element& el);

			//@brief   : add a group to the svg
			//@param gr: element to add
			void add(const Group& gr);

			//@brief         : write the SVG to an ostream
			//@param fileName: name of file to write to
			void write(std::string fileName) const;

			//@brief    : write svg to an ostream
			//@brief os : ostream to write to
			//@brief svg: svg to write
			friend std::ostream& operator<<(std::ostream& os, const SVG& svg) {return os << svg.ss.str() << "</svg>\n";}

		private:
			double            width  ;//canvas width  in pixels
			double            height ;//canvas height in pixels
			double            view[4];//viewbox as minx, miny, width, height
			std::stringstream ss     ;//svg (without closing)
	};


	////////////////////////////////////////////////////////////////////////
	//               Elements, Groups, and Transformations                //
	////////////////////////////////////////////////////////////////////////

	struct Transform;//spatial transformations

	//an SVG is a collection of elements
	struct Element {
		//@brief   : convert the element to an XML string
		//@param os: ostream to write to
		virtual std::ostream& write(std::ostream& os) const = 0;

		//@brief   : write the element to an ostream
		//@brief os: ostream to write to
		//@brief el: element to write
		friend std::ostream& operator<<(std::ostream& os, const Element& el) {return el.write(os);}

		//@brief : get the string representation of the element
		//@return: XML string
		std::string to_string() const;

		//@brief    : transform the element
		//@param trs: transformation to apply
		//@return   : transformed element (in single component group)
		Group transform(const Transform trs) const;

		//@brief  : add an additional translation
		//@param x: x translation
		//@param y: y translation
		//@return : this
		Group translate(const double x, const double y) const;

		//@brief       : add an additional scaling
		//@param factor: scale factor (1 = no scaling)
		//@return      : this
		Group scale(const double factor) const;

		//@brief      : add an additional rotation
		//@param angle: rotation angle in degrees
		//@return     : this
		//@note       : positive angles are clockwise rotations
		Group rotate(const double angle) const;

		//@brief   : create a group with a second element
		//@param el: other element to add
		//@return  : group of {this, el}
		Group add(const Element& el) const;

		//the element class needs to be abstract
		virtual ~Element() = 0;
	};

	//a group of elements (potentially with a transformation applied)
	struct Group {
		//@brief   : add an element to the group
		//@param el: element to add
		//@return  : this
		Group& add(const Element& el);

		//@brief    : add an additional transformation
		//@param trs: transform to add
		//@return   : this
		Group& transform(const Transform trs);

		//@brief  : add an additional translation
		//@param x: x translation
		//@param y: y translation
		//@return : this
		Group& translate(const double x, const double y);

		//@brief       : add an additional scaling
		//@param factor: scale factor (1 = no scaling)
		//@return      : this
		Group& scale(const double factor);

		//@brief      : add an additional rotation
		//@param angle: rotation angle in degrees
		//@return     : this
		//@note       : positive angles are clockwise rotations
		Group& rotate(const double angle);

		//@brief   : convert the group to an XML string
		//@param os: ostream to write to
		std::ostream& write(std::ostream& os, std::string prefix = "\t") const;

		std::vector<std::string> elements  ;//string representation of elements in group
		std::vector<Transform  > transforms;//list of transformations to apply
	};

	//a spatial transformation
	//currently only translation, isotropic scaling, and rotation are supported
	struct Transform {
		//@brief  : set the rotation angle in degrees
		//@param x: x translation
		//@param y: y translation
		//@return : this
		Transform& translate(const double x, const double y);

		//@brief       : set the rotation angle in degrees
		//@param factor: scale factor (1 = no scaling)
		//@return      : this
		Transform& scale(const double factor);

		//@brief      : set the rotation angle in degrees
		//@param angle: rotation angle in degrees
		//@return     : this
		//@note       : positive angles are clockwise rotations
		Transform& rotate(const double angle);

		//@brief  : construct a transformation from a translation
		//@param x: x translation
		//@param y: y translation
		//@return : translation transform
		static Transform Translate(const double x, const double y);

		//@brief       : construct a transformation from a translation
		//@param factor: scale factor (1 = no scaling)
		//@return      : translation transform
		static Transform Scale(const double factor);

		//@brief      : construct a transformation from a translation
		//@param angle: rotation angle in degrees
		//@return     : translation transform
		//@note       : positive angles are clockwise rotations
		static Transform Rotate(const double angle);

		//@brief   : convert the element to an XML string
		//@param os: ostream to write to
		std::ostream& write(std::ostream& os) const;

		//@brief    : write the transform to an ostream
		//@brief os : ostream to write to
		//@brief trs: transform to write
		friend std::ostream& operator<<(std::ostream& os, const Transform& trs) {return trs.write(os);}

		private:
			double  tx ;//x translation
			double  ty ;//y translation
			double  scl;//scale
			double  rot;//rotation in degrees
			uint8_t flg;//transform flags
	};

	////////////////////////////////////////////////////////////////////////
	//                         Element Attributes                         //
	////////////////////////////////////////////////////////////////////////

	//@brief: class to represent an svg color
	struct Color {
		bool   none  ;//flag for the special 'none' color
		double rgb[3];//rgb values as 0->1

		//@brief: construct an empty color ('none')
		Color() : none(true) {}

		//@brief  : construct a color from an rgb triplet
		//@param r: fractional red   value [0,1]
		//@param g: fractional green value [0,1]
		//@param b: fractional blue  value [0,1]
		Color(const double r, const double g, const double b);

		//@brief : set the color to no color
		//@return: *this
		Color& setNone();

		//@brief  : set the color
		//@param r: fractional red   value [0,1]
		//@param g: fractional green value [0,1]
		//@param b: fractional blue  value [0,1]
		//@return : *this
		Color& setRGB(const double r, const double g, const double b);

		//@brief   : write color to an ostream
		//@param os: ostream to write to
		std::ostream& write(std::ostream& os) const;

		//@brief     : write color to an ostream
		//@brief os : ostream to write to
		//@brief clr: color to write
		friend std::ostream& operator<<(std::ostream& os, const Color& clr) {return clr.write(os);}
	};

	//@brief: class to represent an svg stroke
	struct Stroke {
		double width;//width of stroke
		Color  color;//color of stroke

		//@brief  : construct a stroke from a width + color
		//@param w: width of stroke
		//@param c: color of stroke
		//@note   : defaults to 1 wide black
		Stroke(const double w = 1, const Color c = Color(0, 0, 0)) : color(c), width(w) {}

		//@brief  : construct a stroke from a color only
		//@param c: color of stroke
		Stroke(const Color c) : Stroke(1, c) {}

		//@brief  : set stroke width
		//@param w: width to set
		//@return : this
		Stroke& setWidth(const double w);

		//@brief  : set stroke color
		//@param c: color to set
		//@return : this
		Stroke& setColor(const Color c);

		//@brief  : set the color
		//@param r: fractional red   value [0,1]
		//@param g: fractional green value [0,1]
		//@param b: fractional blue  value [0,1]
		//@return : *this
		Stroke& setColor(const double r, const double g, const double b) {return setColor(Color(r, g, b));}

		//@brief   : write stroke to an ostream
		//@param os: ostream to write to
		std::ostream& write(std::ostream& os) const {return os << "stroke=\"" << color << "\" stroke-width=\"" << width << "\" ";}

		//@brief    : write stroke to an ostream
		//@param os : ostream to write to
		//@param stk: stroke to write
		friend std::ostream& operator<<(std::ostream& os, const Stroke& stk) {return stk.write(os);}
	};

	//@brief: class to represent an svg fill
	struct Fill {
		Color color;//for now a fill is just a color (but it could be e.g. a gradient in the future)

		//@brief: default fill (none)
		Fill() {}

		//@brief  : construct fill from a color
		//@param c: color
		Fill(const Color c) : color(c) {}

		//@brief  : set fill color
		//@param c: color to set
		Fill& setColor(const Color c);

		//@brief  : set fill color
		//@param c: color to set
		Fill& setColor(const double r, const double g, const double b) {return setColor(Color(r, g, b));}

		//@brief  : set fill color
		//@param c: color to set
		Fill& setNone() {return setColor(Color());}

		//@brief   : write fill to an ostream
		//@param os: ostream to write to
		std::ostream& write(std::ostream& os) const {return os << "fill=\"" << color << "\" ";}

		//@brief   : write fill to an ostream
		//@param os: ostream to write to
		//@param fl: fill to write
		friend std::ostream& operator<<(std::ostream& os, const Fill& fl) {return fl.write(os);}
	};

	////////////////////////////////////////////////////////////////////////
	//                               Images                               //
	////////////////////////////////////////////////////////////////////////
	struct Image : public Element {
		enum PixelType {Gray = 1, GrayAlpha = 2, RGB = 3, RGBA = 4};//allowed pixel types

		double               x     ;//x coordinate of origin
		double               y     ;//y coordinate of origin
		size_t               width ;//width in pixels
		size_t               height;//height in pixels
		PixelType            type  ;//pixel type
		std::vector<uint8_t> buff  ;//raw data buffer

		//@brief  : construct an image
		//@param w: width  of image in pixels
		//@param h: height of image in pixels
		//@param p: pixel type
		Image(size_t w, size_t h, PixelType p = Gray) : x(0), y(0), width(w), height(h), type(p), buff(w * h * (size_t)p, 0x00) {}

		//@brief   : write complete image tag to an ostream
		//@param os: ostream to write to
		std::ostream& write(std::ostream& os) const;
	};

	////////////////////////////////////////////////////////////////////////
	//                            Basic Shapes                            //
	////////////////////////////////////////////////////////////////////////

	struct Shape : public Element {
		Fill   fill  ;
		Stroke stroke;

		//@brief   : write complete shape tag to an ostream
		//@param os: ostream to write to
		std::ostream& write(std::ostream& os) const;

		//@brief   : write shape tag contents to an ostream
		//@param os: ostream to write to
		virtual std::ostream& writeElement(std::ostream& os) const = 0;

		//@brief  : set the shape's fill
		//@param f: fill to set
		//@return : this
		Shape& setFill(const Fill& f);

		//@brief  : set the shape's fill
		//@param r: red [0, 1]
		//@param g: red [0, 1]
		//@param b: red [0, 1]
		//@return : this
		Shape& setFill(const double r, const double g, const double b) {return setFill(Color(r, g, b));}

		//@brief  : remove the shape's fill
		//@return : this
		Shape& removeFill() {return setFill(Fill(Color()));}

		//@brief  : set the shape's stroke
		//@param s: stroke to set
		//@return : this
		Shape& setStroke(const Stroke& s);

		//@brief  : set the shape's stroke
		//@param w: width of stroke
		//@param c: color of stroke
		//@return : this
		Shape& setStroke(const double w, const Color c) {return setStroke(Stroke(w, c));}

		//@brief  : set the shape's stroke width
		//@param w: width of stroke
		//@return : this
		Shape& setStrokeWidth(const double w) {return setStroke(w, stroke.color);}

		//@brief  : set the shape's stroke color
		//@param c: color of stroke
		//@return : this
		Shape& setStrokeColor(const Color c) {return setStroke(stroke.width, c);}

		//@brief  : set the shape's stroke color
		//@param c: color of stroke
		//@return : this
		Shape& setStrokeColor(const double r, const double g, const double b) {return setStroke(stroke.width, Color(r, g, b));}

		//@brief  : remove the stroke
		//@return : this
		Shape& removeStroke() {return setStrokeColor(Color());}
	};

	struct Path : public Shape {
		//@brief: construct an empty path
		Path() {}

		//@preif  : construct a path starting at specified coordinates
		//@param x: x coordinate of start
		//@param y: y coordinate of start
		Path(const double x, const double y) {moveTo(x, y, true);}

		//@brief    : start a new sub path at the given coordinates (no stroke in between)
		//@param x  : x coordinate of new sub path start
		//@param y  : y coordinate of new sub path start
		//@param abs: true/false if x/y are absolute/relative (to previous position)
		//@param abs: true/false if x/y are absolute/relative (to previous position)
		//@note     : 'M' / 'm'
		void moveTo(const double x, const double y, const bool abs);

		//@brief: close the current subpath
		//@note : 'Z' / 'z'
		void close();

		//@brief    : draw a straight line from the current position to the given position
		//@param x  : x coordinate to draw line to
		//@param y  : y coordinate to draw line to
		//@param abs: true/false if x/y are absolute/relative (to previous position)
		//@note     : 'L' / 'l'
		void lineTo(const double x, const double y, const bool abs);

		//@brief    : draw a horizontal line from the current position to the given position
		//@param x  : x coordinate to draw line to
		//@param abs: true/false if x/y are absolute/relative (to previous position)
		//@note     : 'H' / 'h'
		void hLineTo(const double x, const bool abs);

		//@brief    : draw a vertical line from the current position to the given position
		//@param y  : y coordinate to draw line to
		//@param abs: true/false if x/y are absolute/relative (to previous position)
		//@note     : 'V' / 'v'
		void vLineTo(const double y, const bool abs);

		//@brief    : draw a cubic bezier curve
		//@param x1 : x coordinate of first control point (for start point)
		//@param y1 : y coordinate of first control point (for start point)
		//@param x2 : x coordinate of second control point (for end point)
		//@param y2 : y coordinate of second control point (for end point)
		//@param x  : x coordinate of end point
		//@param y  : y coordinate of end point
		//@param abs: true/false if x/y are absolute/relative (to previous position)
		//@note     : 'C' / 'c'
		void curveTo(const double x1, const double y1, const double x2, const double y2, const double x, const double y, const bool abs);

		//@brief    : draw a cubic bezier curve
		//@param x2 : x coordinate of second control point (for end point)
		//@param y2 : y coordinate of second control point (for end point)
		//@param x  : x coordinate of end point
		//@param y  : y coordinate of end point
		//@param abs: true/false if x/y are absolute/relative (to previous position)
		//@note     : x1/y1 are reflection of the second control point on the previous command relative to the current control point
		//@note     : 'S' / 's'
		void smoothCurveTo(const double x2, const double y2, const double x, const double y, const bool abs);

		//@brief    : draw a quadratic bezier curve
		//@param x1 : x coordinate of control point
		//@param y1 : y coordinate of control point
		//@param x  : x coordinate of end point
		//@param y  : y coordinate of end point
		//@param abs: true/false if x/y are absolute/relative (to previous position)
		//@note     : 'Q' / 'q'
		void quadTo(const double x1, const double y1, const double x, const double y, const bool abs);

		//@brief    : draw a quadratic bezier curve
		//@param x  : x coordinate of end point
		//@param y  : y coordinate of end point
		//@param abs: true/false if x/y are absolute/relative (to previous position)
		//@note     : x1/y1 are reflection of the second control point on the previous command relative to the current control point
		//@note     : 'T' / 't'
		void smoothQuadTo(const double x, const double y, const bool abs);

		//@brief    : draw an elliptical arc
		//@param rx : x axis of ellipse
		//@param ry : y axis of ellipse
		//@param rot: rotation of axis in degrees (relative to current coordinate system)
		//@param lrg: large arc flag
		//@param swp: sweep flag
		//@param x  : x coordinate of end point
		//@param y  : y coordinate of end point
		//@param abs: true/false if x/y are absolute/relative (to previous position)
		//@note     : 'A' / 'a'
		void ellipticalArc(const double rx, const double ry, const double rot, const int lrg, const int swp, const double x, const double y, const bool abs);

		//@brief   : write path to an ostream
		//@param os: ostream to write to
		std::ostream& writeElement(std::ostream& os) const;

		private:
			struct ControlPoint {
				char   c;
				double params[6];
			};

			std::vector<ControlPoint> pts;
	};

	//rectangle
	struct Rect : public Shape {
		double x     ;//x origin
		double y     ;//y origin
		double width ;//width  of rectangle
		double height;//height of rectangle
		double rx    ;//x radius of ellipse for corner round of
		double ry    ;//y radius of ellipse for corner round of

		//@brief  : construct a rectangle
		//@param x: x origin
		//@param y: y origin
		//@param w: width
		//@param h: height
		Rect(const double x, const double y, const double w, const double h) : x(x), y(y), width(w), height(h), rx(0), ry(0) {}

		//@brief   : write rect to an ostream
		//@param os: ostream to write to
		std::ostream& writeElement(std::ostream& os) const;
	};

	//circle
	struct Circle : public Shape {
		double cx;//x center
		double cy;//y center
		double r ;//radius

		//@brief   : construct a circle
		//@param r : radius
		//@param cx: x center
		//@param cy: y center
		Circle(const double r, const double cx = 0, const double cy = 0) : cx(cx), cy(cy), r(r) {}

		//@brief   : write circle to an ostream
		//@param os: ostream to write to
		std::ostream& writeElement(std::ostream& os) const;
	};

	//ellipse
	struct Ellipse : public Shape {
		double cx;//x center
		double cy;//y center
		double rx;//x radius
		double ry;//y radius

		//@brief   : construct an ellipse
		//@param rx: x radius
		//@param ry: y radius
		//@param cx: x center
		//@param cy: y center
		Ellipse(const double rx, const double ry, const double cx = 0, const double cy = 0) : cx(cx), cy(cy), rx(rx), ry(ry) {}

		//@brief   : write ellipse to an ostream
		//@param os: ostream to write to
		std::ostream& writeElement(std::ostream& os) const;
	};

	//line
	struct Line : public Shape {
		double x1;//x start coordinate
		double y1;//y start coordinate
		double x2;//x end coordinate
		double y2;//y end coordinate

		//@brief   : construct a line
		//@param x1: x start coordinate
		//@param y1: y start coordinate
		//@param x2: x end coordinate
		//@param y2: y end coordinate
		Line(const double x1, const double y1, const double x2, const double y2) : x1(x1), y1(y1), x2(x2), y2(y2) {}

		//@brief   : write line to an ostream
		//@param os: ostream to write to
		std::ostream& writeElement(std::ostream& os) const;
	};

	//polyline
	struct Polyline : public Shape {
		std::vector< std::pair<double, double> > pts;//xy coordinates of points

		//@brief  : add a new point to the line
		//@param x: x coordinate of point to add
		//@param y: y coordinate of point to add
		void add(const double x, const double y) {pts.emplace_back(x, y);}

		//@brief   : write polyline to an ostream
		//@param os: ostream to write to
		std::ostream& writeElement(std::ostream& os) const;

		//@brief: get the name of this object
		virtual std::string name() const {return "polyline";}
	};

	//polygon (just a closed polyline)
	struct Polygon : public Polyline {
		//@brief: get the name of this object
		std::string name() const {return "polygon";}

		//@brief  : construct a regular polygon
		//@param n: number of sides (must be at least 3)
		//@param r: distance from origin to points
		//@param x: x origin
		//@param y: y origin
		static Polygon Regular(const size_t n, const double r, const double x = 0, const double y = 0);
	};


	////////////////////////////////////////////////////////////////////////
	//                                Text                                //
	////////////////////////////////////////////////////////////////////////

	struct Text : public Element {
		enum class Style      {Normal, Italic, Oblique};
		enum class Weight     {Normal, Bold, Bolder, Lighter};
		enum class Decoration {None = 0, Underline = 1, Overline = 2, Throughline = 4, Blink = 8};//bitmask

		double               x     ;//x coordinate of origin
		double               y     ;//y coordinate of origin
		std::string          font  ;//font family
		double               size  ;//font size
		Style                style ;//font style
		Weight               weight;//font weight
		Decoration           decor ;//font decoration
		std::string          text  ;//string

		//@brief  : construct text
		Text() : x(0), y(0), font("Helvetica"), size(12), style(Style::Normal), weight(Weight::Normal), decor(Decoration::None), text() {}

		//@brief   : write text tag to an ostream
		//@param os: ostream to write to
		std::ostream& write(std::ostream& os) const;
	};
}

////////////////////////////////////////////////////////////////////////
//                       Implementation Details                       //
////////////////////////////////////////////////////////////////////////

#include <fstream>
#include <cmath>
#include <stdexcept>

#include "util/base64.hpp"

#define MINIZ_NO_STDIO
#define MINIZ_NO_TIME
#define MINIZ_NO_ZLIB_APIS
#include "miniz/miniz.c"

#include "constants.hpp"

namespace svg {

	////////////////////////////////////////////////////////////////////////
	//                           Element Members                           //
	////////////////////////////////////////////////////////////////////////

	//@brief : get the string representation of the element
	//@return: XML string
	std::string Element::to_string() const {
		std::stringstream ss;
		write(ss);
		return ss.str();
	}

	//@brief    : transform the element
	//@param trs: transformation to apply
	//@return   : transformed element (in single component group)
	Group Element::transform(const Transform trs) const {
		Group grp;
		return grp.add(*this).transform(trs);
	}

	//@brief  : add an additional translation
	//@param x: x translation
	//@param y: y translation
	//@return : this
	Group Element::translate(const double x, const double y) const {
		return transform(Transform::Translate(x, y));
	}

	//@brief       : add an additional scaling
	//@param factor: scale factor (1 = no scaling)
	//@return      : this
	Group Element::scale(const double factor) const {
		return transform(Transform::Scale(factor));
	}

	//@brief      : add an additional rotation
	//@param angle: rotation angle in degrees
	//@return     : this
	//@note       : positive angles are clockwise rotations
	Group Element::rotate(const double angle) const {
		return transform(Transform::Rotate(angle));
	}

	//@brief   : create a group with a second element
	//@param el: other element to add
	//@return  : group of {this, el}
	Group Element::add(const Element& el) const {
		Group grp;
		return grp.add(*this).add(el);
	}

	Element::~Element() {}

	////////////////////////////////////////////////////////////////////////
	//                           Group Members                            //
	////////////////////////////////////////////////////////////////////////

	//@brief   : add an element to the group
	//@param el: element to add
	//@return  : this
	Group& Group::add(const Element& el) {
		elements.push_back(el.to_string());
		return *this;
	}

	//@brief    : add an additional transformation
	//@param trs: transform to add
	//@return   : this
	Group& Group::transform(const Transform trs) {
		transforms.push_back(trs);
		return *this;
	}

	//@brief  : add an additional translation
	//@param x: x translation
	//@param y: y translation
	//@return : this
	Group& Group::translate(const double x, const double y) {
		return transform(Transform::Translate(x, y));
	}

	//@brief       : add an additional scaling
	//@param factor: scale factor (1 = no scaling)
	//@return      : this
	Group& Group::scale(const double factor) {
		return transform(Transform::Scale(factor));
	}

	//@brief      : add an additional rotation
	//@param angle: rotation angle in degrees
	//@return     : this
	//@note       : positive angles are clockwise rotations
	Group& Group::rotate(const double angle) {
		return transform(Transform::Rotate(angle));
	}

	//@brief   : convert the group to an XML string
	//@param os: ostream to write to
	std::ostream& Group::write(std::ostream& os, std::string prefix) const {
		if(transforms.empty()) {
			//open group
			os << prefix << "<g>\n";

			//loop over elements writing
			for(const std::string& el : elements) os << prefix << '\t' << el << '\n';

			//close group
			os << prefix << "</g>\n";
		} else {
			//loop over transforms opening a group
			for(size_t i = 0; i < transforms.size(); i++) {
				os << prefix;
				for(size_t j = 0; j < i; j++) os << '\t';
				os << "<g " << transforms[i] << ">\n";//open group
			}

			//loop over elements writing
			for(const std::string& el : elements) {
				os << prefix;
				for(size_t j = 1; j <= transforms.size(); j++) os << '\t';
				os << el << '\n';
			}

			//loop over transforms closing groups
			for(size_t i = transforms.size() - 1; i < transforms.size(); i--) {
				os << prefix;
				for(size_t j = 0; j < i; j++) os << '\t';
				os << "</g>\n";//close group
			}
		}
		return os;
	}

	////////////////////////////////////////////////////////////////////////
	//                         Transform Members                          //
	////////////////////////////////////////////////////////////////////////

	//@brief  : set the rotation angle in degrees
	//@param x: x translation
	//@param y: y translation
	//@return : this
	Transform& Transform::translate(const double x, const double y) {
		tx = x;
		ty = y;
		flg |= 0x01;
		return *this;
	}

	//@brief       : set the rotation angle in degrees
	//@param factor: scale factor (1 = no scaling)
	//@return      : this
	Transform& Transform::scale(const double factor) {
		scl = factor;
		flg |= 0x02;
		return *this;
	}

	//@brief      : set the rotation angle in degrees
	//@param angle: rotation angle in degrees
	//@return     : this
	//@note       : positive angles are clockwise rotations
	Transform& Transform::rotate(const double angle) {
		rot = angle;
		flg |= 0x04;
		return *this;
	}

	//@brief  : construct a transformation from a translation
	//@param x: x translation
	//@param y: y translation
	//@return : translation transform
	Transform Transform::Translate(const double x, const double y) {
		return Transform().translate(x, y);
	}

	//@brief       : construct a transformation from a translation
	//@param factor: scale factor (1 = no scaling)
	//@return      : translation transform
	Transform Transform::Scale(const double factor) {
		return Transform().scale(factor);
	}

	//@brief      : construct a transformation from a translation
	//@param angle: rotation angle in degrees
	//@return     : translation transform
	//@note       : positive angles are clockwise rotations
	Transform Transform::Rotate(const double angle) {
		return Transform().rotate(angle);
	}

	//@brief   : convert the element to an XML string
	//@param os: ostream to write to
	std::ostream& Transform::write(std::ostream& os) const {
		if(flg > 0) {//only both if we have a transform
			os << "transform=\" ";
			if(flg & 0x01) os << "translate(" << tx << ',' << ty << ") ";
			if(flg & 0x02) os << "scale(" << scl << ") ";
			if(flg & 0x04) os << "rotate(" << rot << ") ";
			os << "\" ";
		}
		return os;
	}

	////////////////////////////////////////////////////////////////////////
	//                            SVG members                             //
	////////////////////////////////////////////////////////////////////////

	//@brief  : create an SVG
	//@param w: width  of SVG in pixels
	//@param h: height of SVG in pixels
	SVG::SVG(const double w, const double h) : width(w), height(h) {
		view[0] = 0;
		view[1] = 0;
		view[2] = w;
		view[3] = h;

		ss << "<?xml version=\"1.0\" encoding=\"UTF-8\" ?>\n";
		ss << "<svg width=\"" << width << "\" height=\"" << height
		   << "\" viewBox=\"" << view[0] << ' ' << view[1] << ' ' << view[2] << ' ' << view[3]
		   << "\" xmlns=\"http://www.w3.org/2000/svg\" "
		   << "xmlns:xlink=\"http://www.w3.org/1999/xlink\">\n";
	}

	//@brief  : copy constructor
	//@param s: svg to copy
	SVG::SVG(const SVG& s) {
		width  = s.width ;
		height = s.height;
		std::copy(s.view, s.view + 4, view);
		ss << s.ss.str();
	}

	//@brief   : add an element to the svg
	//@param el: element to add
	void SVG::add(const Element& el) {
		ss << '\t' << el << '\n';
	}

	//@brief   : add a group to the svg
	//@param gr: element to add
	void SVG::add(const Group& gr) {
		gr.write(ss, "\t");
	}

	//@brief         : write the SVG to an ostream
	//@param fileName: name of file to write to
	void SVG::write(std::string fileName) const {
		std::ofstream os(fileName);
		os << *this;
	}

	////////////////////////////////////////////////////////////////////////
	//                         Element Attributes                         //
	////////////////////////////////////////////////////////////////////////

	//@brief  : construct a color from an rgb triplet
	//@param r: fractional red   value [0,1]
	//@param g: fractional green value [0,1]
	//@param b: fractional blue  value [0,1]
	Color::Color(const double r, const double g, const double b) {
		none = false;
		rgb[0] = r;
		rgb[1] = g;
		rgb[2] = b;
	}

	//@brief : set the color to no color
	//@return: *this
	Color& Color::setNone() {
		none = true;
		return *this;
	}

	//@brief  : set the color
	//@param r: fractional red   value [0,1]
	//@param g: fractional green value [0,1]
	//@param b: fractional blue  value [0,1]
	//@return : *this
	Color& Color::setRGB(const double r, const double g, const double b) {
		rgb[0] = r;
		rgb[1] = g;
		rgb[2] = b;
		return *this;
	}

	//@brief   : write color to an ostream
	//@param os: ostream to write to
	std::ostream& Color::write(std::ostream& os) const {
		if(none) {
			os << "none";
		} else {
			os << "rgb(" << rgb[0]*255 << ',' << rgb[1]*255 << ',' << rgb[2]*255 << ')';
		}
		return os;
	}

	//@brief  : set stroke width
	//@param w: width to set
	//@return : this
	Stroke& Stroke::setWidth(const double w) {
		width = w;
		return *this;
	}

	//@brief  : set stroke color
	//@param c: color to set
	//@return : this
	Stroke& Stroke::setColor(const Color c) {
		color = c;
		return *this;
	}

	//@brief  : set fill color
	//@param c: color to set
	Fill& Fill::setColor(const Color c) {
		color = c;
		return *this;
	}

	////////////////////////////////////////////////////////////////////////
	//                               Images                               //
	////////////////////////////////////////////////////////////////////////

	//@brief   : write complete image tag to an ostream
	//@param os: ostream to write to
	std::ostream& Image::write(std::ostream& os) const {
		//sanity check
		const int chan = (int) type;//are there 1, 2, 3 or 4 samples per pixel
		if(buff.size() != width * height * chan) throw std::logic_error("image shape doesn't match buffer size");

		//convert to png in memory
		size_t pngSize = 0;
		const mz_uint compr = MZ_BEST_COMPRESSION;//compression level [0,10]
		const mz_bool flip  = MZ_FALSE;//flip the image?
		void *pPNG_data = tdefl_write_image_to_png_file_in_memory_ex((void*)buff.data(), (int)width, (int)height, chan, &pngSize, compr, flip);
		if(!pPNG_data) throw std::runtime_error("failed to create PNG image");

		//base64 encode
		os << "<image x = \"" << x << "\" y = \"" << y
		   << "\" width=\"" << width << "\" height=\"" << height << "\" "
		   << "xlink:href=\"data:image/png;base64,";
		base64::encode((char*)pPNG_data, pngSize, os);//write png using base64 encoding
		mz_free(pPNG_data);//cleanup memory allocated by png creation
		return os << "\"/>";
	}

	////////////////////////////////////////////////////////////////////////
	//                            Basic Shapes                            //
	////////////////////////////////////////////////////////////////////////

	//@brief   : write complete shape tag to an ostream
	//@param os: ostream to write to
	std::ostream& Shape::write(std::ostream& os) const {
		os << "<";
		writeElement(os);
		return os << fill << stroke << "/>"; 
	}

	//@brief  : set the shape's fill
	//@param f: fill to set
	//@return : this
	Shape& Shape::setFill(const Fill& f) {
		fill = f;
		return *this;
	}

	//@brief  : set the shape's stroke
	//@param s: stroke to set
	//@return : this
	Shape& Shape::setStroke(const Stroke& s) {
		stroke = s;
		return *this;
	}

	//@brief    : start a new sub path at the given coordinates (no stroke in between)
	//@param x  : x coordinate of new sub path start
	//@param y  : y coordinate of new sub path start
	//@param abs: true/false if x/y are absolute/relative (to previous position)
	//@param abs: true/false if x/y are absolute/relative (to previous position)
	//@note     : 'M' / 'm'
	void Path::moveTo(const double x, const double y, const bool abs) {
		ControlPoint p;
		p.c = abs ? 'M' : 'm';
		p.params[0] = x;
		p.params[1] = y;
		pts.push_back(p);
	}

	//@brief: close the current subpath
	//@note : 'Z' / 'z'
	void Path::close() {
		ControlPoint p;
		p.c = 'z';//case doesn't matter for z
		pts.push_back(p);
	}

	//@brief    : draw a straight line from the current position to the given position
	//@param x  : x coordinate to draw line to
	//@param y  : y coordinate to draw line to
	//@param abs: true/false if x/y are absolute/relative (to previous position)
	//@note     : 'L' / 'l'
	void Path::lineTo(const double x, const double y, const bool abs) {
		ControlPoint p;
		p.c = abs ? 'L' : 'l';
		p.params[0] = x;
		p.params[1] = y;
		pts.push_back(p);
	}

	//@brief    : draw a horizontal line from the current position to the given position
	//@param x  : x coordinate to draw line to
	//@param abs: true/false if x/y are absolute/relative (to previous position)
	//@note     : 'H' / 'h'
	void Path::hLineTo(const double x, const bool abs) {
		ControlPoint p;
		p.c = abs ? 'H' : 'h';
		p.params[0] = x;
		pts.push_back(p);
	}

	//@brief    : draw a vertical line from the current position to the given position
	//@param y  : y coordinate to draw line to
	//@param abs: true/false if x/y are absolute/relative (to previous position)
	//@note     : 'V' / 'v'
	void Path::vLineTo(const double y, const bool abs) {
		ControlPoint p;
		p.c = abs ? 'V' : 'v';
		p.params[0] = y;
		pts.push_back(p);
	}

	//@brief    : draw a cubic bezier curve
	//@param x1 : x coordinate of first control point (for start point)
	//@param y1 : y coordinate of first control point (for start point)
	//@param x2 : x coordinate of second control point (for end point)
	//@param y2 : y coordinate of second control point (for end point)
	//@param x  : x coordinate of end point
	//@param y  : y coordinate of end point
	//@param abs: true/false if x/y are absolute/relative (to previous position)
	//@note     : 'C' / 'c'
	void Path::curveTo(const double x1, const double y1, const double x2, const double y2, const double x, const double y, const bool abs) {
		ControlPoint p;
		p.c = abs ? 'C' : 'c';
		p.params[0] = x1;
		p.params[1] = y1;
		p.params[2] = x2;
		p.params[3] = y2;
		p.params[4] = x ;
		p.params[5] = y ;
		pts.push_back(p);
	}

	//@brief    : draw a cubic bezier curve
	//@param x2 : x coordinate of second control point (for end point)
	//@param y2 : y coordinate of second control point (for end point)
	//@param x  : x coordinate of end point
	//@param y  : y coordinate of end point
	//@param abs: true/false if x/y are absolute/relative (to previous position)
	//@note     : x1/y1 are reflection of the second control point on the previous command relative to the current control point
	//@note     : 'S' / 's'
	void Path::smoothCurveTo(const double x2, const double y2, const double x, const double y, const bool abs) {
		ControlPoint p;
		p.c = abs ? 'S' : 's';
		p.params[0] = x2;
		p.params[1] = y2;
		p.params[2] = x ;
		p.params[3] = y ;
		pts.push_back(p);
	}

	//@brief    : draw a quadratic bezier curve
	//@param x1 : x coordinate of control point
	//@param y1 : y coordinate of control point
	//@param x  : x coordinate of end point
	//@param y  : y coordinate of end point
	//@param abs: true/false if x/y are absolute/relative (to previous position)
	//@note     : 'Q' / 'q'
	void Path::quadTo(const double x1, const double y1, const double x, const double y, const bool abs) {
		ControlPoint p;
		p.c = abs ? 'Q' : 'q';
		p.params[0] = x1;
		p.params[1] = y1;
		p.params[2] = x ;
		p.params[3] = y ;
		pts.push_back(p);
	}

	//@brief    : draw a quadratic bezier curve
	//@param x  : x coordinate of end point
	//@param y  : y coordinate of end point
	//@param abs: true/false if x/y are absolute/relative (to previous position)
	//@note     : x1/y1 are reflection of the second control point on the previous command relative to the current control point
	//@note     : 'T' / 't'
	void Path::smoothQuadTo(const double x, const double y, const bool abs) {
		ControlPoint p;
		p.c = abs ? 'T' : 't';
		p.params[0] = x;
		p.params[1] = y;
		pts.push_back(p);
	}

	//@brief    : draw an elliptical arc
	//@param rx : x axis of ellipse
	//@param ry : y axis of ellipse
	//@param rot: rotation of axis in degrees (relative to current coordinate system)
	//@param lrg: large arc flag
	//@param swp: sweep flag
	//@param x  : x coordinate of end point
	//@param y  : y coordinate of end point
	//@param abs: true/false if x/y are absolute/relative (to previous position)
	//@note     : 'A' / 'a'
	void Path::ellipticalArc(const double rx, const double ry, const double rot, const int lrg, const int swp, const double x, const double y, const bool abs) {
		if(!((lrg == 0 || lrg == 1) && (swp == 0 || swp == 1))) throw std::runtime_error("elliptical arc flags must be 0 or 1");
		ControlPoint p;
		p.c = abs ? 'A' : 'a';
		p.params[0] = rx ;
		p.params[1] = ry ;
		p.params[2] = rot;
		p.params[4] = x  ;
		p.params[5] = y  ;
		((int*)&p.params[3])[0] = lrg * 0x01 + swp * 0x02;
		pts.push_back(p);
	}

	//@brief   : write rect to an ostream
	//@param os: ostream to write to
	std::ostream& Path::writeElement(std::ostream& os) const {
		os << "path d=\"";
		for(const ControlPoint& p : pts) {
			os << p.c;
			switch(p.c) {
				// 0 parameter commands
				case 'Z':
				case 'z':
					break;

				// 1 parameter commands
				case 'H':
				case 'h':
				case 'V':
				case 'v':
					os << p.params[0] << ' ';
					break;

				// 2 parameter commands
				case 'M':
				case 'm':
				case 'L':
				case 'l':
				case 'T':
				case 't':
					os << p.params[0] << ' ' << p.params[1] << ' ';
					break;

				// 4 parameter commands
				case 'S':
				case 's':
				case 'Q':
				case 'q':
					os << p.params[0] << ' ' << p.params[1] << ' ';
					os << p.params[2] << ' ' << p.params[3] << ' ';
					break;

				// 6 parameter commands
				case 'C':
				case 'c':
					os << p.params[0] << ' ' << p.params[1] << ' ' << p.params[2] << ' ' << p.params[3] << ' ' << p.params[4] << ' ' << p.params[5] << ' ';
					break;

				// special cases
				case 'A':
				case 'a':
					os << p.params[0] << ' ' << p.params[1] << ' ' << p.params[2] << ' ';
					switch(((int*)&p.params[3])[0]) {
						case 0 : os << "0 0 "; break;
						case 1 : os << "1 0 "; break;
						case 2 : os << "0 1 "; break;
						case 3 : os << "1 1 "; break;
						default: throw std::logic_error("bad ellipse flag");
					}
					os << p.params[4] << ' ' << p.params[5] << ' ';
					break;

				default: throw std::logic_error("bad path flag");
			}
		}
		return os << "\" ";
	}

	//@brief   : write rect to an ostream
	//@param os: ostream to write to
	std::ostream& Rect::writeElement(std::ostream& os) const {
		return os << "rect "
		          << "x=\""      << x     << "\" y=\""      << y      << "\" "
		          << "width=\""  << width << "\" height=\"" << height << "\" "
		          << "rx=\""     << rx    << "\" ry=\""     << ry     << "\" ";
	}

	//@brief   : write circle to an ostream
	//@param os: ostream to write to
	std::ostream& Circle::writeElement(std::ostream& os) const {
		return os << "circle cx=\"" << cx << "\" cy=\"" << cy << "\" " << "r=\"" << r << "\" ";
	}

	//@brief   : write ellipse to an ostream
	//@param os: ostream to write to
	std::ostream& Ellipse::writeElement(std::ostream& os) const {
		return os << "ellipse cx=\"" << cx << "\" cy=\"" << cy << "\" rx=\"" << rx << "\" ry=\"" << ry << "\" ";
	}

	//@brief   : write line to an ostream
	//@param os: ostream to write to
	std::ostream& Line::writeElement(std::ostream& os) const {
		return os << "line x1=\"" << x1 << "\" y1=\"" << y1 << "\" x2=\"" << x2 << "\" y2=\"" << y2 << "\" ";
	}

	//@brief   : write poly to an ostream
	//@param os: ostream to write to
	std::ostream& Polyline::writeElement(std::ostream& os) const {
		os << name() << " points=\"";
		for(const std::pair<double, double> pt : pts) os << pt.first << ',' << pt.second << ' ';
		return os << "\" ";
	}

	//@brief  : construct a regular polygon
	//@param n: number of sides (must be at least 3)
	//@param r: distance from origin to points
	//@param x: x origin
	//@param y: y origin
	Polygon Polygon::Regular(const size_t n, const double r, const double x, const double y) {
		Polygon ply;
		for(size_t i = 0; i < n; i++) {
			const double theta = emsphinx::Constants<double>::pi2 * i / n + emsphinx::Constants<double>::pi;
			ply.add(std::sin(theta) * r + x, std::cos(theta) * r + y);
		}
		return ply;
	}


	////////////////////////////////////////////////////////////////////////
	//                                Text                                //
	////////////////////////////////////////////////////////////////////////
	
	//@brief   : write text tag to an ostream
	//@param os: ostream to write to
	std::ostream& Text::write(std::ostream& os) const {
		os << "<text x=\"" << x << "\" y=\"" << y << "\" ";
		os << "font-family=\"" << font << "\" ";
		os << "font-size=\"" << size << "\" ";
		os << "font-style=\"";
		switch(style) {
			case Style::Normal : os << "normal" ; break;
			case Style::Italic : os << "italic" ; break;
			case Style::Oblique: os << "oblique"; break;
		}
		os << "\" font-weight=\"";
		switch(weight) {
			case Weight::Normal : os << "normal" ; break;
			case Weight::Bold   : os << "bold"   ; break;
			case Weight::Bolder : os << "bolder" ; break;
			case Weight::Lighter: os << "lighter"; break;
		}
		os << "\" text-decoration=\"";
		if((int)decor & (int)Decoration::Underline  ) os << "underline "   ;
		if((int)decor & (int)Decoration::Overline   ) os << "overline "    ;
		if((int)decor & (int)Decoration::Throughline) os << "line-through ";
		if((int)decor & (int)Decoration::Blink      ) os << "blink "       ;
		return os << "\" >" << text << "</text>";
	}
}

#endif//_SVG_HPP_
