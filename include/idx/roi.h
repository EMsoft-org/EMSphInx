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

#ifndef _ROI_H_
#define _ROI_H_

#include <vector>
#include <string>

namespace emsphinx {

	enum class DrawMode {
		Rectangle,
		Ellipse  ,
		Polygon  ,
	};

	class RoiSelection {
		DrawMode                           mode  = DrawMode::Rectangle;//type of ROI selection
		bool                               mkSq  = false              ;//should the aspect ratio be fixed (rectangle / ellipse only)
		std::vector< std::pair<int, int> > pts                        ;//selection coordinates
		int                                curPt = -1                 ;//point currently being manipulated (or -1 for none)
		bool                               inv   = false              ;//is the ROI selection inverted (exclude selected instead of include)

		public:
			//drawing mode independent

			//@brief   : check if coordinates are near an existing point
			//@param x : coordinates to check
			//@param y : coordinates to check
			//@param r2: square radius cutoff
			//@return  : index of nearby point (-1 if none are nearby)
			int nearPt(const int x, const int y, const int r2) const;

			//@brief: determine if a point is currently selected
			bool hasSelection() const {return -1 == curPt;}

			//@brief   : try to select a point
			//@param x : coordinates to try selecting near
			//@param y : coordinates to try selecting near
			//@param r2: square radius cutoff
			//@return  : true if a point was selected
			bool trySelect(const int x, const int y, const int r2) {curPt = nearPt(x, y, r2); return -1 != curPt;}

			//@brief: clear current shape
			//@return: true if anything changed, false otherwise
			bool clear();

			//@brief : get the drawing mode
			//@return: drawing mode
			DrawMode getMode() const {return mode;}

			//@brief : chamge the drawing mode
			//@return: true if the mode changed, false otherwise
			bool changeMode(const DrawMode dm);

			//@brief  : set fixed aspect ratio
			//@param b: true to fix aspect ratio, false to free
			void setFixedAspect(const bool b) {mkSq = b;}

			//@brief : get current coordinate list
			//@return: read only access to list
			const std::vector< std::pair<int, int> >& getPts() const {return pts;}

			//@brief : check if there is a complete shape
			//@return: true if there is a complete shape, false otherwise
			const bool hasShape() const {return pts.empty() ? false : DrawMode::Polygon == mode ? (pts.front() == pts.back() && pts.size() > 2u && curPt + 1u < pts.size()) : true;}

			//@brief   : translate
			//@param dX: x shift
			//@param dY: y shift
			void translatePoints(const int dX, const int dY) {for(std::pair<int, int>& pt : pts) {pt.first += dX, pt.second += dY;}}

			//@brief    : update coordinates of a point
			//@param idx: point to update
			//@param x  : new x coordiante
			//@param y  : new y coordiante
			void movePoint(const size_t idx, const int x, const int y);

			//rectangle building

			//@brief  : start drawing a new rectangle
			//@param x: x origin of new rectangle
			//@param y: y origin of new rectangle
			void startRect(const int x, const int y);

			//@brief    : start drawing a new rectangle
			//@param ori: origin of new rectangle
			void finishRect() {curPt = -1;}

			//@brief  : resize a rectangle
			//@param x: new x coordinate
			//@param y: new y coordinate
			//@return : true if anything changed, false otherwise
			bool resizeRect(const int x, const int y);

			//@brief    : stop resizing a rectangle
			//@param pt : new coordinate
			void ungrabRect() {curPt = -1; if(pts.empty() ? true : pts.front() == pts.back()) clear();}

			//polygon building

			//@brief  : start drawing a polygon
			//@param x: x origin of polygon
			//@param y: y origin of polygon
			void startPoly(const int x, const int y);

			//@brief: check if we are currently building a polygon
			bool buildingPoly() const {return pts.empty() ? false : ( 2 == pts.size() ? true : pts.back() != pts.front() );}

			//@brief   : add a new point to an existing polygon
			//@param x : x coordinate of new point to add
			//@param y : y coordinate of new point to add
			//@param r2: cutoff radius^2 for 2 points being the same
			//@note    : finishes the polygon if pt is sufficiently close to the first point
			//@return  : true if polygon was closed, false otherwise
			bool addPolyPt(const int x, const int y, const int r2);

			//@brief : remove currently selected point from an existing polygon
			//@return: true if point was removed, false otherwise
			bool remPolyPt();

			//@brief      : duplicated currently selected point in an existing polygon
			//@param delta: offset from previous point
			//@return     : true if a point was duplicate, false otherwise
			bool dupPolyPt(const int delta);

			//@brief   : finish drawing a polygon
			//@return  : true if polygon was closed, false otherwise
			bool finishPoly();

			//@brief   : move a polygon vertex
			//@param x : new x coordinate
			//@param y : new y coordinate
			//@return  : true if polygon was closed, false otherwise
			bool resizePoly(const int x, const int y);

			//@brief    : stop resizing a polygon
			//@param pt : new coordinate
			void ungrabPoly() {if(!buildingPoly()) curPt = -1;}

			//@brief  : build an image mask for inside / outside the selection area
			//@param w: image width
			//@param h: image height
			std::vector<char> buildMask(const size_t w, const size_t h) const;

			//@brief  : check if a point is inside or outside the selection
			//@param x: x coordinate to check
			//@param y: y coordinate to check
			//@return : true/false if (x,y) is inside/outside the shape (false if no shape)
			bool inside(const int x, const int y) const;

			//@brief : convert an ROI to a string
			//@return: string representation
			std::string to_string() const;

			//@brief    : parse an ROI from a string
			//@param str: string representation
			void from_string(std::string str);

			//@brief : get if the mask is inverted or not
			//@return: true if the ROI describes the exluded section, false if it describes the included section
			bool getInv() const {return inv;}

			//@brief  : set if the mask is inverted or not
			//@param i: true if the ROI should describe the exluded section, false if it should describe the included section
			void setInv(const bool& i) {inv = i;}

	};

}//emsphinx

#include <limits>
#include <sstream>
#include <cmath>

namespace emsphinx {

	//@brief   : check if coordinates are near an existing point
	//@param x : coordinates to check
	//@param y : coordinates to check
	//@param r2: square radius cutoff
	//@return  : index of nearby point (-1 if none are nearby)
	int RoiSelection::nearPt(const int x, const int y, const int r2) const {
		if(pts.empty()) return -1;
		//loop over existing points finding closest one
		size_t nrPt = 0;
		int rMin = std::numeric_limits<int>::max();
		for(size_t i = 0; i < pts.size(); i++) {
			const int dx = pts[i].first  - x;
			const int dy = pts[i].second - y;
			const int rr = dx*dx + dy*dy;
			if(rr < rMin) {
				nrPt = i;
				rMin = rr;
			}
		}
		return rMin < r2 ? int(nrPt) : -1;
	}

	//@brief : clear current shape
	//@return: true if anything changed, false otherwise
	bool RoiSelection::clear() {
		inv = false;
		if(!pts.empty()) {
			pts.clear();
			curPt = -1;
			return true;
		}
		return false;
	}

	//@brief : chamge the drawing mode
	//@return: true if the mode changed, false otherwise
	bool RoiSelection::changeMode(const DrawMode dm) {
		if(mode != dm) {
			clear();
			mode = dm;
			return true;
		}
		return false;
	}

	//@brief    : update coordinates of a point
	//@param idx: point to update
	//@param x  : new x coordiante
	//@param y  : new y coordiante
	void RoiSelection::movePoint(const size_t idx, const int x, const int y) {
		if(pts.size() > 3 && !buildingPoly() && (0 == idx || pts.size() - 1 == idx)) {
			pts[0].first  = x;
			pts[0].second = y;
			pts.back() = pts.front();
		} else {
			pts[idx].first  = x;
			pts[idx].second = y;
		}
	}

	//rectangle building

	//@brief  : start drawing a new rectangle
	//@param x: x origin of new rectangle
	//@param y: y origin of new rectangle
	void RoiSelection::startRect(const int x, const int y) {
		pts.assign(2, std::pair<int, int>(x, y));
		curPt = 1;
	}

	//@brief  : resize a rectangle
	//@param x: new x coordinate
	//@param y: new y coordinate
	//@return : true if anything changed, false otherwise
	bool RoiSelection::resizeRect(const int x, const int y) {
		if(!pts.empty() && -1 != curPt) {
			std::pair<int, int> pt(x, y);
			if(mkSq) {
				int dx = x - pts[1-curPt].first ;
				int dy = y - pts[1-curPt].second;
				if(std::abs(dx) > std::abs(dy)) {
					dy  = dy < 0 ? -1 : 1;
					dy *= std::abs(dx);
				} else {
					dx  = dx < 0 ? -1 : 1;
					dx *= std::abs(dy);
				}
				pt = pts[1-curPt];
				pt.first  += dx;
				pt.second += dy;
			}
			pts[curPt] = pt;
			return true;
		}
		return false;
	}

	//polygon building

	//@brief  : start drawing a polygon
	//@param x: x origin of polygon
	//@param y: y origin of polygon
	void RoiSelection::startPoly(const int x, const int y) {
		pts.assign(2, std::pair<int, int>(x, y));
		curPt = 1;
	}

	//@brief   : add a new point to an existing polygon
	//@param x : x coordinate of new point to add
	//@param y : y coordinate of new point to add
	//@param r2: cutoff radius^2 for 2 points being the same
	//@note    : finishes the polygon if pt is sufficiently close to the first point
	//@return  : true if polygon was closed, false otherwise
	bool RoiSelection::addPolyPt(const int x, const int y, const int r2) {
		if(buildingPoly()) {
			//check if this point is close enough to first point
			int dx = x - pts.front().first ;
			int dy = y - pts.front().second;
			const bool finish = pts.size() > 2 && (dx*dx + dy * dy) < r2;

			if(finish) {//close polygon if near first point
				finishPoly();
				curPt = -1;
				return true;
			} else {//otherwise add new point
				dx = pts.back().first  - pts[pts.size() - 2].first ;
				dy = pts.back().second - pts[pts.size() - 2].second;
				const bool dup = (dx*dx + dy * dy) < r2;//is the last and 2nd to last point the same?
				if(!dup) pts.push_back(std::pair<int, int>(x, y));//don't add duplicate points (eg on up/down of mouse click)
				curPt = int(pts.size()) - 1;
				return false;
			}
		}
		return false;
	}

	//@brief : remove currently selected point from an existing polygon
	//@return: true if point was removed, false otherwise
	bool RoiSelection::remPolyPt() {
		if(-1 == curPt) return false;
		if(buildingPoly() && pts.size() > 2 ) {//currently building and can remove
			if(0 != curPt) std::swap(pts[curPt], pts[curPt-1]);
			pts.erase(pts.begin() + curPt);
			if(-1 != curPt) --curPt;
			return true;
		} else if(pts.size() > 4) {//completed and can remove
			pts.erase(pts.begin() + curPt);
			if(0 == curPt) pts.back() = pts.front();
			else if(pts.size() == curPt) pts.front() = pts.back();
			curPt = -1;
			return true;
		}
		return false;
	}

	//@brief      : duplicated currently selected point in an existing polygon
	//@param delta: offset from previous point
	//@return     : true if a point was duplicate, false otherwise
	bool RoiSelection::dupPolyPt(const int delta) {
		if(-1 == curPt) return false;//no point selected
		if(buildingPoly()) return false;//finish building first

		//determine line direction from previous point to current point
		int idxPrev = (curPt == 0) ? int(pts.size()) - 2 : curPt - 1;
		const int dx = pts[curPt].first  - pts[idxPrev].first ;
		const int dy = pts[curPt].second - pts[idxPrev].second;
		const int r2 = dx*dx + dy*dy;
		pts.insert(pts.begin() + curPt, pts[curPt]);
		++curPt;
		if(r2 > 0) {
			double rr = std::sqrt(double(r2));
			pts[curPt].first  += (int)std::round( double(dx * delta) / rr );
			pts[curPt].second += (int)std::round( double(dy * delta) / rr );
		}
		return true;
	}

	//@brief   : finish drawing a polygon
	//@return  : true if polygon was closed, false otherwise
	bool RoiSelection::finishPoly() {
		if(buildingPoly()) {//are we currently building a polygon?
			if(pts.back() != pts.front()) pts.back() = pts.front();//pts.push_back(pts.front());//add closing point if needed (TODO check)
			if(pts.size() < 4) pts.clear();//need at least a triangle
			curPt = -1;
			return true;
		}
		return false;
	}

	//@brief   : move a polygon vertex
	//@param x : new x coordinate
	//@param y : new y coordinate
	//@return  : true if polygon was closed, false otherwise
	bool RoiSelection::resizePoly(const int x, const int y) {
		if(!pts.empty() && -1 != curPt) {
			std::pair<int, int> pt(x, y);

			if(mkSq) {
				int dx = x - pts[curPt - 1].first ;
				int dy = y - pts[curPt - 1].second;
				if(std::abs(dx) > std::abs(dy)) dy = 0; else dx = 0;
				pt.first  = pts[curPt - 1].first  + dx;
				pt.second = pts[curPt - 1].second + dy;
			}
			pts[curPt] = pt;
			if(0 == curPt) pts.back() = pt;//keep closed if needed
			return true;
		}
		return false;
	}

	//@brief   : check if point c is above, below, or on the line defined by a ==> b
	//@param ax: x coordinate of point a
	//@param ay: y coordinate of point a
	//@param bx: x coordinate of point b
	//@param by: y coordinate of point b
	//@param cx: x coordinate of point c
	//@param cy: y coordinate of point c
	//@return  : a value >0, <0, ==0 for above, below, and on respectively
	//@note    : this isn't robust against roundoff errors if switched to floating point
	int orient2(int ax, int ay, int bx, int by, int cx, int cy) {return ( (bx - ax) * (cy - ay) - (cx -  ax) * (by - ay) );}

	//@brief  : build an image mask for inside / outside the selection area
	//@param w: image width
	//@param h: image height
	std::vector<char> RoiSelection::buildMask(const size_t w, const size_t h) const {
		//handle empty shape
		std::vector<char> mask(w * h, inv ? 1 : 0);//empty mask
		if(pts.empty()) return mask;//no points

		//start by getting bounding box of shape
		int xMin = pts.front().first , yMin = pts.front().second;
		int xMax = pts.front().first , yMax = pts.front().second;
		for(const std::pair<int, int>& pt : pts) {
			if(pt.first  < xMin) xMin = pt.first ;
			if(pt.second < yMin) yMin = pt.second;
			if(pt.first  > xMax) xMax = pt.first ;
			if(pt.second > yMax) yMax = pt.second;
		}

		//keep indices inside image bounds
		xMin = std::max(0, std::min(xMin, int(w)));
		yMin = std::max(0, std::min(yMin, int(h)));
		xMax = std::max(0, std::min(xMax, int(w)));
		yMax = std::max(0, std::min(yMax, int(h)));
		if(xMin == xMax || yMin == yMax) return mask;//no area inside image

		//now handle shapes specifically
		switch(mode) {
			case DrawMode::Rectangle://bounding box is mask
				for(int j = yMin; j < yMax; j++) std::fill(mask.begin() + j * w + xMin, mask.begin() + j * w + xMax, inv ? 0 : 1);
			break;

			case DrawMode::Ellipse  : {//do ellipse check for pixels within mask
				//compute center of ellipse once
				int x0 = pts.front().first  + pts.back().first ;//2 * x center
				int y0 = pts.front().second + pts.back().second;//2 * y center
				int aa = pts.front().first  - pts.back().first ;//2 * x axis length
				int bb = pts.front().second - pts.back().second;//2 * y axis length
				aa *= aa;//4 * a^2
				bb *= bb;//4 * b^2
				for(int j = yMin; j < yMax; j++) {
					int dy = j*2 - y0;//2*y distance from center
					for(int i = xMin; i < xMax; i++) {
						int dx = i*2 - x0;//2*x distance from center
						double rr = double(dx * dx) / aa + double(dy * dy) / bb;
						if(rr <= 1.0) mask[j * w + i] = inv ? 0 : 1;
					}
				}
			} break;

			case DrawMode::Polygon  ://compute winding number of pixels within mask
				for(int j = yMin; j < yMax; j++) {
					for(int i = xMin; i < xMax; i++) {

						//compute winding number of i,j
						//logic from http://geomalgorithms.com/a03-_inclusion.html
						int wn = 0;
						for(size_t k = 1; k < pts.size(); k++) {
							if(pts[k-1].second <= j) {// start j <= j
								if(pts[k].second > j)// an upward crossing
									if(orient2(pts[k-1].first , pts[k-1].second, pts[k].first , pts[k].second, i, j) > 0)  // P left of  edge
										++wn;            // have  a valid up intersect
							} else {// start j > j (no test needed)
								if(pts[k].second <= j)     // a downward crossing
									if(orient2(pts[k-1].first , pts[k-1].second, pts[k].first , pts[k].second, i, j) < 0)  // P right of  edge
										--wn;            // have  a valid down intersect
							}
						}
						if(wn != 0) mask[j * w + i] = inv ? 0 : 1;
					}
				}
			break;
		}
		return mask;
	}

	//@brief  : check if a point is inside or outside the selection
	//@param x: x coordinate to check
	//@param y: y coordinate to check
	//@return : true/false if (x,y) is inside/outside the shape (false if no shape)
	bool RoiSelection::inside(const int x, const int y) const {
		if(pts.empty()) return inv;//no shape

		//now handle shapes specifically
		switch(mode) {
			case DrawMode::Rectangle: {//bounding box is mask
				const int xL = std::min(pts.front().first , pts.back().first );
				const int yB = std::min(pts.front().second, pts.back().second);
				const int xR = std::max(pts.front().first , pts.back().first );
				const int yT = std::max(pts.front().second, pts.back().second);
				return (x > xL && x < xR && y > yB && y < yT) ? !inv : inv;
			} break;

			case DrawMode::Ellipse  : {//do ellipse check for pixels within mask
				const int xL = std::min(pts.front().first , pts.back().first );
				const int yB = std::min(pts.front().second, pts.back().second);
				const int xR = std::max(pts.front().first , pts.back().first );
				const int yT = std::max(pts.front().second, pts.back().second);
				if(! ( x > xL && x < xR && y > yB && y < yT) ) return inv;//outside bounding box


				const int x0 = xR + xL;//2*x center
				const int y0 = yT + yB;//2*y center
				const int a  = xR - xL;//2*x axis length
				const int b  = yT - yB;//2*y axis length
				const int dx = x * 2 - x0;
				const int dy = y * 2 - y0;
				const double rr = double(dx * dx) / (a * a) + double(dy * dy) / (b * b);
				return (rr <= 1.0) ? !inv : inv;
			} break;

			case DrawMode::Polygon  ://compute winding number of pixels within mask
				//compute winding number of x, y
				//logic from http://geomalgorithms.com/a03-_inclusion.html
				int wn = 0;
				for(size_t i = 1; i < pts.size(); i++) {
					if(pts[i-1].second <= y) {// start y <= y
						if(pts[i].second > y)// an upward crossing
							if(orient2(pts[i-1].first , pts[i-1].second, pts[i].first , pts[i].second, x, y) > 0)  // P left of  edge
								++wn;            // have  a valid up intersect
					} else {// start y > y (no test needed)
						if(pts[i].second <= y)     // a downward crossing
							if(orient2(pts[i-1].first , pts[i-1].second, pts[i].first , pts[i].second, x, y) < 0)  // P right of  edge
								--wn;            // have  a valid down intersect
					}
				}
				return (wn != 0) ? !inv : inv;
			break;
		}
		return inv;
	}

	//@brief : convert an ROI to a string
	//@return: string representation
	std::string RoiSelection::to_string() const {
		if(!hasShape()) return "";
		std::ostringstream ss;
		if(inv) ss << 'i';
		switch(mode) {
			case DrawMode::Ellipse  : ss << 'e';
			case DrawMode::Rectangle:
				if(2 != pts.size()) return "";
				ss <<                pts[0].first << ", " <<                 pts[0].second << ", ";
				ss << pts[1].first - pts[0].first << ", " << pts[1].second - pts[0].second;
			break;

			case DrawMode::Polygon  :
				for (size_t i = 0; i < pts.size(); i++) {
					ss << pts[i].first << ", " << pts[i].second;
					if (pts.size() != i + 1) ss << ", ";
				}
			break;
		}
		return ss.str();
	}

	//@brief    : parse an ROI from a string
	//@param str: string representation
	void RoiSelection::from_string(std::string str) {
		clear();
		if(str.empty()) return;
		if("0" == str) return;
		std::istringstream ss(str);
		mode = DrawMode::Rectangle;
		if('i' == ss.peek()) {
			inv = true;
			char skip = ss.get();
		}
		if('e' == ss.peek()) {
			mode = DrawMode::Ellipse;
			char skip = ss.get();
		}

		char c = ',';
		int i;
		std::vector<int> v;
		while(ss >> i) {
			v.push_back(i);
			if(ss >> c) {
				if(',' != c) throw std::runtime_error("expected ',' between coordinates in ROI string");
			}
		}
		if(0 != v.size() % 2) throw std::runtime_error("odd number of points in ROI string");
		if(4 == v.size()) {
			pts.resize(2);
			pts[0].first = v[0]       ; pts[0].second = v[1]       ;
			pts[1].first = v[0] + v[2]; pts[1].second = v[1] + v[3];
		} else {//polygon
			if(DrawMode::Ellipse == mode) throw std::runtime_error("too many points for ellipse in ROI string");
			mode = DrawMode::Polygon;
			pts.resize(v.size() / 2);
			for(size_t i = 0; i < pts.size(); i++) {
				pts[i].first = v[2*i]; pts[i].second = v[2*i+1];
			}
			if(pts.front() != pts.back()) throw std::runtime_error("polygon not closed in ROI string");
		}
	}

}//emsphinx


#endif//_ROI_H_
