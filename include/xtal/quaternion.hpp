/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                     *
 * Copyright (c) 2019, William C. Lenthe                               *
 * All rights reserved.                                                *
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
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef _quaternion_h_
#define _quaternion_h_

#include <string>
#include <ostream>

namespace xtal {

	////////////////////////////////////////////////////////////////////////
	//                      Quaternion Wrapper Class                      //
	////////////////////////////////////////////////////////////////////////

	//POD struct instead of class + packed check allows casting from array of doubles to array of quats ([w1,x1,y1,z1,w2,x2,y2,z2,w3,x3,y3,z3] -> [q1,q2,q3])
	template<typename Real>
	struct Quat {
		Real w, x, y, z;

		//@brief: default constructor
		Quat();

		//@brief   : construct from values
		//@param vw: value of w
		//@param vx: value of x
		//@param vy: value of y
		//@param vz: value of z
		Quat(const Real vw, const Real vx, const Real vy, const Real vz);

		//@brief : construct a quaternion of zeros
		//@return: (0, 0, 0, 0)
		static Quat Zero    ();

		//@brief : construct an identity quaternion
		//@return: (1, 0, 0, 0)
		static Quat Identity();

		//@brief : get a pointer to the underlying data as {w, x, y, z}
		//@return: pointer to data
		Real      *const data()      ;
		Real const*const data() const;

		//unary operations
		Quat  conj       ()             const;
		Quat  inv        ()             const;
		Quat  neg        ()             const;
		Quat  cAbs       ()             const;
		Quat  expl       ()             const;
		Quat  normalize  ()             const;
		Real  norm       ()             const;
		Real  norm2      ()             const;
		Quat  operator-  ()             const;

		//binary quaternion scalar operations
		Quat& operator+=(const Real& s)      ;
		Quat& operator-=(const Real& s)      ;
		Quat& operator*=(const Real& s)      ;
		Quat& operator/=(const Real& s)      ;
		Quat  operator+ (const Real& s) const;
		Quat  operator- (const Real& s) const;
		Quat  operator* (const Real& s) const;
		Quat  operator/ (const Real& s) const;

		//binary quaternion quaternion operations
		Real  dot       (const Quat& q) const;
		Quat& operator+=(const Quat& q)      ;
		Quat& operator-=(const Quat& q)      ;
		Quat& operator*=(const Quat& q)      ;
		Quat& operator/=(const Quat& q)      ;
		Quat  operator+ (const Quat& q) const;
		Quat  operator- (const Quat& q) const;
		Quat  operator* (const Quat& q) const;
		Quat  operator/ (const Quat& q) const;
		bool  operator< (const Quat& q) const;
		bool  operator==(const Quat& q) const;

		//vector rotation
		void rotateVector(Real const * const vIn, Real * const vOut) const;

		//string conversion
		std::string to_string(const size_t dig = std::numeric_limits<Real>::digits10, const bool pos = false) const;

		//@brief   : formatted output
		//@param os: location to write formatted quaternion
		//@param qu: quaternion or write
		//@return  : os with quaternion written to it
		template <typename T> friend std::ostream& operator<<(std::ostream& os, const Quat<T>& qu);
	};

	////////////////////////////////////////////////////////////////////////////////////
	//       underlying quaternion math functions (operate on pointer to wxyz)        //
	////////////////////////////////////////////////////////////////////////////////////

	namespace quat {

		////////////////////////////////////////////////////////////////////////
		//                    Unary Quaternion Operations                     //
		////////////////////////////////////////////////////////////////////////

		//@brief     : urnary quaternion operations with a quaternion result
		//@param qIn : quaternion
		//@param qOut: location to write operator(qIn)
		template <typename Real> void conj     (Real const * const qIn                          , Real * const qOut);//conjugate
		template <typename Real> void inv      (Real const * const qIn                          , Real * const qOut);//inverse
		template <typename Real> void neg      (Real const * const qIn                          , Real * const qOut);//negate
		template <typename Real> void cAbs     (Real const * const qIn                          , Real * const qOut);//element wise absolute value
		template <typename Real> void expl     (Real const * const qIn                          , Real * const qOut);//explement (positive w equivalent)
		template <typename Real> void normalize(Real const * const qIn                          , Real * const qOut);//qIn / norm(qIn)

		//@brief     : urnary quaternion operations with a scalar result
		//@param qIn : quaternion
		//@return    : result
		template <typename Real> Real norm     (Real const * const qIn                                             );//magnitude
		template <typename Real> Real norm2    (Real const * const qIn                                             );//squared magnitude

		////////////////////////////////////////////////////////////////////////
		//                    Binary Quaternion Operations                    //
		////////////////////////////////////////////////////////////////////////

		//@brief     : binary operations between a quaternion and a scalar
		//@param qIn : input quaternion
		//@param s   : scalar to element wise operate with
		//@param qOut: location to write qIn (operator) s
		template <typename Real> void scalarAdd(Real const * const qIn , Real const         s   , Real * const qOut);
		template <typename Real> void scalarSub(Real const * const qIn , Real const         s   , Real * const qOut);
		template <typename Real> void scalarMul(Real const * const qIn , Real const         s   , Real * const qOut);
		template <typename Real> void scalarDiv(Real const * const qIn , Real const         s   , Real * const qOut);

		//@brief     : binary element wise operations between 2 quaternions with a scalar result
		//@param qIn1: first quaternion
		//@param qIn2: second quaternion
		//@return    : qIn1 (operator) qIn2
		template <typename Real> Real dot      (Real const * const qIn1, Real const * const qIn2                   );
		template <typename Real> bool less     (Real const * const qIn1, Real const * const qIn2                   );
		template <typename Real> bool equal    (Real const * const qIn1, Real const * const qIn2                   );

		//@brief     : binary element wise operations between 2 quaternions with a quaternion result
		//@param qIn1: first quaternion
		//@param qIn2: second quaternion
		//@param qOut: location to write qIn1 (operator) qIn2
		template <typename Real> void add      (Real const * const qIn1, Real const * const qIn2, Real * const qOut);
		template <typename Real> void sub      (Real const * const qIn1, Real const * const qIn2, Real * const qOut);
		template <typename Real> void mul      (Real const * const qIn1, Real const * const qIn2, Real * const qOut);//proper quaternion multiplication
		template <typename Real> void div      (Real const * const qIn1, Real const * const qIn2, Real * const qOut);//proper quaternion division: such that qu * qr / qr == qu

		////////////////////////////////////////////////////////////////////////
		//                     Other Quaternion Functions                     //
		////////////////////////////////////////////////////////////////////////

		//@brief     : actively rotate a vector by a quaternion (q * v * q.conj())
		//@param q   : quaternion to rotate by
		//@param vIn : vector to rotate
		//@param vOut: location to write rotated vector
		template <typename Real> void rotateVector(Real const * const q, Real const * const vIn, Real * const vOut);

		//@brief    : convert quaternion to a nicely formatted string
		//@param qu : quaternion to convert to string
		//@param dig: significant digits to print with
		//@param pos: should a '+' be used before positive numbers
		//@return   : nicely formatted string as w x y z
		template <typename Real> std::string to_string(Real const * const qu, const size_t dig = std::numeric_limits<Real>::digits10, const bool pos = false);
	}
}

////////////////////////////////////////////////////////////////////////
//                          Implementations                           //
////////////////////////////////////////////////////////////////////////

#include "constants.hpp"

#include <algorithm>
#include <functional>
#include <numeric>
#include <cmath>
#include <iomanip>
#include <sstream>

namespace xtal {

	////////////////////////////////////////////////////////////////////////
	//                      Quaternion Class Members                      //
	////////////////////////////////////////////////////////////////////////

	//@brief: default constructor
	template <typename Real>
	Quat<Real>::Quat() {
		static_assert(sizeof(Quat) == 4 * sizeof(Real), "Quaternion struct must be packed"                   );
		static_assert(std::is_floating_point<Real>::value   , "Quaternion must be templated on floating point type");
	}

	//@brief   : construct from values
	//@param vw: value of w
	//@param vx: value of x
	//@param vy: value of y
	//@param vz: value of z
	template <typename Real> Quat<Real>::Quat(const Real vw, const Real vx, const Real vy, const Real vz) : w(vw), x(vx), y(vy), z(vz) {}

	//@brief : construct a quaternion of zeros
	//@return: (0, 0, 0, 0)
	template <typename Real> Quat<Real> Quat<Real>::Zero    () {return Quat(0, 0, 0, 0);}

	//@brief : construct an identity quaternion
	//@return: (1, 0, 0, 0)
	template <typename Real> Quat<Real> Quat<Real>::Identity() {return Quat(1, 0, 0, 0);}

	//@brief : get a pointer to the underlying data as {w, x, y, z}
	//@return: pointer to data
	template <typename Real> Real      *const Quat<Real>::data ()       {return (Real      *const)this;}
	template <typename Real> Real const*const Quat<Real>::data () const {return (Real const*const)this;}

	//unary operations
	template <typename Real> Quat<Real> Quat<Real>::conj       ()              const {Quat q; quat::conj     (data(),           q.data()); return q;    }
	template <typename Real> Quat<Real> Quat<Real>::inv        ()              const {Quat q; quat::inv      (data(),           q.data()); return q;    }
	template <typename Real> Quat<Real> Quat<Real>::neg        ()              const {Quat q; quat::neg      (data(),           q.data()); return q;    }
	template <typename Real> Quat<Real> Quat<Real>::cAbs       ()              const {Quat q; quat::cAbs     (data(),           q.data()); return q;    }
	template <typename Real> Quat<Real> Quat<Real>::expl       ()              const {Quat q; quat::expl     (data(),           q.data()); return q;    }
	template <typename Real> Quat<Real> Quat<Real>::normalize  ()              const {Quat q; quat::normalize(data(),           q.data()); return q;    }
	template <typename Real>      Real  Quat<Real>::norm       ()              const { return quat::norm     (data()                    );              }
	template <typename Real>      Real  Quat<Real>::norm2      ()              const { return quat::norm2    (data()                    );              }
	template <typename Real> Quat<Real> Quat<Real>::operator-  ()              const { return       neg      (                          );              }

	//binary quaternion scalar operations
	template <typename Real> Quat<Real>& Quat<Real>::operator+=(const Real& s)       {        quat::scalarAdd(data(), s       ,   data()); return *this;}
	template <typename Real> Quat<Real>& Quat<Real>::operator-=(const Real& s)       {        quat::scalarSub(data(), s       ,   data()); return *this;}
	template <typename Real> Quat<Real>& Quat<Real>::operator*=(const Real& s)       {        quat::scalarMul(data(), s       ,   data()); return *this;}
	template <typename Real> Quat<Real>& Quat<Real>::operator/=(const Real& s)       {        quat::scalarDiv(data(), s       ,   data()); return *this;}
	template <typename Real> Quat<Real>  Quat<Real>::operator+ (const Real& s) const {Quat q; quat::scalarAdd(data(), s       , q.data()); return q;    }
	template <typename Real> Quat<Real>  Quat<Real>::operator- (const Real& s) const {Quat q; quat::scalarSub(data(), s       , q.data()); return q;    }
	template <typename Real> Quat<Real>  Quat<Real>::operator* (const Real& s) const {Quat q; quat::scalarMul(data(), s       , q.data()); return q;    }
	template <typename Real> Quat<Real>  Quat<Real>::operator/ (const Real& s) const {Quat q; quat::scalarDiv(data(), s       , q.data()); return q;    }

	//binary quaternion quaternion operations
	template <typename Real>      Real   Quat<Real>::dot       (const Quat& q) const { return quat::dot      (data(), q.data()          );              }
	template <typename Real> Quat<Real>& Quat<Real>::operator+=(const Quat& q)       {        quat::add      (data(), q.data(),   data()); return *this;}
	template <typename Real> Quat<Real>& Quat<Real>::operator-=(const Quat& q)       {        quat::sub      (data(), q.data(),   data()); return *this;}
	template <typename Real> Quat<Real>& Quat<Real>::operator*=(const Quat& q)       {        quat::mul      (data(), q.data(),   data()); return *this;}
	template <typename Real> Quat<Real>& Quat<Real>::operator/=(const Quat& q)       {        quat::div      (data(), q.data(),   data()); return *this;}
	template <typename Real> Quat<Real>  Quat<Real>::operator+ (const Quat& q) const {Quat r; quat::add      (data(), q.data(), r.data()); return r;    }
	template <typename Real> Quat<Real>  Quat<Real>::operator- (const Quat& q) const {Quat r; quat::sub      (data(), q.data(), r.data()); return r;    }
	template <typename Real> Quat<Real>  Quat<Real>::operator* (const Quat& q) const {Quat r; quat::mul      (data(), q.data(), r.data()); return r;    }
	template <typename Real> Quat<Real>  Quat<Real>::operator/ (const Quat& q) const {Quat r; quat::div      (data(), q.data(), r.data()); return r;    }
	template <typename Real> bool        Quat<Real>::operator< (const Quat& q) const { return quat::less     (data(), q.data()           );             }
	template <typename Real> bool        Quat<Real>::operator==(const Quat& q) const { return quat::equal    (data(), q.data()           );             }

	//vector rotation
	template <typename Real> void Quat<Real>::rotateVector(Real const * const vIn, Real * const vOut) const {quat::rotateVector(data(), vIn, vOut);}

	//string conversion
	template <typename Real> std::string Quat<Real>::to_string(const size_t dig, const bool pos) const {return quat::to_string(data(), dig, pos);}

	//@brief   : formatted output
	//@param os: location to write formatted quaternion
	//@param qu: quaternion or write
	//@return  : os with quaternion written to it
	template <typename Real> std::ostream& operator<<(std::ostream& os, const Quat<Real>& qu) {return os << qu.to_string(6);}

	////////////////////////////////////////////////////////////////////////////////////
	//       underlying quaternion math functions (operate on pointer to wxyz)        //
	////////////////////////////////////////////////////////////////////////////////////
	
	namespace quat {
		////////////////////////////////////////////////////////////////////////
		//                    Unary Quaternion Operations                     //
		////////////////////////////////////////////////////////////////////////

		//@brief     : quaternion conjugate
		//@param qIn : quaternion
		//@param qOut: location to write conjugate(qIn)
		template <typename Real>
		void conj(Real const * const qIn, Real * const qOut) {
			qOut[0] = qIn[0];//w is unchanged
			std::transform(qIn+1, qIn+4, qOut+1, std::negate<Real>());//x, y, and z are negated
		}

		//@brief     : quaternion inverse
		//@param qIn : quaternion
		//@param qOut: location to write inverse(qIn)
		template <typename Real>
		void inv(Real const * const qIn, Real * const qOut) {
			conj(qIn, qOut);
			scalarDiv(qOut, norm2(qOut), qOut);//inv(q) = conj(q) / norm(q)
		}

		//@brief     : quaternion negate
		//@param qIn : quaternion
		//@param qOut: location to write -qIn
		template <typename Real>
		void neg (Real const * const qIn, Real * const qOut) {
			std::transform(qIn, qIn+4, qOut, std::negate<Real>());
		}

		//@brief     : element wise absolute value
		//@param qIn : quaternion
		//@param qOut: location to write element wise abs(qIn)
		template <typename Real>
		void cAbs(Real const * const qIn, Real * const qOut) {
			std::transform(qIn, qIn+4, qOut, static_cast<Real(*)(Real)>(&std::fabs));
		}

		//@brief     : quaternion explement
		//@param qIn : quaternion
		//@param qOut: location to write operator(qIn)
		//@note      : complementary angles sum to 90, supplementary to 180, and explary to 360
		template <typename Real>
		void expl(Real const * const qIn, Real * const qOut) {
			if(std::signbit(qIn[0])) {
				neg(qIn, qOut);
			} else {
				std::copy(qIn, qIn+4, qOut);
			}
		}

		//@brief     : normalize a quaternion
		//@param qIn : quaternion to normalize
		//@param qOut: location to write normalized quaternion
		template <typename Real>
		void normalize(Real const * const qIn                          , Real * const qOut) {
			scalarDiv(qIn, norm(qIn), qOut);
		}

		//@brief     : quaternion norm (absolute value) 
		//@param qIn : quaternion
		//@return    : norm
		template <typename Real>
		Real norm     (Real const * const qIn                                             ) {
			return std::sqrt(norm2(qIn));
		}
		//@brief     : quaternion norm^2 (absolute value) 
		//@param qIn : quaternion
		//@return    : norm^2
		template <typename Real>
		Real norm2    (Real const * const qIn                                             ) {
			return dot(qIn, qIn);
		}

		////////////////////////////////////////////////////////////////////////
		//                    Binary Quaternion Operations                    //
		////////////////////////////////////////////////////////////////////////

		//@brief     : element wise add to a quaternion
		//@param qIn : input quaternion
		//@param s   : scalar to element wise operate with
		//@param qOut: location to write qIn (operator) s
		template <typename Real>
		void scalarAdd(Real const * const qIn , Real const         s   , Real * const qOut) {
			std::transform(qIn , qIn +4,       qOut, [&s](const Real i){return i+s;});
		}

		//@brief     : element wise subtract from to a quaternion
		//@param qIn : input quaternion
		//@param s   : scalar to element wise operate with
		//@param qOut: location to write qIn (operator) s
		template <typename Real>
		void scalarSub(Real const * const qIn , Real const         s   , Real * const qOut) {
			std::transform(qIn , qIn +4,       qOut, [&s](const Real i){return i-s;});
		}

		//@brief     : element wise multiply with a quaternion
		//@param qIn : input quaternion
		//@param s   : scalar to element wise operate with
		//@param qOut: location to write qIn (operator) s
		template <typename Real>
		void scalarMul(Real const * const qIn , Real const         s   , Real * const qOut) {
			std::transform(qIn , qIn +4,       qOut, [&s](const Real i){return i*s;});
		}

		//@brief     : element wise divide into a quaternion
		//@param qIn : input quaternion
		//@param s   : scalar to element wise operate with
		//@param qOut: location to write qIn (operator) s
		template <typename Real>
		void scalarDiv(Real const * const qIn , Real const         s   , Real * const qOut) {
			std::transform(qIn , qIn +4,       qOut, [&s](const Real i){return i/s;});
		}

		//@brief     : element wise add 2 quaternions
		//@param qIn1: first quaternion
		//@param qIn2: second quaternion
		//@param qOut: location to write qIn1 + qIn2
		template <typename Real>
		void add      (Real const * const qIn1, Real const * const qIn2, Real * const qOut) {
			std::transform(qIn1, qIn1+4, qIn2, qOut, std::plus <Real>());
		}

		//@brief     : element wise subtract 2 quaternions
		//@param qIn1: first quaternion
		//@param qIn2: second quaternion
		//@param qOut: location to write qIn1 - qIn2
		template <typename Real>
		void sub      (Real const * const qIn1, Real const * const qIn2, Real * const qOut) {
			std::transform(qIn1, qIn1+4, qIn2, qOut, std::minus<Real>());
		}

		//@brief     : multiply two quaternions
		//@param qIn1: first quaternion
		//@param qIn2: second quaternion
		//@param qOut: location to write qIn1 * qIn2
		template <typename Real>
		void mul      (Real const * const qIn1, Real const * const qIn2, Real * const qOut) {
			const Real w = qIn1[0]*qIn2[0] - qIn1[1]*qIn2[1] -  qIn1[2]*qIn2[2] - qIn1[3]*qIn2[3]        ;
			const Real x = qIn1[0]*qIn2[1] + qIn1[1]*qIn2[0] + (qIn1[2]*qIn2[3] - qIn1[3]*qIn2[2]) * pijk;
			const Real y = qIn1[0]*qIn2[2] + qIn1[2]*qIn2[0] + (qIn1[3]*qIn2[1] - qIn1[1]*qIn2[3]) * pijk;
			const Real z = qIn1[0]*qIn2[3] + qIn1[3]*qIn2[0] + (qIn1[1]*qIn2[2] - qIn1[2]*qIn2[1]) * pijk;
			qOut[0] = w;
			qOut[1] = x;
			qOut[2] = y;
			qOut[3] = z;
		}

		//@brief     : divde two quaternions
		//@param qIn1: first quaternion
		//@param qIn2: second quaternion
		//@param qOut: location to write qIn1 / qIn2
		template <typename Real>
		void div      (Real const * const qIn1, Real const * const qIn2, Real * const qOut) {
			// qIn1 / qIn2 == qIn1 * conj(qIn2) / norm2(qIn2)
			const Real v =   qIn2[0]*qIn2[0] + qIn2[1]*qIn2[1] +  qIn2[2]*qIn2[2] + qIn2[3]*qIn2[3]              ; // norm2(qIn2)
			const Real w = ( qIn1[0]*qIn2[0] + qIn1[1]*qIn2[1] +  qIn1[2]*qIn2[2] + qIn1[3]*qIn2[3]        ) / v;
			const Real x = (-qIn1[0]*qIn2[1] + qIn1[1]*qIn2[0] - (qIn1[2]*qIn2[3] - qIn1[3]*qIn2[2]) * pijk) / v;
			const Real y = (-qIn1[0]*qIn2[2] + qIn1[2]*qIn2[0] - (qIn1[3]*qIn2[1] - qIn1[1]*qIn2[3]) * pijk) / v;
			const Real z = (-qIn1[0]*qIn2[3] + qIn1[3]*qIn2[0] - (qIn1[1]*qIn2[2] - qIn1[2]*qIn2[1]) * pijk) / v;
			qOut[0] = w;
			qOut[1] = x;
			qOut[2] = y;
			qOut[3] = z;
		}

		//@brief     : dot product of 2 quaternions (sum of element wise product)
		//@param qIn1: first quaternion
		//@param qIn2: second quaternion
		//@return    : qIn1.qIn2
		template <typename Real>
		Real dot (Real const * const qIn1, Real const * const qIn2                   ) {
			return std::inner_product(qIn1, qIn1+4, qIn2, Real(0));
		}

		//@brief     : check if a quaternion is less than another
		//@param qIn1: first quaternion
		//@param qIn2: second quaternion
		//@return    : qIn1 < qIn2
		//@note      : lexicographical compare
		template <typename Real>
		bool less     (Real const * const qIn1, Real const * const qIn2                   ) {
			return std::lexicographical_compare(qIn1, qIn1 + 4, qIn2, qIn2 + 4);
		}

		//@brief     : check if two quaternions are equal
		//@param qIn1: first quaternion
		//@param qIn2: second quaternion
		//@return    : qIn1 == qIn2
		template <typename Real>
		bool equal    (Real const * const qIn1, Real const * const qIn2                   ) {
			return std::equal(qIn1, qIn1 + 4, qIn2);
		}

		////////////////////////////////////////////////////////////////////////
		//                     Other Quaternion Functions                     //
		////////////////////////////////////////////////////////////////////////

		//@brief     : actively rotate a vector by a quaternion (q * v * q.conj())
		//@param q   : quaternion to rotate by
		//@param vIn : vector to rotate
		//@param vOut: location to write rotated vector
		template <typename Real>
		void rotateVector(Real const * const q, Real const * const vIn, Real * const vOut) {
			//q * v [stored locally since vIn/Out could overlap and vOut could only have 3 elements]
			const Real w = -(q[1]*vIn[0] +  q[2]*vIn[1] + q[3]*vIn[2])       ;
			const Real x =   q[0]*vIn[0] + (q[2]*vIn[2] - q[3]*vIn[1]) * pijk;
			const Real y =   q[0]*vIn[1] + (q[3]*vIn[0] - q[1]*vIn[2]) * pijk;
			const Real z =   q[0]*vIn[2] + (q[1]*vIn[1] - q[2]*vIn[0]) * pijk;

			//(q * v) * q.conj() [stored locally in case q and vOut overlap]
			const Real vx = x*q[0] - w*q[1] + (z*q[2] - y*q[3]) * pijk;
			const Real vy = y*q[0] - w*q[2] + (x*q[3] - z*q[1]) * pijk;
			const Real vz = z*q[0] - w*q[3] + (y*q[1] - x*q[2]) * pijk;

			//save output 
			vOut[0] = vx;
			vOut[1] = vy;
			vOut[2] = vz;
		}

		//@brief    : convert quaternion to a nicely formatted string
		//@param qu : quaternion to convert to string
		//@param dig: significant digits to print with
		//@param pos: should a '+' be used before positive numbers
		//@return   : nicely formatted string as w x y z
		template <typename Real>
		std::string to_string(Real const * const qu, const size_t dig, const bool pos) {
			std::stringstream ss;
			ss << std::fixed       ;//decimal (vs scientific) output
			ss << std::showpoint   ;//always show the decimal for floating point numbers
			ss << std::left        ;//fill unused width to the right
			ss << std::setfill('0');//pad with 0 instead of spaces
			if(pos) ss << std::showpos;
			if(!pos && !std::signbit(qu[0])) ss << ' ';//leave space so that '-' doesn't cause misalignment
			ss << std::setw(dig) << qu[0] << ' ';
			if(!pos && !std::signbit(qu[1])) ss << ' ';//leave space so that '-' doesn't cause misalignment
			ss << std::setw(dig) << qu[1] << ' ';
			if(!pos && !std::signbit(qu[2])) ss << ' ';//leave space so that '-' doesn't cause misalignment
			ss << std::setw(dig) << qu[2] << ' ';
			if(!pos && !std::signbit(qu[3])) ss << ' ';//leave space so that '-' doesn't cause misalignment
			ss << std::setw(dig) << qu[3];
			return ss.str();
		}
	}
}

#endif//_quaternion_h_
