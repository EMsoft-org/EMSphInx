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

#ifndef _NML_H_
#define _NML_H_

#include <string>
#include <istream>
#include <vector>
#include <map>

#define NML_USE_H5
#ifdef NML_USE_H5
	#include "H5Cpp.h"
#endif

namespace nml {
	enum class Type {
		None  ,
		Bool  ,
		Int   ,
		Double,
		String
	};

	//a (bad) variant class as a place holder until I can use std::variant<bool, int, float, std::string> (c++17)
	//could probably use a union here but the memory overhead should be minimal right now
	struct Variant {
		//@brief: construct an empty variant
		Variant() {clear();}

		//@brief: construct a variant from a value
		template <typename T> Variant(const T& v) {set(v);}

		//@brief: clear the stored type
		void clear() {t = Type::None; f = true;}

		//@brief : get the the stored value
		//@return: stored value
		//@note  : throws if the return type isn't the same as (or nicely castable to) the stored type (e.g. can't get<double> if t != Double)
		bool        getBool  ();
		int         getInt   ();
		double      getDouble();
		std::string getString();

		//@brief  : set the stored value and update the type if needed
		//@param v: value to store
		template <typename T> void set(const T& v);

		//@brief: get the stored type
		Type getType() const {return t;}

		//@brief : check if the value has been read via get() since setting
		//@return: true/false if the value has / hasn't been read
		//@note  : always returns true if Type is None
		bool used() const {return f;}

		//@brief : convert the currently held value to a string representatino
		//@return: string representation
		std::string to_string() const;

		private:
			//storage for value
			union {
				bool   b;
				int    i;
				double d;
			} u;//at least putting these in a union save a bit of memory
			std::string s;

			Type t;//flag for stored value type
			bool f;//flag for if the value has been read since setting
	};

	//a token is a name value pair with the value a vector of variants
	struct Value : public std::vector<Variant> {
		//@brief  : get the type of the currently stored value
		Type getType() const {return this->empty() ? Type::None : this->front().getType();}

		//@brief : check if the current value has been read
		//@return: true if the current value was read since setting, false otherwise
		bool wasUsed() const;// {return std::all_of(this->cbegin(), this->cend(), [](const Variant& v) {return v.used();});}
	};

	//a namelist is a collection of key/value pairs (with case insensitive keys)
	class NameList : private std::map<std::string, Value> {
		std::vector<std::string> fileLines;//original file

		//@brief  : convert a string to lower case
		//@param s: string to convert
		//@return : lowercase string
		static std::string ToLower(std::string s);// {std::transform(s.begin(), s.end(), s.begin(), [](char c){return std::tolower(c);}); return s;}

		//@brief    : find the first token with a given key
		//@param nm : key to search for
		//@return   : token with specified key (throws if not found)
		//@note     : write acces isn't public facing
		Value      & at(const std::string nm);

		public:
			//@brief: construct an empty namelist
			NameList() {}

			//@brief   : construct a namelist from a file
			//@param nm: file name to read from
			NameList(const std::string nm) {read(nm);}

			//@brief   : read a namelist file
			//@param is: istream to read namelist from
			//@param cm: comment characters (ignore lines if the first character after white space is in cm)
			void read(std::istream& is, std::string cm = "!");

			//@brief   : read a namelist file
			//@param nm: filename of namelist to read
			//@param cm: comment characters (ignore lines if the first character after white space is in cm)
			void read(std::string nm, std::string cm = "!");

			//@brief    : find the first token with a given key
			//@param nm : key to search for
			//@return   : token with specified key (throws if not found)
			Value const& at(const std::string nm) const;

			//@brief    : attempt to parse a bool from a list of namelist parameters
			//@param nm : key to search for
			//@return   : parsed bool (throws if not found)
			bool        getBool  (std::string nm) {return at(nm)[0].getBool  ();}

			//@brief    : attempt to parse a int from a list of namelist parameters
			//@param nm : key to search for
			//@return   : parsed int (throws if not found)
			int         getInt   (std::string nm) {return at(nm)[0].getInt   ();}

			//@brief    : attempt to parse a double from a list of namelist parameters
			//@param nm : key to search for
			//@return   : parsed double (throws if not found)
			double      getDouble(std::string nm) {return at(nm)[0].getDouble();}

			//@brief    : attempt to parse a string from a list of namelist parameters
			//@param nm : key to search for
			//@return   : parsed string (throws if not found)
			std::string getString(std::string nm) {return at(nm)[0].getString();}

			//@brief    : attempt to parse a list fo bool from a list of namelist parameters
			//@param nm : key to search for
			//@return   : parsed bools (throws if not found)
			std::vector<bool       > getBools  (std::string nm);

			//@brief    : attempt to parse a list fo int from a list of namelist parameters
			//@param nm : key to search for
			//@return   : parsed ints (throws if not found)
			std::vector<int        > getInts   (std::string nm);

			//@brief    : attempt to parse a list fo double from a list of namelist parameters
			//@param nm : key to search for
			//@return   : parsed doubles (throws if not found)
			std::vector<double     > getDoubles(std::string nm);

			//@brief    : attempt to parse a list fo string from a list of namelist parameters
			//@param nm : key to search for
			//@return   : parsed strings (throws if not found)
			std::vector<std::string> getStrings(std::string nm);

			//@brief : check if there are unused tokens
			//@return: false if all tokens have been parse (via get___), false otherwise
			bool fullyParsed() const;

			//@brief : get a list of unused tokens
			//@return: names of all unused tokens
			std::string unusedTokens() const;

		#ifdef NML_USE_H5
			//@brief    : write all parsed tokens to an hdf file
			//@param grp: hgf group to write values to
			void writeParameters(H5::Group grp);

			//@brief    : write the raw file to a dataset
			//@param grp: hgf group to create dataset in
			//@param nm : name of dataset to write
			void writeFile(H5::Group grp, std::string nm) const;
		#endif
	};
}

////////////////////////////////////////////////////////////////////////////////
//                           Implementation Details                           //
////////////////////////////////////////////////////////////////////////////////

#include <fstream>
#include <sstream>
#include <algorithm>
#include <cctype>

namespace nml {

	bool        Variant::getBool  () {if(t == Type::Bool  ) {f = true; return u.b;}                                            throw std::runtime_error("stored type isn't boolean");}
	int         Variant::getInt   () {if(t == Type::Int   ) {f = true; return u.i;}                                            throw std::runtime_error("stored type isn't integer");}
	double      Variant::getDouble() {if(t == Type::Double) {f = true; return u.d;} if(t == Type::Int) {f = true; return u.i;} throw std::runtime_error("stored type isn't double" );}//cast ints to doubles silently
	std::string Variant::getString() {if(t == Type::String) {f = true; return   s;}                                            throw std::runtime_error("stored type isn't string" );}
	// std::string Variant::getString() {f = true; return to_string();}

	// template <> bool        Variant::get() {if(t == Type::Bool  ) {f = true; return u.b;} throw std::runtime_error("stored type isn't boolean");}
	// template <> int         Variant::get() {if(t == Type::Int   ) {f = true; return u.i;} throw std::runtime_error("stored type isn't integer");}
	// template <> double      Variant::get() {if(t == Type::Double) {f = true; return u.d;} throw std::runtime_error("stored type isn't double" );}//can get int as double
	// template <> std::string Variant::get() {if(t == Type::String) {f = true; return   s;} throw std::runtime_error("stored type isn't string" );}

	//@brief  : set the stored value and update the type if needed
	//@param v: value to store
	template <> void Variant::set(const bool       & v) {u.b = v; t = Type::Bool  ; f = false;}
	template <> void Variant::set(const int        & v) {u.i = v; t = Type::Int   ; f = false;}
	template <> void Variant::set(const double     & v) {u.d = v; t = Type::Double; f = false;}
	template <> void Variant::set(const std::string& v) {  s = v; t = Type::String; f = false;}
	template <typename T> void Variant::set(const T& v) {
		static_assert(std::is_same<T, bool       >::value ||
		              std::is_same<T, int        >::value ||
		              std::is_same<T, double     >::value ||
		              std::is_same<T, std::string>::value, "T must be bool, int, double, or string");
		throw std::logic_error("Variant::set() shouldn't compile for this type");
	}

	//@brief : convert the currently held value to a string representatino
	//@return: string representation
	std::string Variant::to_string() const {
		if(t == Type::String) return s;
		std::stringstream ss;
		if(t == Type::Bool  ) ss << (u.b ? ".TRUE." : ".FALSE.");//convert bool   to string
		if(t == Type::Int   ) ss <<  u.i                        ;//convert int    to string
		if(t == Type::Double) ss <<  u.d                        ;//convert double to string
		return ss.str();
	}

	//@brief :check if the current value has been read
	//@return: true if the current value was read since setting, false otherwise
	bool Value::wasUsed() const {
		return std::all_of(this->cbegin(), this->cend(), [](const Variant& v){return v.used();});
	}

	namespace detail {
		//@brief: check if a string can be converted to the templated type
		//@param str: string to check
		//@param val: location to write parsed value
		//@return   : true if a value was successfully parsed, false otherwise
		template<typename T>
		bool tryParse(std::string str, T& val) {
			std::istringstream iss(str);//wrap string with istringstream
			iss >> std::noskipws >> val;//try to parse as 'T' without skipping white space
			return iss.eof() && !iss.fail();//make sure there wasn't an error and we consumed the entire string
		}
	}

	//@brief  :convert a string to lower case
	//@param s: string to convert
	//@return : lowercase string
	std::string NameList::ToLower(std::string s) {
		std::transform(s.begin(), s.end(), s.begin(), [](char c){return std::tolower(c);});
		return s;
	}

	//@brief   : read a namelist file
	//@param is: istream to read namelist from
	//@param cm: comment characters (ignore lines if the first character after white space is in cm)
	void NameList::read(std::istream& is, std::string cm) {
		this->clear();
		char space, delim;
		std::string line, key, tok;
		size_t lineNum = 0;
		std::getline(is, line);//for now just skip the first line
		fileLines.push_back(line);
		if(line.empty() ? false : std::string::npos != line.find('=')) throw std::runtime_error("namelist files cannot have key value pairs in the first line");
		++lineNum;
		bool gotComma = true;//was there a comma on the previous line
		while(std::getline(is, line)) {//loop over lines while the input stream is good
			//handle empty and comment lines
			++lineNum;
			fileLines.push_back(line);
			if(line.empty()) continue;//skip empty lines
			if(std::string::npos != cm.find_first_of(line.front())) continue;//skip comment lines
			if(0 == line.compare(" /")) continue;//final line

			//wrap line as input string stream
			std::istringstream iss(line);

			//make sure there was a comma before this token (or it was the first token) and check if there is one afterwards
			if(!gotComma) throw std::runtime_error("missing comma between previous entry and namelist line " + std::to_string(lineNum) + " \"" + line + "\"");
			gotComma = false;//clear flag in case we reverse all the way through the line
			iss.seekg(-1,std::ios::end);//move to end of line
			const std::streamoff len = iss.tellg();//get length
			for(std::streamoff i = 0; i <= len; i++) {//reverse loop over iss
				if(std::isspace(iss.peek())) {//keep reversing through whitespace
					iss.unget();
				} else {//this is an actual characater
					gotComma = ',' == iss.peek();//check if the last non-white space character is a comma
					iss.seekg(0);//return to begining
					iss.clear();//clear eof bit
					break;//we're done
				}
			}

			//extract " key = " and make sure it isn't a duplicate
			if(!(iss >> std::noskipws >> space >> key >> std::skipws >> delim >> std::noskipws)) throw std::runtime_error("error parsing line '" + line + "' from name list");
			if(space != ' ') throw std::runtime_error("missing leading space in namelist line " + std::to_string(lineNum) + " \"" + line + "\"");
			key = ToLower(key);
			if(this->find(key) != this->end()) throw std::runtime_error("key \"" + key + "\" was defined twice in the name list");
			if(delim != '=') throw std::runtime_error("bad delimeter (expected '=') in namelist line " + std::to_string(lineNum) + " \"" + line + "\"");
			while(std::isspace(iss.peek()) && iss) iss.get();//skip all whitespce after '=' before value

			//now extract comma delimited tokens from the reaminging input stream
			//in the future (c++14/17) this probably be updated with std::quoted
			Value val;
			if('\'' == iss.peek()) {//we have strings to extract
				iss.get();//skip '
				while(iss.good()) {
					//traverse looking for next unescaped '
					char c;
					tok.clear();//create an empty string to accumlate into
					bool found = false;
					while(iss >> c) {//extract one character from the input stream
						if('\'' == c) {//we may have found a closing (unescaped quote)
							if(tok.empty() ? true : '\\' != tok.back()) {
								found = true;//closing found
								break;//we're done
							} else {
								tok.pop_back();//remove escapement
							}
						}
						tok += c;//accumulate characters into token
					}

					//make sure we found the end of the string and save
					if(!found) throw std::runtime_error("no closing quote for a token in line " + std::to_string(lineNum) + " \"" + line + "\"");
					val.push_back(tok);//save token

					//advance to the start of the next string
					if(iss >> std::skipws >> delim) {
						if(',' != delim) throw std::runtime_error("unexpected delimiter between strings in line " + std::to_string(lineNum) + " \"" + line + "\"");
						if(iss >> std::skipws >> delim) {
							if('\'' != delim) throw std::runtime_error("unexpected delimiter string opening in line " + std::to_string(lineNum) + " \"" + line + "\"");
						}
					}
				}
			} else {//we have string representations of booleans, integers, or floating point numbers to extract
				while(std::getline(iss, tok, ',')) {//tokenize
					//clean up token string
					tok.erase(std::remove_if(tok.begin(), tok.end(), [](const char& c){return std::isspace(c);}), tok.end());//trim white space
					std::transform(tok.begin(), tok.end(), tok.begin(), [](char c){return std::tolower(c);});//convert to lower case

					//attempt to parse
					if(!tok.empty()) {
						int vI;
						double vD;
						if(0 == tok.compare(".true.")) {//first check fortran style logical
							val.push_back(true);
						} else if(0 == tok.compare(".false.")) {
							val.push_back(false);
						} else if(detail::tryParse<double>(tok, vD)) {//if it isn't a logical it must be a number
							if(detail::tryParse<int>(tok, vI)) {//strings parse-able as int are a subset of string parse-able as fp
								val.push_back(vI);
							} else {
								val.push_back(vD);
							}
						} else {
							throw std::runtime_error("couldn't parse token \"" + tok + "\" from line " + std::to_string(lineNum) + " \"" + line + "\" as bool, int, or float (strings must be in single quotes, e.g. key = 'value')");
						}
					}
				}

				//check for a mix of parsed types
				bool hasBool = false, hasInt = false, hasDouble = false;
				for(const Variant& v : val) {
					switch(v.getType()) {
						case Type::Bool  : hasBool   = true; break;
						case Type::Int   : hasInt    = true; break;
						case Type::Double: hasDouble = true; break;
						case Type::None  : //intentional fall through
						case Type::String: throw std::logic_error("parsed bool/int/double to string or none");
					}
				}

				//check for mixed types
				if((hasInt || hasDouble) && hasBool) throw std::runtime_error("line " + std::to_string(lineNum) + " \"" + line + "\" has a mix of numbers and booleans");

				//int + double ==> double
				if(hasInt && hasDouble) {//!hasBool
					for(Variant& v : val) {
						if(Type::Int == v.getType()) {//cast this int to a double for consistency
							v.set<double>(v.getInt());
						}
					}
				}
			}

			//save key/value pair
			// if(1 != val.size()) throw std::runtime_error("only scalar inputs are currently supported");
			this->operator[](key) = val;
		}
	}

	//@brief   : read a namelist file
	//@param nm: filename of namelist to read
	//@param cm: comment characters (ignore lines if the first character after white space is in cm)
	void NameList::read(std::string nm, std::string cm) {
		std::ifstream is(nm);
		read(is, cm);
	}

	//@brief    : find the first token with a given key
	//@param nm : key to search for
	//@return   : token with specified key (throws if not found)
	//@note     : write acces isn't public facing
	Value      & NameList::at(const std::string nm)       {
		try {
			return std::map<std::string, Value>::at(ToLower(nm));
		} catch (std::out_of_range&) {
			throw std::runtime_error("couldn't find `" + nm + "' in namelist");
		}
	}

	//@brief    : find the first token with a given key
	//@param nm : key to search for
	//@return   : token with specified key (throws if not found)
	Value const& NameList::at(const std::string nm) const {
		try {
			return std::map<std::string, Value>::at(ToLower(nm));
		} catch (std::out_of_range&) {
			throw std::runtime_error("couldn't find `" + nm + "' in namelist");
		}
	}

	//@brief    : attempt to parse a list fo bool from a list of namelist parameters
	//@param nm : key to search for
	//@return   : parsed bools (throws if not found)
	std::vector<bool       > NameList::getBools  (std::string nm) {
		std::vector<bool       > ret;
		Value& key = at(nm);
		for(size_t i = 0; i < key.size(); i++) ret.push_back(key[i].getBool  ());
		return ret;

	}

	//@brief    : attempt to parse a list fo int from a list of namelist parameters
	//@param nm : key to search for
	//@return   : parsed ints (throws if not found)
	std::vector<int        > NameList::getInts   (std::string nm) {
		std::vector<int        > ret;
		Value& key = at(nm);
		for(size_t i = 0; i < key.size(); i++) ret.push_back(key[i].getInt   ());
		return ret;

	}

	//@brief    : attempt to parse a list fo double from a list of namelist parameters
	//@param nm : key to search for
	//@return   : parsed doubles (throws if not found)
	std::vector<double     > NameList::getDoubles(std::string nm) {
		std::vector<double     > ret;
		Value& key = at(nm);
		for(size_t i = 0; i < key.size(); i++) ret.push_back(key[i].getDouble());
		return ret;

	}

	//@brief    : attempt to parse a list fo string from a list of namelist parameters
	//@param nm : key to search for
	//@return   : parsed strings (throws if not found)
	std::vector<std::string> NameList::getStrings(std::string nm) {
		std::vector<std::string> ret;
		Value& key = at(nm);
		for(size_t i = 0; i < key.size(); i++) ret.push_back(key[i].getString());
		return ret;

	}

	//@brief : check if there are unused tokens
	//@return: false if all tokens have been parse (via get___), false otherwise
	bool NameList::fullyParsed() const {
		return std::all_of(this->cbegin(), this->cend(), [](const std::pair<std::string, Value>& p){return p.second.wasUsed();});		
	}

	//@brief : get a list of unused tokens
	//@return: names of all unused tokens
	std::string NameList::unusedTokens() const {
		std::string unused;
		for(std::map<std::string, Value>::const_iterator iter = this->cbegin(); iter != this->cend(); ++iter) {//loop over key/value pairs
			if(!iter->second.wasUsed()) {//check if this value was used
				if(!unused.empty()) unused += ',';//seperate tokens
				unused += iter->first;//add token key
			}
		}
		return unused;
	}

#ifdef NML_USE_H5
	//@brief    : write all parsed tokens to an hdf file
	//@param grp: hgf group to write values to
	void NameList::writeParameters(H5::Group grp) {
		//loop over key/value pairs writing values of anything that was actually parsed
		for(std::map<std::string, Value>::const_iterator iter = this->cbegin(); iter != this->cend(); ++iter) {
			if(iter->second.wasUsed()) {//check if this value was used
				switch(iter->second.getType()) {
					case Type::None  : throw std::logic_error("cannot write None typed variant to HDF file");
					
					case Type::Bool  : {
						std::vector<bool       > vBool = getBools  (iter->first);
						std::vector<hbool_t    > vals(vBool.begin(), vBool.end());
						hsize_t dims[1] = {vals.size()};
						grp.createDataSet(iter->first, H5::PredType::NATIVE_HBOOL   , H5::DataSpace(1, dims)).write(vals.data(), H5::PredType::NATIVE_HBOOL  );
					} break;

					case Type::Int   : {
						std::vector<int        > vals = getInts   (iter->first);
						hsize_t dims[1] = {vals.size()};
						grp.createDataSet(iter->first, H5::PredType::NATIVE_INT     , H5::DataSpace(1, dims)).write(vals.data(), H5::PredType::NATIVE_INT    );
					} break;

					case Type::Double: {
						std::vector<double     > vals = getDoubles(iter->first);
						hsize_t dims[1] = {vals.size()};
						grp.createDataSet(iter->first, H5::PredType::NATIVE_DOUBLE , H5::DataSpace(1, dims)).write(vals.data(), H5::PredType::NATIVE_DOUBLE );
					} break;

					case Type::String: {
						std::vector<std::string> vStr = getStrings(iter->first);
						std::vector<char const*> vals;
						for(const std::string& s : vStr) vals.push_back(s.data());
						hsize_t dims[1] = {vals.size()};
						grp.createDataSet(iter->first, H5::StrType(0, H5T_VARIABLE), H5::DataSpace(1, dims)).write(vals.data(), H5::StrType(0, H5T_VARIABLE));
					} break;
				}
			}
		}

		//write an attribute with any unused values
		std::string unused = unusedTokens();
		if(!unused.empty()) grp.createAttribute("Unused Tokens", H5::StrType(0, H5T_VARIABLE), H5::DataSpace(H5S_SCALAR)).write(H5::StrType(0, H5T_VARIABLE), unused);
	}

	//@brief    : write the raw file to a dataset
	//@param grp: hgf group to create dataset in
	//@param nm : name of dataset to write
	void NameList::writeFile(H5::Group grp, std::string nm) const {
		std::vector<char const*> vals;
		for(const std::string& s : fileLines) vals.push_back(s.data());
		hsize_t dims[1] = {vals.size()};
		grp.createDataSet(nm, H5::StrType(0, H5T_VARIABLE), H5::DataSpace(1, dims)).write(vals.data(), H5::StrType(0, H5T_VARIABLE));
	}

#endif

}

#endif//_NML_H_
