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

#ifndef _SYSNAMES_H_
#define _SYSNAMES_H_

#include <string>

//@brief : get the name of the logged in user
//@return: user name
std::string getUserName();

//@brief : get the name of the computer (often the host name)
//@return: computer name
std::string getComputerName();

//@brief    : create a folder if needed needed
//@param dir: directory to create
//@return   : full path to createed (or existing) folder with terminating '\' or '/'
std::string buildFolder(std::string dir);

//@brief : get the path separator character
//@return: '\\' for windows, '/' for linux
char getPathSep();

//@brief : get the the folder to write data shared by applications and users
//@return: windows - C:\ProgramData\ (need this text to silence compiler warnings about "\ \n")
//         apple   - ~/Library/Application Support/ (should be '/Library/Preferences' but that isn't useful without privileges)
//         linux   - ~/.local/share/                (should be '/etc' but that isn't useful without privileges)
std::string getSharedDataDir();

//@brief : get the the folder to write data shared by applications for current user
//@return: folder to write data shared across users
//@note  : windows - C:\users\{user}\AppData\Local on windows
//         apple - ~/Library/Application Support/
//         linux - ~/.local/share (should be /usr/lib but that isn't useful without privileges)
std::string getUserDataDir();

//@brief    : get the the folder to write data for current application shared by all users
//@param app: application name
//@return   : getSharedDataDir()/app
std::string getAppDataDir(std::string app = "EMSphInx");

//@brief    : get the the folder to write data for current application and current user
//@param app: application name
//@return   : getUserDataDir()/app
std::string getUserAppDataDir(std::string app = "EMSphInx");

//@brief    : get the the folder to write data for current application shared by all users
//@param app: application name
//@return   : getSharedDataDir()/app
std::string getAppDataDir(std::string app) {return buildFolder(getSharedDataDir() + app);}

//@brief    : get the the folder to write data for current application and current user
//@param app: application name
//@return   : getUserDataDir()/app
std::string getUserAppDataDir(std::string app) {return buildFolder(getUserDataDir() + app);}

//@brief     : check if a file exists
//@param name: file name to check existance of
//@return    : true if the file exists, false otherwise
bool fileExists(std::string name);

//@brief         : determine the size of a file in bytes
//@param filename: name of file to get size of
//@return        : size of file in bytes (-1 if the file doesn't exist)
std::int64_t fileSize(const char* name);
std::int64_t fileSize(std::string name) {return fileSize(name.c_str());}

////////////////////////////////////////////////////////////////////////////////
//                           Implementation Details                           //
////////////////////////////////////////////////////////////////////////////////

#include <stdexcept>

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
	//limit windows includes
	#ifndef NOMINMAX//don't redefine
		#define NOMINMAX//the min/max macros can cause issues with std::min/max
	#endif
	#ifndef WIN32_LEAN_AND_MEAN//don't redefine
		#define WIN32_LEAN_AND_MEAN//we don't need much, don't bother including it all the extra bits
	#endif

	//define windows version (for GetFileSizeEx)
	#ifndef WINVER//don't redefine
		#define WINVER 0x0501//this is windows xp so it shouldn't be asking for too much...
	#endif
	#include <windows.h>
	#include <Lmcons.h>
	#include <shlobj.h>

	//@brief : get the name of the logged in user
	//@return: user name
	std::string getUserName() {
		char uNm[UNLEN+1];
		DWORD szUser = sizeof(uNm);
		if(GetUserNameA(uNm, &szUser)) return std::string(uNm);//or GetUserNameEx 
		throw std::runtime_error("failed to GetUserName");
	}

	//@brief : get the name of the computer (often the host name)
	//@return: computer name
	std::string getComputerName() {
		char hNm[UNLEN+1];
		DWORD szHost = sizeof(hNm);
		if(GetComputerNameA(hNm, &szHost)) return std::string(hNm);
		throw std::runtime_error("failed to GetComputerName");
	}

	//@brief    : create the specified directory if needed and return name
	//@param dir: directory to create
	//@return   : full path to created (or existing) folder with terminating '\'
	std::string buildFolder(const std::string dir) {
		if(CreateDirectoryA(dir.c_str(), NULL) || ERROR_ALREADY_EXISTS == GetLastError()) return dir + '\\';//make sure the directory exists and return
		throw std::runtime_error("failed to CreateDirectory");
	}

	//@brief : get the path separator character
	//@return: '\\' for windows, '/' for linux
	char getPathSep() {return '\\';}

	//@brief : get the the folder to write data shared by applications and users
	//@return: windows - C:\ProgramData\ (need this text to silence compiler warnings about "\ \n")
	//         apple   - ~/Library/Application Support/ (should be '/Library/Preferences' but that isn't useful without privileges)
	//         linux   - ~/.local/share/                (should be '/etc' but that isn't useful without privileges)
	std::string getSharedDataDir() {
		char dir[MAX_PATH];
		if(SHGetFolderPathA(NULL, CSIDL_COMMON_APPDATA, NULL, 0, dir)) throw std::runtime_error("SHGetFolderPathA failed");//this is C:\ProgramData (analogous to /var)
		return buildFolder(dir);
	}

	//@brief : get the the folder to write data shared by applications for current user
	//@return: folder to write data shared across users
	//@note  : windows - C:\users\{user}\AppData\Local on windows
	//         apple - ~/Library/Application Support/
	//         linux - ~/.local/share (should be /usr/lib but that isn't useful without privileges)
	std::string getUserDataDir() {
		char dir[MAX_PATH];
		if(SHGetFolderPathA(NULL, CSIDL_LOCAL_APPDATA, NULL, 0, dir)) throw std::runtime_error("SHGetFolderPathA failed");//this is C:\users\{user}\AppData\Local (analogus to /usr)
		return buildFolder(std::string(dir));
	}

	//@brief     : check if a file exists
	//@param name: file name to check existance of
	//@return    : true if the file exists, false otherwise
	bool fileExists(std::string name) {
		const DWORD attrib = GetFileAttributesA(name.c_str());//get info about the file
		return (attrib != INVALID_FILE_ATTRIBUTES && !(attrib & FILE_ATTRIBUTE_DIRECTORY));//check if the file exists
	}

	//@brief     : determine the size of a file in bytes
	//@param name: name of file to get size of
	//@return    : size of file in bytes (-1 if the file doesn't exist)
	std::int64_t fileSize(const char* name) {
		//first make sure the file exists
		const DWORD attrib = GetFileAttributesA(name);//get info about the file
		const bool fileExists = (attrib != INVALID_FILE_ATTRIBUTES && !(attrib & FILE_ATTRIBUTE_DIRECTORY));//check if the file exists
		if(!fileExists) return -1;

		//if it does exist get the size
		WIN32_FILE_ATTRIBUTE_DATA fad;
		if(!GetFileAttributesExA(name, GetFileExInfoStandard, &fad)) throw std::runtime_error("failed to get file attributes");
		LARGE_INTEGER size;
		size.HighPart = fad.nFileSizeHigh;
		size.LowPart = fad.nFileSizeLow;
		return (std::int64_t) size.QuadPart;
	}

#elif __APPLE__ || __linux__ || __unix__ || defined(_POSIX_VERSION)

	#include <unistd.h>
	#include <pwd.h>
	#include <fcntl.h>
	#include <sys/types.h>
	#include <sys/param.h>
	#include <sys/stat.h>
	#include <stdlib.h>
	
	//___64 is now deprecated on mac but should still be preferred on linux
	#if __APPLE__
		#include <Availability.h>
		#if __MAC_OS_X_VERSION_MAX_ALLOWED >= 1050//1050 is __MAC_10_5, the first version that stat64/mmap64 are deprecated
			#define __64_FUNCS_DEPRECATED_
		#endif
	#endif
	#ifdef __64_FUNCS_DEPRECATED_//___64 functions are deprecated (e.g. just 'stat' is the 64 bit version and 32 bit versions don't exist)
		#define STAT  stat
	#else//explicitly use 64 bit functions
		#define STAT  stat64
	#endif

	//@brief : get the name of the logged in user
	//@return: user name
	std::string getUserName() {
		struct passwd *p = getpwuid(getuid());//try to get password info (or geteuid if you want the effective user, getpwuid_r for thread save)
		if(NULL == p) throw std::runtime_error("failed to getpwuid_r");
		std::string nm(p->pw_gecos);//get full user name
		return nm.empty() ? std::string(p->pw_name) : nm;//fall back to login name
	}

	//@brief : get the name of the computer (often the host name)
	//@return: computer name
	std::string getComputerName() {
		char hNm[MAXHOSTNAMELEN];
		if(gethostname(hNm, MAXHOSTNAMELEN)) throw std::runtime_error("failed to gethostname");
		return std::string(hNm);
	}

	//@brief    : create a folder if needed needed
	//@param dir: directory to create
	//@return   : full path to createed (or existing) folder with terminating '\' or '/'
	std::string buildFolder(std::string dir) {
		struct stat sb;
		if(stat(dir.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode)) return dir + '/';//already exists
		int ret = mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
		if(0 != ret) throw std::runtime_error("failed to create directory " + dir);
		return dir + '/';
	}

	//@brief : get the path separator character
	//@return: '\\' for windows, '/' for linux
	char getPathSep() {return '/';}

	//@brief : get the the folder to write data shared by applications and users
	//@return: windows - C:\ProgramData\ (need this text to silence compiler warnings about "\ \n")
	//         apple   - ~/Library/Application Support/ (should be '/Library/Preferences' but that isn't useful without privileges)
	//         linux   - ~/.local/share/                (should be '/etc' but that isn't useful without privileges)
	std::string getSharedDataDir() {
		//get the data home
		char* env = getenv("XDG_DATA_HOME");//first try to get $XDG_DATA_HOME
		if(NULL != env) return buildFolder(env);//got XDG_DATA_HOME
		env = getenv("HOME");//fall back to $HOME/.local/share
	#if __APPLE__
		const std::string par = "/Library";
		const std::string cld = "Application Support";
	#else
		const std::string par = "/.local";
		const std::string cld = "share";
	#endif 
		if(NULL != env) return buildFolder(buildFolder(std::string(env) + par) + cld);//got HOME

		//if there we no environment variables fall back to passwd
		struct passwd *p = getpwuid(getuid());//try to get password info (or geteuid if you want the effective user)
		if(NULL != p) return buildFolder(buildFolder(std::string(p->pw_dir) + par) + cld);//got passwd

		//if there was no passwd we've run out of options
		throw std::runtime_error("failed to getpwuid_r");
	}

	//@brief : get the the folder to write data shared by applications for current user
	//@return: folder to write data shared across users
	//@note  : windows - C:\users\{user}\AppData\Local on windows
	//         apple - ~/Library/Application Support/
	//         linux - ~/.local/share (should be /usr/lib but that isn't useful without privileges)
	std::string getUserDataDir() { return getSharedDataDir(); }
	
	//@brief     : check if a file exists
	//@param name: file name to check existance of
	//@return    : true if the file exists, false otherwise
	bool fileExists(std::string name) {
		struct stat fileStat;
		return stat(name.c_str(), &fileStat) >= 0;//try to get information about the file without opening it
	}

	//@brief     : determine the size of a file in bytes
	//@param name: name of file to get size of
	//@return    : size of file in bytes (-1 if the file doesn't exist)
	std::int64_t fileSize(const char* name) {
		struct STAT fileStat;
		const bool fileExists = STAT(name, &fileStat) >= 0;//try to get information about the file without opening it
		return fileExists ? fileStat.st_size : -1;
	}

#else

	static_assert(false, "user name and host name not implemented for this os");

#endif

#endif//_SYSNAMES_H_
