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

#include <wx/app.h>

#include "wx/IndexingFrame.h"

class MyApp: public wxApp {
	IndexingFrame* frame;

	bool OnInit();

	public:
		DECLARE_EVENT_TABLE()

		//@brief: event handle for OSX about menu
		void OnAbout(wxCommandEvent& evt) {frame->showAbout();}
};

IMPLEMENT_APP(MyApp)

BEGIN_EVENT_TABLE(MyApp, wxApp)
	EVT_MENU(wxID_ABOUT, MyApp::OnAbout)
END_EVENT_TABLE()

#include "sphinx.xpm"

bool MyApp::OnInit() {
	frame = new IndexingFrame(NULL);

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
	static const bool setIcon = true;
#elif __APPLE__ || __linux__ || __unix__ || defined(_POSIX_VERSION)
	#if __APPLE__
		static const bool setIcon = false;//already handled by  bundle
	#else
		static const bool setIcon = true;
	#endif
#endif
	if(setIcon) frame->SetIcon( wxICON(sphinx) );
	
	frame->Show();

	return true;
}

