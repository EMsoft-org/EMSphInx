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

#ifndef _IDX_FRAME_H
#define _IDX_FRAME_H

#include <wx/artprov.h>
#include <wx/xrc/xmlres.h>
#include <wx/string.h>
#include <wx/bitmap.h>
#include <wx/image.h>
#include <wx/icon.h>
#include <wx/menu.h>
#include <wx/gdicmn.h>
#include <wx/font.h>
#include <wx/colour.h>
#include <wx/settings.h>
#include <wx/sizer.h>
#include <wx/gauge.h>
#include <wx/button.h>
#include <wx/statusbr.h>
#include <wx/frame.h>
#include <wx/splitter.h>
#include <wx/progdlg.h>
#include <wx/msgdlg.h>

#include "ImagePanel.h"
#include "EbsdSummaryPanel.h"
#include "EbsdNamelistWizard.h"
#include "BibtexDialog.h"

#include "modality/ebsd/nml.hpp"
#include "modality/ebsd/idx.hpp"
#include "util/timer.hpp"

///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
/// Class IndexingFrame
///////////////////////////////////////////////////////////////////////////////

BEGIN_DECLARE_EVENT_TYPES()
	DECLARE_EVENT_TYPE(wxEVT_IndexingThread, -1)//type for thread to communicate with frame
END_DECLARE_EVENT_TYPES()

class IndexingFrame : public wxFrame {
	private:
		std::thread                  idxThd ;//thread for indexing
		std::atomic_flag             flg    ;//flag to tell worker thread to stop
		emsphinx::ebsd::Namelist     nml    ;//namelist
		wxImage                      imOrig ;//original image (e.g. iq)
		wxImage                      imPaint;//ipf painted image

		enum IdxStatus {
			IdxError = -  1,
			IdxExit  = -  2,
			IdxCmplt =  101,
		};

	protected:
		wxMenuBar       * m_menubar ;
		wxMenu          * m_menuFile;
		wxMenu          * m_menuHelp;
		wxSplitterWindow* m_split   ;
		wxImagePanel    * m_imPan   ;
		EbsdSummaryPanel* m_sumPan  ;
		wxGauge         * m_prog    ;
		wxButton        * m_btn     ;
		wxStatusBar     * m_statBar ;

		//@brief: menu item functions
		virtual void OnFileOpen  ( wxCommandEvent& event ) { loadNamelist(); }
		virtual void OnFileLoad  ( wxCommandEvent& event ) { saveNamelist(); }
		virtual void OnFileWizard( wxCommandEvent& event ) { runWizard   (); }
		virtual void OnHelpAbout ( wxCommandEvent& event ) { showAbout   (); }
		virtual void OnHelpRefs  ( wxCommandEvent& event ) { showRefs(true); }
		virtual void OnHelpHelp  ( wxCommandEvent& event ) { showHelp    (); }

		void WizardClosed(wxCloseEvent& event);

		//@brief: action for start/stop button click
		virtual void OnBtn( wxCommandEvent& event ) {
			if("Start" == m_btn->GetLabel()) startIdx();
			else stopIdx();
		}

		//@brief: action for indexing thread updates
		void OnThread(wxCommandEvent& event) {
			//update the status bar and get event type
			m_statBar->SetStatusText(event.GetString());
			int i = event.GetInt();

			//handle event type
			if(IdxError == i) {//an error
				//nothing else to do for now (wait for exit message)
				wxMessageDialog msgDlg(this, event.GetString(), "Error Indexing");
				msgDlg.ShowModal();//make sure user sees message (since it will be overwritten with "Indexing Stopped")
			} else if(IdxCmplt == i || IdxExit == i) {//(un)successful indexing completion
				idxThd.join();//wait for the thread to wrap up
				m_sumPan->EnableEditing(true);//enable the parameter panel
				m_btn->SetLabel("Start");//change stop button to a start button again
				m_btn->Enable(true);\
				if(IdxCmplt == i) showRefs(false);
			} else {//status update
				m_btn->Enable(true);
				m_prog->SetValue(event.GetInt());//update progress bar
				m_imPan->setImage(imPaint);//update image
				m_imPan->Refresh();//repaint
				m_imPan->paintNow();
			}
		}

	public:
		//@brief: load settings from a namelist file
		void loadNamelist();

		//@brief: save the current settings to a namelist file
		void saveNamelist();

		//@brief: run the indexing wizard
		void runWizard();

		//@brief: run the about window
		void showAbout();

		//@brief      : run the references
		//@param force: should the references be shown even if the dialog was previously silenced
		void showRefs(const bool force);

		//@brief: show the program help
		void showHelp();

		void startIdx();

		void stopIdx();

		void setImage(wxImage im) {m_imPan->setImage(im); imOrig = im; imPaint = im.Copy(); m_imPan->Refresh();}

		IndexingFrame( wxWindow* parent, wxWindowID id = wxID_ANY, const wxString& title = "EMSphInx", const wxPoint& pos = wxDefaultPosition, const wxSize& size = wxSize( 720,540 ), long style = wxDEFAULT_FRAME_STYLE|wxTAB_TRAVERSAL );

		~IndexingFrame();

		wxDECLARE_EVENT_TABLE();
};

BEGIN_EVENT_TABLE(IndexingFrame, wxFrame)
    EVT_COMMAND(wxID_ANY, wxEVT_IndexingThread, IndexingFrame::OnThread)
END_EVENT_TABLE()

///////////////////////////////////////////////////////////////////////////
// C++ code generated with wxFormBuilder (version Oct 26 2018)
// http://www.wxformbuilder.org/
//
// PLEASE DO *NOT* EDIT THIS FILE!
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////

#include <wx/filedlg.h>
#include <wx/aboutdlg.h>
#include "constants.hpp"
#include "sht_file.hpp"

DEFINE_EVENT_TYPE(wxEVT_IndexingThread)

//@brief: load settings from a namelist file
void IndexingFrame::loadNamelist() {
	wxFileDialog opDlg(this, _("Load Namelist file"), "", "", "Namelist files (*.nml)|*.nml", wxFD_OPEN|wxFD_FILE_MUST_EXIST);
	if(opDlg.ShowModal() == wxID_CANCEL) return;//cancel

	//read nml and parse
	try {
		std::ifstream is(opDlg.GetPath().ToStdString());//open file
		std::string str((std::istreambuf_iterator<char>(is)), std::istreambuf_iterator<char>());//read entire thing into string
		std::string warning = nml.from_string(str);//parse namelist
		m_sumPan->setNamelist(nml);
		if(!warning.empty()) {//display any warnings
			wxMessageDialog msgDlg(this, "unused parameters: " + warning, "Warnings Parsing Namelist");
			msgDlg.ShowModal();//make sure user sees message (since it will be overwritten with "Indexing Stopped")
		}
	} catch (std::exception& e) {
		nml.clear();
		wxMessageDialog msgDlg(this, e.what(), "Error Reading Namelist");
		msgDlg.ShowModal();//make sure user sees message (since it will be overwritten with "Indexing Stopped")
	}
}

//@brief: save the current settings to a namelist file
void IndexingFrame::saveNamelist() {
	wxFileDialog svDlg(this, _("Save Namelist file"), "", "", "Namelist files (*.nml)|*.nml", wxFD_SAVE|wxFD_OVERWRITE_PROMPT);
	if(svDlg.ShowModal() == wxID_CANCEL) return;//cancel
	std::ofstream os(svDlg.GetPath().ToStdString());
	os << nml.to_string();
}

//@brief: run the indexing wizard
void IndexingFrame::runWizard() {
	EbsdNamelistWizard* wizard = new EbsdNamelistWizard(this);
	wizard->setNamelist(nml);
	wizard->Show();
	wizard->Connect( wxEVT_CLOSE_WINDOW, wxCloseEventHandler( IndexingFrame::WizardClosed ), NULL, this );
}

void IndexingFrame::WizardClosed(wxCloseEvent& event) {
	EbsdNamelistWizard* wizard = (EbsdNamelistWizard*) event.GetEventObject();
	if(wizard->isFinished()) {
		m_sumPan->setNamelist(wizard->getNamelist());
		nml = wizard->getNamelist();
		wxImage im = wizard->getMap();
		setImage(im);
	} else {
	}
	event.Skip();
}

//@brief: run the about window
void IndexingFrame::showAbout() {
	wxAboutDialogInfo aboutInfo;
	aboutInfo.SetName("EMSphInx");
	aboutInfo.SetVersion(emsphinx::GitBranch + "-" + emsphinx::GitHash);//the ":" in Version seems to cause display issues...
	aboutInfo.SetDescription(_("Spherical Harmonics EBSD Indexing"));
	aboutInfo.SetCopyright("\
EMSphInx: (C) 2019 De Graef Group, Carnegie Mellon University\n\
crystallography: (C) 2015-2019 William C. Lenthe\n\
utilities: (C) 2018-2019 William C. Lenthe\n\
wxWidgets: (C) 1998-2005 Julian Smart, Robert Roebling et al\n\
FFTW: (C) 2003, 2007-11 Matteo Frigo, Massachusetts Institute of Technology\n\
miniz: (C) 2013-2014 RAD Game Tools and Valve Software\n\
miniz: (C) 2010-2014 Rich Geldreich and Tenacious Software LLC\n\
");
	aboutInfo.SetWebSite("https://github.com/EMsoft-org/Emsphinx");
	aboutInfo.AddDeveloper("William C. Lenthe");
	aboutInfo.AddArtist("William C. Lenthe");
	aboutInfo.SetLicense("\
GPL-2.0 (https://www.gnu.org/licenses/old-licenses/gpl-2.0.html)\n\
with BSD-3 components, see the license file for details\n\
Interested in a commercial license? Contact:\n\
\n\
Center for Technology Transfer and Enterprise Creation\n\
4615 Forbes Avenue, Suite 302\n\
Pittsburgh, PA 15213\n\
\n\
phone: 412.268.7393\n\
email: innovation@cmu.edu\n\
website: https://www.cmu.edu/cttec/");
    wxAboutBox(aboutInfo);
}

//@brief      : run the references
//@param force: should the references be shown even if the dialog was previously silenced
void IndexingFrame::showRefs(const bool force) {
	const std::string silFil = getUserAppDataDir() + "silenced";//flag file for silencing
	const bool silenced = wxFileExists(silFil);
	if(force || !silenced) {
		//start with references from software
		std::string refs = "indexing:\n" + emsphinx::references::LentheIdx
		                 + "\nmaster pattern file format:\n" + emsphinx::references::LentheFile
		                 + "\npsuedo-symmetry:\n" + emsphinx::references::LenthePS;

		//accumulate references from selected master pattern files
		wxArrayString mpList = m_sumPan->getMP();
		std::string acmDoi, acmXtal;
		for(size_t i = 0; i < mpList.GetCount(); i++) {
			std::ifstream is(mpList[i].ToStdString(), std::ios::in | std::ios::binary);//open file for reading
			if(is) {//file actually exists
				//read file and get file level doi
				sht::File file;
				file.read(is);
				std::string doi = file.header.doi.substr(0, file.header.doiLen());
				if(!doi.empty()) acmDoi += doi + '\n';//accumulate doi

				//now loop over crystals getting structure references
				for(const sht::CrystalData& xtal : file.mpData.xtals) {
					std::string xref = xtal.refs.substr(0, xtal.refsLen());
					if(!xref.empty()) acmXtal += xref + '\n';//accumulate doi
				}
			}
		}

		//add file SHT file references
		if(!acmDoi .empty()) refs += "\nmaster patterns:\n" + acmDoi;
		if(!acmXtal.empty()) refs += "\ncrystal structures:\n" + acmXtal;

		//show references
		BibtexDialog bDlg(refs, this);
		bDlg.Silence(silenced);
		bDlg.ShowModal();
		if(bDlg.Silence() && !silenced) {// unsilenced --> silenced
			std::ofstream os(silFil);
		} else if(!bDlg.Silence() && silenced) {// silenced --> unsilenced
			remove(silFil.c_str());
		}
	}
}

//@brief: show the program help
void IndexingFrame::showHelp() {
	if(!wxLaunchDefaultBrowser("https://emsphinx.readthedocs.io/")) {
		wxMessageDialog msgDlg(this, "Please visit https://emsphinx.readthedocs.io/ for documentation", "Error Launching Browser");
	}
}

void IndexingFrame::startIdx() {
	//update namelist
	m_sumPan->EnableEditing(false);//disable the parameter panel
	if(m_sumPan->updateNamelist(&nml)) nml.namelistFile.clear();

	//make sure we have a scan (most parameters are sanity checked during indexing, but rescaling an empty image throws)
	if(nml.scanDims[0] < 1 || nml.scanDims[1] < 1) {
		wxMessageDialog msgDlg(this, "Empty scan (0 pixels)", "Error Indexing");
		msgDlg.ShowModal();//make sure user sees message (since it will be overwritten with "Indexing Stopped")
		m_sumPan->EnableEditing(true);//disable the parameter panel
		return;
	}

	//make sure the image is the right size
	bool ok = imOrig.IsOk();
	if( !imOrig.IsOk() || imOrig.GetWidth() != nml.scanDims[0] || imOrig.GetHeight() != nml.scanDims[1]) {
		imOrig = wxImage (nml.scanDims[0], nml.scanDims[1]);
		imOrig.SetAlpha();
		unsigned char* pAlpha = imOrig.GetAlpha();
		std::fill(pAlpha, pAlpha + nml.scanDims[0] * nml.scanDims[1], 0xFF);
		setImage(imOrig);
	}

	//update gui
	imPaint = imOrig.Copy();
	m_btn->SetLabel("Stop");//change start button to stop button
	m_btn->Enable(false);

	//start indexing thread
	flg.test_and_set();
	idxThd = std::thread([&](){
		//create event for dispatching
		wxCommandEvent evt(wxEVT_IndexingThread);
		evt.SetInt(0);
		evt.SetString("Initializing Indexers");
		wxPostEvent(this, evt);//update image + progress bar
		try {
			emsphinx::ebsd::IndexingData<double> idxData(nml, (char*)imPaint.GetData());//read inputs and build indexers
			ThreadPool pool(idxData.threadCount);//build thread pool
			try {
				//queue parallel indexing
				time_t tmStart = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
				Timer t;
				size_t batches = idxData.pat->numPat() / nml.batchSize;//how many batches are needed
				if(batches * nml.batchSize < idxData.pat->numPat()) ++batches;//extra batch for leftovers
				for(size_t i = 0; i < batches; i++) {//loop over batches
					const size_t start = i * nml.batchSize;//first pattern
					const size_t end = std::min(start + nml.batchSize, idxData.pat->numPat());//last pattern
					pool.schedule(std::bind(idxData.workItem, std::placeholders::_1));//queue indexing
				}

				//wait for completion posting updates periodically
				const std::chrono::milliseconds uptFreq(250);//milliseconds between updates (too long -> unresponsive, too short -> resource intensive), 0.25 ~human reaction time
				while(!pool.waitAll(uptFreq)) {
					if(!flg.test_and_set()) {
						//post error status update
						evt.SetString("Stop Requested");
						wxPostEvent(this, evt);
						pool.clear();//clear all unstarted items
						evt.SetString("Indexing Stopped");
						evt.SetInt(0);
						wxPostEvent(this, evt);//update image + progress bar
						evt.SetInt(IdxExit);
						wxPostEvent(this, evt);//mark exit
						return;
					}

					//get the time elapsed since we started (without resetting the reference point)
					const double elapsed = t.poll(false);
					const double percent = double(idxData.idxCtr) / idxData.numIdx;
					const double rate = elapsed / percent;//estimate seconds per %
					const double remaining = rate * (1.0 - percent);//estimate seconds left

					//print update
					std::ostringstream ss;
					Timer::PrintSeconds(elapsed  , ss);
					ss << " elapsed, ";
					Timer::PrintSeconds(remaining, ss);
					ss << " remaining   ";
					evt.SetString(ss.str());
					evt.SetInt(std::min((int)(percent * 100 + 0.5), 100));
					wxPostEvent(this, evt);
				}
				
				//compute indexing time and save
				const double total = t.poll();
				time_t tmEnd = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
				idxData.save(tmStart, tmEnd, total);
			
				//post completion status
				evt.SetString("Indexing Complete");
				evt.SetInt(100);
				wxPostEvent(this, evt);//post status update with 100% completion
				evt.SetInt(IdxCmplt);
				wxPostEvent(this, evt);//post actual completion message
				return;
			 } catch (std::exception& e) {//clear items and post
				pool.clear();//clear all unstarted items
				//post error message
				evt.SetString(e.what());
				evt.SetInt(IdxError);
				wxPostEvent(this, evt);
			}
			//pool goes out of scope here (waiting for completion)
		} catch (std::exception& e) {
			//post error message
			evt.SetString(e.what());
			evt.SetInt(IdxError);
			wxPostEvent(this, evt);
		} catch (...) {
			evt.SetString("Unhandled exception type during indexing, please contact developers");
			evt.SetInt(IdxError);
			wxPostEvent(this, evt);
		}

		//post exit message
		evt.SetString("Indexing Stopped");
		evt.SetInt(0);
		wxPostEvent(this, evt);//update image + progress bar
		evt.SetInt(IdxExit);
		wxPostEvent(this, evt);//mark exit
	});
}

void IndexingFrame::stopIdx() {
	flg.clear();//tell any running indexers to stop
	m_sumPan->EnableEditing(true);//enable the parameter panel
	m_btn->Enable(false);//disable the start/stop button
}

IndexingFrame::IndexingFrame( wxWindow* parent, wxWindowID id, const wxString& title, const wxPoint& pos, const wxSize& size, long style ) : wxFrame( parent, id, title, pos, size, style ) {
	this->SetSizeHints( wxDefaultSize, wxDefaultSize );

	//build menus
	m_menubar  = new wxMenuBar( 0 );
	m_menuFile = new wxMenu();
	m_menuHelp = new wxMenu();

	//build menu items
	wxMenuItem* m_menuFileOpen   = new wxMenuItem( m_menuFile, wxID_ANY  , wxString( wxT("Open..."     ) ) + wxT('\t') + wxT("ctrl+o"), wxString( wxT("load a namelist file"  ) ), wxITEM_NORMAL );
	wxMenuItem* m_menuFileSaveAs = new wxMenuItem( m_menuFile, wxID_ANY  , wxString( wxT("Save As..."  ) ) + wxT('\t') + wxT("ctrl+s"), wxString( wxT("export a namelist file") ), wxITEM_NORMAL );
	wxMenuItem* m_menuFileWizard = new wxMenuItem( m_menuFile, wxID_ANY  , wxString( wxT("Wizard..."   ) ) + wxT('\t') + wxT("ctrl+w"), wxString( wxT("launch nml builder"    ) ), wxITEM_NORMAL );
	wxMenuItem* m_menuHelpAbout  = new wxMenuItem( m_menuHelp, wxID_ABOUT, wxString( wxT("About"    ) )                            , wxString( wxT("about this software"   ) ), wxITEM_NORMAL );
	wxMenuItem* m_menuHelpRefs   = new wxMenuItem( m_menuHelp, wxID_ANY  , wxString( wxT("Citations...") )                            , wxString( wxT("relevant literature"   ) ), wxITEM_NORMAL );
	wxMenuItem* m_menuHelpHelp   = new wxMenuItem( m_menuHelp, wxID_HELP , wxString( wxT("Help..."     ) )                            , wxString( wxT("documentation browser" ) ), wxITEM_NORMAL );

	//set menu item bitmaps
	m_menuFileOpen  ->SetBitmap( wxArtProvider::GetBitmap( wxART_FILE_OPEN   , wxART_MENU ) );
	m_menuFileSaveAs->SetBitmap( wxArtProvider::GetBitmap( wxART_FILE_SAVE_AS, wxART_MENU ) );
	m_menuFileWizard->SetBitmap( wxNullBitmap                                               );
	m_menuHelpAbout ->SetBitmap( wxArtProvider::GetBitmap( wxART_INFORMATION , wxART_MENU ) );
	m_menuHelpRefs  ->SetBitmap( wxNullBitmap                                               );
	m_menuHelpHelp  ->SetBitmap( wxArtProvider::GetBitmap( wxART_HELP_BOOK   , wxART_MENU ) );

	//assemble file menu
	m_menuFile->Append( m_menuFileOpen   );
	m_menuFile->Append( m_menuFileSaveAs );
	m_menuFile->Append( m_menuFileWizard );
	m_menubar ->Append( m_menuFile, wxT("File") );

	//assemble help menu
	m_menuHelp->Append( m_menuHelpAbout );
	m_menuHelp->Append( m_menuHelpRefs  );
	m_menuHelp->Append( m_menuHelpHelp  );
	m_menubar ->Append( m_menuHelp, wxT("Help") );

	this->SetMenuBar( m_menubar );

	//create splitter window
	m_split = new wxSplitterWindow( this, wxID_ANY, wxDefaultPosition, wxDefaultSize, wxSP_3DSASH|wxSP_LIVE_UPDATE );
	m_split->SetSashGravity( 1 );//put all growth into left side
	m_split->SetMinimumPaneSize( 50 );

	//create image + indexing parameters and add to splitter
	m_imPan  = new wxImagePanel    ( m_split, wxID_ANY,              wxDefaultPosition, wxDefaultSize, wxTAB_TRAVERSAL );
	m_sumPan = new EbsdSummaryPanel( m_split, wxID_ANY,              wxDefaultPosition, wxDefaultSize, wxTAB_TRAVERSAL );
	m_split->SplitVertically(m_imPan, m_sumPan);

	//create sizers and assemble
	wxBoxSizer* vSizer = new wxBoxSizer( wxVERTICAL   );
	wxBoxSizer* hSizer = new wxBoxSizer( wxHORIZONTAL );
	vSizer->Add( m_split, 1, wxEXPAND, 5 );
	vSizer->Add( hSizer , 0, wxEXPAND, 5 );

	//build and assemble controls
	m_prog   = new wxGauge         ( this, wxID_ANY, 100        , wxDefaultPosition, wxDefaultSize, wxGA_HORIZONTAL );
	m_btn    = new wxButton        ( this, wxID_ANY, wxT("Start"), wxDefaultPosition, wxDefaultSize, 0 );
	m_prog->SetValue( 0 );
	m_prog->SetRange( 100 );
	hSizer->Add( m_prog  , 1, wxEXPAND, 0 );
	hSizer->Add( m_btn   , 0, wxEXPAND, 0 );

	//layout
	this->SetSizer( vSizer );
	this->Layout();
	m_statBar = this->CreateStatusBar( 1, wxSTB_ELLIPSIZE_END|wxSTB_SIZEGRIP, wxID_ANY );
	this->Centre( wxBOTH );

	// Connect Events
	m_btn->Connect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( IndexingFrame::OnBtn ), NULL, this );
	m_menuFile->Bind(wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler( IndexingFrame::OnFileOpen   ), this, m_menuFileOpen  ->GetId());
	m_menuFile->Bind(wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler( IndexingFrame::OnFileLoad   ), this, m_menuFileSaveAs->GetId());
	m_menuFile->Bind(wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler( IndexingFrame::OnFileWizard ), this, m_menuFileWizard->GetId());
	m_menuHelp->Bind(wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler( IndexingFrame::OnHelpAbout  ), this, wxID_ABOUT);
	m_menuHelp->Bind(wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler( IndexingFrame::OnHelpRefs   ), this, m_menuHelpRefs  ->GetId());
	m_menuHelp->Bind(wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler( IndexingFrame::OnHelpHelp   ), this, wxID_HELP);
}

IndexingFrame::~IndexingFrame() {
	m_btn->Disconnect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( IndexingFrame::OnBtn ), NULL, this );
}

#endif//_IDX_FRAME_H

