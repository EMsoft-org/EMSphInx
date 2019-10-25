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

#ifndef _MP_SEL_H_
#define _MP_SEL_H_

#include <wx/artprov.h>
#include <wx/xrc/xmlres.h>
#include <wx/string.h>
#include <wx/checklst.h>
#include <wx/gdicmn.h>
#include <wx/font.h>
#include <wx/colour.h>
#include <wx/settings.h>
#include <wx/bmpbuttn.h>
#include <wx/bitmap.h>
#include <wx/image.h>
#include <wx/icon.h>
#include <wx/button.h>
#include <wx/sizer.h>
#include <wx/statbox.h>
#include <wx/srchctrl.h>
#include <wx/panel.h>
#include <wx/rearrangectrl.h>

#include "wx/ValidityPanel.h"
#include "wx/MasterFileList.hpp"

///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
/// Class MasterPatternSelectPanel
///////////////////////////////////////////////////////////////////////////////
class MasterPatternSelectPanel : public ValidityPanel {

	MasterFileList * m_selLst  ;
	wxBitmapButton * m_btnUp   ;
	wxBitmapButton * m_btnDn   ;
	wxBitmapButton * m_btnAdd  ;
	wxBitmapButton * m_btnDel  ;
	wxBitmapButton * m_btnFlt  ;
	wxBitmapButton * m_btnDir  ;
	wxSearchCtrl   * m_mpSearch;
	MasterFileList * m_libLst  ;

	wxArrayString m_hidNames;//hidden patterns from known master pattern list (not matching current sort)

	void UpClicked      ( wxCommandEvent& event ) { m_selLst->MoveSelected(true );}//move the selected item in the top box up
	void DownClicked    ( wxCommandEvent& event ) { m_selLst->MoveSelected(false);}//move the selected item in the top box down
	void DelClicked     ( wxCommandEvent& event ) { m_libLst->RemoveSelected(); }//remove all selected items in the bottom box
	void FltClicked     ( wxCommandEvent& event );//update filters
	void SearchChanged  ( wxCommandEvent& event ) { m_libLst->setSearch(event.GetString()); }
	void IdxUnchecked   ( wxListEvent   & event );//move top box items to bottom box on untick
	void KnownChecked   ( wxListEvent   & event );//move bottom box items to top box on tick
	void BrowseClicked  ( wxCommandEvent& event );//search for new master pattern files
	void DirClicked     ( wxCommandEvent& event );//add a folder of files to the bottom box

	protected:

		void ClearSearch() {m_mpSearch->SetValue("");}

		//@brief    : add a pattern file to the control
		//@param str: pattern file to add
		//@param top: add to top box (use) instead of bottom (library)
		//@return   : status string
		std::string AddPatternFile(wxString str, bool top = false);

	public:

		//@brief : check if the panel is currently valid
		//@return: "" if valid, reason invalid otherwise
		std::string validMsg() const;

		//@brief : get a list of all known master pattern files
		//@return: full paths to known pattern files
		std::vector<wxString> getLibrary() const;

		//@brief    : set list of known master pattern files
		//@param lib: full paths to known pattern files
		void setLibrary(const std::vector<wxString>& lib) {for(const wxString& str : lib) AddPatternFile(str);}

		std::vector<std::string> getSelected() const;

		MasterPatternSelectPanel( wxWindow* parent, wxWindowID id = wxID_ANY, const wxPoint& pos = wxDefaultPosition, const wxSize& size = wxSize( 300,300 ), long style = wxTAB_TRAVERSAL, const wxString& name = wxEmptyString );
		~MasterPatternSelectPanel();
};

///////////////////////////////////////////////////////////////////////////
// C++ code generated with wxFormBuilder (version Oct 26 2018)
// http://www.wxformbuilder.org/
//
// PLEASE DO *NOT* EDIT THIS FILE!
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////

#include "SinglePanelDialog.h"
#include "MasterFileFilterPanel.h"

void MasterPatternSelectPanel::FltClicked( wxCommandEvent& event ) {
	//start by getting existings filters
	std::pair<float, float> kv, tlt;
	std::pair<int, int> sg;
	ElementMask el;
	bool ext;
	m_libLst->getFilters(kv, tlt, el, ext, sg);
	if(MasterFileList::AllKv () == kv ) {kv .first = kv .second = NAN;}
	if(MasterFileList::AllTlt() == tlt) {tlt.first = tlt.second = NAN;}
	SinglePanelDialog<MasterFileFilterPanel> dlg(this, wxID_ANY, "Select Filters");
	dlg.getPanel()->setBounds(kv, tlt, el, ext, sg);
	if(wxID_OK == dlg.ShowModal()) {
		dlg.getPanel()->getBounds(kv, tlt, el, ext, sg);
		if(kv .first != kv .first) kv  = MasterFileList::AllKv ();
		if(tlt.first != tlt.first) tlt = MasterFileList::AllTlt();
		m_libLst->setFilters(kv, tlt, el, ext, sg);
	}
} 

//move top box items to bottom box on untick
void MasterPatternSelectPanel::IdxUnchecked( wxListEvent& event ) {
	m_libLst->AddItem(m_selLst->GetItem(event.GetIndex()), false);
	m_selLst->RemoveItem(event.GetIndex());
	testValid();
}

//move bottom box items to top box on tick
void MasterPatternSelectPanel::KnownChecked( wxListEvent& event ) {
	if(wxEVT_LIST_ITEM_ACTIVATED == event.GetEventType()) {
		std::vector<long> sel = m_libLst->GetSelection();
		for(const long& i : sel) m_selLst->AddItem(m_libLst->GetItem(i), true );
		m_libLst->RemoveSelected();
	} else {//check box event
		m_selLst->AddItem(m_libLst->GetItem(event.GetIndex()), true );
		m_libLst->RemoveItem(event.GetIndex());
	}
	testValid();
}

//@brief: search for new master pattern files
void MasterPatternSelectPanel::BrowseClicked  ( wxCommandEvent& event ) {
	//browse for file
	wxFileDialog openFileDlg(this, _("Add Master Pattern File"), wxEmptyString, wxEmptyString, "SHT files (*.sht)|*.sht", wxFD_OPEN|wxFD_FILE_MUST_EXIST);
	if(openFileDlg.ShowModal() == wxID_CANCEL) return;

	//do the file loading
	std::string msg = AddPatternFile(openFileDlg.GetPath(), true);
	if(!msg.empty()) {
		wxMessageDialog msgDlg(this, msg, "Error Adding File");
		msgDlg.ShowModal();
	}
}

#include <wx/dir.h>
#include <wx/dirdlg.h>

//@brief: add a folder of files to the bottom box
void MasterPatternSelectPanel::DirClicked     ( wxCommandEvent& event ) {
	//browse for folder
	wxDirDialog openDirDlg(this, _("Add Master Pattern Folder"), wxEmptyString, wxDD_DIR_MUST_EXIST);
	if(openDirDlg.ShowModal() == wxID_CANCEL) return;

	//loop over all files in folder
	wxArrayString files;
	wxDir::GetAllFiles(openDirDlg.GetPath(), &files, "*.sht", wxDIR_DIRS | wxDIR_FILES);
	std::vector<wxString> added;
	for(const wxString& str : files) {
		if(AddPatternFile(str, false).empty()) added.push_back(str);
	}
}

//@brief : check if the panel is currently valid
//@return: "" if valid, reason invalid otherwise
std::string MasterPatternSelectPanel::validMsg() const {
	size_t numIt = m_selLst->GetItemCount();
	if(0 == numIt) return "at least one master pattern required (double click/enter to move)";
	else if (1 == numIt) return "";
	MasterFile mf0 = m_selLst->GetItem(0);
	for(size_t i = 1; i < numIt; i++) {
		MasterFile mfi = m_selLst->GetItem(i);
		if(mfi.kv  != mf0.kv ) return "all master patterns must have the same kv";
		if(mfi.tlt != mf0.tlt) return "all master patterns must have the same tilt";
	}
	return "";
}

//@brief : get a list of all known master pattern files
//@return: full paths to known pattern files
std::vector<wxString> MasterPatternSelectPanel::getLibrary() const {
	std::vector<wxString> sel = m_selLst->getAllFiles();
	std::vector<wxString> knw = m_libLst->getAllFiles();
	knw.insert(knw.end(), sel.begin(), sel.end());
	return knw;
}

std::vector<std::string> MasterPatternSelectPanel::getSelected() const {
	std::vector<wxString> sel = m_selLst->getAllFiles();
	std::vector<std::string> ret;
	for(const wxString& str : sel) ret.push_back(str.ToStdString());
	return ret;
}

MasterPatternSelectPanel::MasterPatternSelectPanel( wxWindow* parent, wxWindowID id, const wxPoint& pos, const wxSize& size, long style, const wxString& name ) : ValidityPanel( parent, id, pos, size, style, name ) {
	//split panel into 2 vertical boxes
	wxBoxSizer      * bMP      = new wxBoxSizer      ( wxVERTICAL );
	wxStaticBoxSizer* sbMP     = new wxStaticBoxSizer( new wxStaticBox( this, wxID_ANY, wxT("Indexing Master Patterns") ), wxVERTICAL );
	wxStaticBoxSizer* sbLibMP  = new wxStaticBoxSizer( new wxStaticBox( this, wxID_ANY, wxT("Master Pattern Library"  ) ), wxVERTICAL );
	bMP->Add( sbMP   , 1, wxEXPAND, 5 );
	bMP->Add( sbLibMP, 1, wxEXPAND, 5 );

	//build sizer for strip of bottons to each panel
	wxBoxSizer     * hSizerMP  = new wxBoxSizer( wxHORIZONTAL );
	wxBoxSizer     * hSizerKwn = new wxBoxSizer( wxHORIZONTAL );

	//elements for top box
	m_selLst   = new MasterFileList( sbMP   ->GetStaticBox()                                                                                                                         );
	m_btnUp    = new wxBitmapButton( sbMP   ->GetStaticBox(), wxID_ANY, wxArtProvider::GetBitmap( wxART_GO_UP      , wxART_BUTTON ), wxDefaultPosition, wxDefaultSize, wxBU_AUTODRAW );
	m_btnDn    = new wxBitmapButton( sbMP   ->GetStaticBox(), wxID_ANY, wxArtProvider::GetBitmap( wxART_GO_DOWN    , wxART_BUTTON ), wxDefaultPosition, wxDefaultSize, wxBU_AUTODRAW );
	m_btnAdd   = new wxBitmapButton( sbMP   ->GetStaticBox(), wxID_ANY, wxArtProvider::GetBitmap( wxART_FILE_OPEN  , wxART_BUTTON ), wxDefaultPosition, wxDefaultSize, wxBU_AUTODRAW );
	m_mpSearch = new wxSearchCtrl  ( sbLibMP->GetStaticBox(), wxID_ANY, wxEmptyString                                              , wxDefaultPosition, wxDefaultSize, 0             );
	m_libLst   = new MasterFileList( sbLibMP->GetStaticBox()                                                                                                                         );
	m_btnDel   = new wxBitmapButton( sbLibMP->GetStaticBox(), wxID_ANY, wxArtProvider::GetBitmap( wxART_DELETE     , wxART_BUTTON ), wxDefaultPosition, wxDefaultSize, wxBU_AUTODRAW );
	m_btnFlt   = new wxBitmapButton( sbLibMP->GetStaticBox(), wxID_ANY, wxArtProvider::GetBitmap( wxART_FIND       , wxART_BUTTON ), wxDefaultPosition, wxDefaultSize, wxBU_AUTODRAW );
	m_btnDir   = new wxBitmapButton( sbLibMP->GetStaticBox(), wxID_ANY, wxArtProvider::GetBitmap( wxART_FOLDER_OPEN, wxART_BUTTON ), wxDefaultPosition, wxDefaultSize, wxBU_AUTODRAW );  

	m_selLst->setSort(MasterFile::CompFunc());
	#ifndef __WXMAC__
	m_mpSearch->ShowSearchButton( true );
	#endif
	m_mpSearch->ShowCancelButton( true );
	m_libLst->EnableMultipleSelection(true);

	//finishin assmebling sizers
	sbMP   ->Add( m_selLst  , 1, wxALL|wxEXPAND, 5 );
	sbMP   ->Add( hSizerMP  , 0,       wxEXPAND, 5 );
	sbLibMP->Add( m_mpSearch, 0, wxALL|wxEXPAND, 5 );
	sbLibMP->Add( m_libLst  , 1, wxALL|wxEXPAND, 5 );
	sbLibMP->Add( hSizerKwn , 0,       wxEXPAND, 5 );

	//assemble indexing button strip
	hSizerMP ->Add( m_btnUp , 0, wxALIGN_CENTER_VERTICAL|wxLEFT , 5 );
	hSizerMP ->Add( m_btnDn , 0, wxALIGN_CENTER_VERTICAL|wxLEFT , 5 );
	hSizerMP ->Add( 0       , 0, 1, wxEXPAND                    , 5 );
	hSizerMP ->Add( m_btnAdd, 0, wxALIGN_CENTER_VERTICAL|wxRIGHT, 5 );

	//assemble library button strip
	hSizerKwn->Add( m_btnDel, 0, wxALIGN_CENTER_VERTICAL|wxRIGHT, 5 );
	hSizerKwn->Add( m_btnFlt, 0, wxALIGN_CENTER_VERTICAL|wxRIGHT, 5 );
	hSizerKwn->Add( 0       , 0, 1, wxEXPAND                    , 5 );
	hSizerKwn->Add( m_btnDir, 0, wxALIGN_CENTER_VERTICAL|wxRIGHT, 5 );

	this->SetSizer( bMP );
	this->Layout();

	// Connect Events
	m_selLst  ->Connect( wxEVT_LIST_ITEM_UNCHECKED   , wxListEventHandler   ( MasterPatternSelectPanel::IdxUnchecked  ), NULL, this );
	m_btnUp   ->Connect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( MasterPatternSelectPanel::UpClicked     ), NULL, this );
	m_btnDn   ->Connect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( MasterPatternSelectPanel::DownClicked   ), NULL, this );
	m_btnAdd  ->Connect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( MasterPatternSelectPanel::BrowseClicked ), NULL, this );
	m_btnDel  ->Connect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( MasterPatternSelectPanel::DelClicked    ), NULL, this );
	m_btnFlt  ->Connect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( MasterPatternSelectPanel::FltClicked    ), NULL, this );
	m_btnDir  ->Connect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( MasterPatternSelectPanel::DirClicked    ), NULL, this );
	m_mpSearch->Connect( wxEVT_COMMAND_TEXT_UPDATED  , wxCommandEventHandler( MasterPatternSelectPanel::SearchChanged ), NULL, this );
	m_libLst  ->Connect( wxEVT_LIST_ITEM_CHECKED     , wxListEventHandler   ( MasterPatternSelectPanel::KnownChecked  ), NULL, this );
	m_selLst  ->Connect( wxEVT_LIST_ITEM_ACTIVATED   , wxListEventHandler   ( MasterPatternSelectPanel::IdxUnchecked  ), NULL, this );
	m_libLst  ->Connect( wxEVT_LIST_ITEM_ACTIVATED   , wxListEventHandler   ( MasterPatternSelectPanel::KnownChecked  ), NULL, this );
}

MasterPatternSelectPanel::~MasterPatternSelectPanel() {
	// Disconnect Events
	m_selLst  ->Disconnect( wxEVT_LIST_ITEM_UNCHECKED   , wxListEventHandler   ( MasterPatternSelectPanel::IdxUnchecked  ), NULL, this );
	m_btnUp   ->Disconnect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( MasterPatternSelectPanel::UpClicked     ), NULL, this );
	m_btnDn   ->Disconnect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( MasterPatternSelectPanel::DownClicked   ), NULL, this );
	m_btnAdd  ->Disconnect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( MasterPatternSelectPanel::BrowseClicked ), NULL, this );
	m_btnDel  ->Disconnect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( MasterPatternSelectPanel::DelClicked    ), NULL, this );
	m_btnFlt  ->Disconnect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( MasterPatternSelectPanel::FltClicked    ), NULL, this );
	m_btnDir  ->Disconnect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( MasterPatternSelectPanel::DirClicked    ), NULL, this );
	m_mpSearch->Disconnect( wxEVT_COMMAND_TEXT_UPDATED  , wxCommandEventHandler( MasterPatternSelectPanel::SearchChanged ), NULL, this );
	m_libLst  ->Disconnect( wxEVT_LIST_ITEM_CHECKED     , wxListEventHandler   ( MasterPatternSelectPanel::KnownChecked  ), NULL, this );
	m_selLst  ->Disconnect( wxEVT_LIST_ITEM_ACTIVATED   , wxListEventHandler   ( MasterPatternSelectPanel::IdxUnchecked  ), NULL, this );
	m_libLst  ->Disconnect( wxEVT_LIST_ITEM_ACTIVATED   , wxListEventHandler   ( MasterPatternSelectPanel::KnownChecked  ), NULL, this );
}

//@brief    : add a pattern file to the control
//@param str: pattern file to add
//@param top: add to top box (use) instead of bottom (library)
//@return   : status string
std::string MasterPatternSelectPanel::AddPatternFile(wxString str, bool top) {
	if(m_selLst->HasItem(str)) return "file already in indexing list";
	if(top) {//
		if(m_libLst->HasItem(str)) return "file already in library (please add from library instead)";
		bool added = m_selLst->AddItem(MasterFile(str), true );
		testValid();
		return added ? "" : "couldn't read master pattern from file";
	} else {
		if(m_libLst->HasItem(str)) return "file already in library";
		// ClearSearch();//it is probably confusing to users if the new file doesn't show up
		return m_libLst->AddItem(MasterFile(str), false) ? "" : "couldn't read master pattern from file";
	}
}

#endif//_MP_SEL_H_
