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

#ifndef _BIBTEX_DLG_H_
#define _BIBTEX_DLG_H_

#include <wx/artprov.h>
#include <wx/string.h>
#include <wx/stattext.h>
#include <wx/gdicmn.h>
#include <wx/textctrl.h>
#include <wx/bitmap.h>
#include <wx/image.h>
#include <wx/icon.h>
#include <wx/button.h>
#include <wx/checkbox.h>
#include <wx/sizer.h>
#include <wx/dialog.h>
#include <wx/clipbrd.h>


///////////////////////////////////////////////////////////////////////////////
/// Class BibtexDialog
///////////////////////////////////////////////////////////////////////////////
class BibtexDialog : public wxDialog {
	private:

	protected:
		wxTextCtrl* m_txt   ;
		wxButton  * m_btnCpy;
		wxCheckBox* m_chk   ;
		wxButton  * m_ok    ;

		void CopyData     ( wxCommandEvent& event ) { 
			if (wxTheClipboard->Open()) {
				wxTheClipboard->SetData( new wxTextDataObject(m_txt->GetValue()) );
				wxTheClipboard->Close();
			}
		}
		void DismissDialog( wxCommandEvent& event ) { Close(); }

	public:

		bool Silence() const {return m_chk->GetValue();}

		BibtexDialog( wxWindow* parent, wxWindowID id = wxID_ANY, const wxString& title = wxEmptyString, const wxPoint& pos = wxDefaultPosition, const wxSize& size = wxSize( 475,325 ), long style = wxDEFAULT_DIALOG_STYLE );
		~BibtexDialog();

};

///////////////////////////////////////////////////////////////////////////

BibtexDialog::BibtexDialog( wxWindow* parent, wxWindowID id, const wxString& title, const wxPoint& pos, const wxSize& size, long style ) : wxDialog( parent, id, title, pos, size, style ) {
	this->SetSizeHints( wxDefaultSize, wxDefaultSize );

	wxBoxSizer* vSizer = new wxBoxSizer( wxVERTICAL   );
	wxBoxSizer* hSizer = new wxBoxSizer( wxHORIZONTAL );

	wxStaticText* staTxt = new wxStaticText( this, wxID_ANY, wxT("If you find this software useful, please consider citing the corresponding paper:"), wxDefaultPosition, wxDefaultSize, 0 ); staTxt->Wrap( -1 );

	std::string key = "@article{EMSphInx,\n\ttitle = \"A spherical harmonic transform approach to the indexing of electron back-scattered diffraction patterns\",\n\tjournal = \"Ultramicroscopy\",\n\tvolume = \"207\",\n\tpages = \"112841\",\n\tyear = \"2019\",\n\tissn = \"0304-3991\",\n\tdoi = \"https://doi.org/10.1016/j.ultramic.2019.112841\",\n\tauthor = \"W.C. Lenthe and S. Singh and M. De Graef\",\n}";
	m_txt    = new wxTextCtrl( this, wxID_ANY, key                                 , wxDefaultPosition, wxDefaultSize, wxTE_DONTWRAP|wxTE_MULTILINE|wxTE_READONLY );
	m_btnCpy = new wxButton  ( this, wxID_ANY, wxT("Copy to Clipboard"            ), wxDefaultPosition, wxDefaultSize, 0                                          );
	m_chk    = new wxCheckBox( this, wxID_ANY, wxT("Don't show this message again"), wxDefaultPosition, wxDefaultSize, 0                                          );
	m_ok     = new wxButton  ( this, wxID_ANY, wxT("Dismiss"                      ), wxDefaultPosition, wxDefaultSize, 0                                          );
	m_btnCpy->SetBitmap( wxArtProvider::GetBitmap( wxART_COPY, wxART_BUTTON ) );
	m_ok->SetDefault();
	
	vSizer->Add( staTxt  , 0, wxALL                          , 5 );
	vSizer->Add( m_txt   , 1, wxALL|wxEXPAND                 , 5 );
	vSizer->Add( m_btnCpy, 0, wxALL|wxALIGN_CENTER_HORIZONTAL, 5 );

	hSizer->AddStretchSpacer();
	hSizer->Add( m_chk   , 0, wxALL|wxALIGN_CENTER_VERTICAL  , 5 );
	hSizer->AddStretchSpacer();
	hSizer->Add( m_ok    , 0, wxALL|wxALIGN_CENTER_VERTICAL  , 5 );
	hSizer->AddStretchSpacer();

	vSizer->Add( hSizer, 0, wxEXPAND, 5 );

	this->SetSizer( vSizer );
	this->Layout();

	this->Centre( wxBOTH );

	// Connect Events
	m_btnCpy->Connect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( BibtexDialog::CopyData      ), NULL, this );
	m_ok    ->Connect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( BibtexDialog::DismissDialog ), NULL, this );
}

BibtexDialog::~BibtexDialog() {
	m_btnCpy->Disconnect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( BibtexDialog::CopyData      ), NULL, this );
	m_ok    ->Disconnect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( BibtexDialog::DismissDialog ), NULL, this );
}


#endif//_BIBTEX_DLG_H_
