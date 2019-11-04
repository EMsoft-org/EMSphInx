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

#ifndef _MP_CVT_DLG_H_
#define _MP_CVT_DLG_H_

#include <wx/string.h>
#include <wx/stattext.h>
#include <wx/filepicker.h>
#include <wx/textctrl.h>
#include <wx/button.h>
#include <wx/sizer.h>
#include <wx/dialog.h>

///////////////////////////////////////////////////////////////////////////////
/// Class MpConvertDialog
///////////////////////////////////////////////////////////////////////////////
class MpConvertDialog : public wxDialog {
	private:

	protected:
		wxFilePickerCtrl* m_fileIn ;
		wxTextCtrl      * m_txtForm;
		wxTextCtrl      * m_txtName;
		wxTextCtrl      * m_txtSSyb;
		wxFilePickerCtrl* m_fileOut;
		wxButton        * m_button ;

		// Virtual event handlers, overide them in your derived class
		virtual void OnInputChanged( wxFileDirPickerEvent& event ) { event.Skip(); }
		virtual void OnFormulaChanged( wxCommandEvent& event ) { event.Skip(); }
		virtual void OnOutputChanged( wxFileDirPickerEvent& event ) { event.Skip(); }
		virtual void OnConvert( wxCommandEvent& event ) { event.Skip(); }


	public:

		MpConvertDialog( wxWindow* parent, wxWindowID id = wxID_ANY, const wxString& title = wxT("Convert EMsoft Master Pattern"), const wxPoint& pos = wxDefaultPosition, const wxSize& size = wxSize( 300,350 ), long style = wxDEFAULT_DIALOG_STYLE );
		~MpConvertDialog();

};

MpConvertDialog::MpConvertDialog( wxWindow* parent, wxWindowID id, const wxString& title, const wxPoint& pos, const wxSize& size, long style ) : wxDialog( parent, id, title, pos, size, style ) {
	this->SetSizeHints( wxDefaultSize, wxDefaultSize );

	wxFlexGridSizer* fgSizerfgSizer = new wxFlexGridSizer( 6, 2, 0, 0 );
	fgSizer->AddGrowableCol( 1 );
	fgSizer->AddGrowableRow( 0 ); fgSizer->AddGrowableRow( 1 );
	fgSizer->AddGrowableRow( 2 ); fgSizer->AddGrowableRow( 3 );
	fgSizer->AddGrowableRow( 4 ); fgSizer->AddGrowableRow( 5 );
	fgSizer->SetFlexibleDirection( wxHORIZONTAL );
	fgSizer->SetNonFlexibleGrowMode( wxFLEX_GROWMODE_ALL );

	wxStaticText* staTxtIn   = new wxStaticText( this, wxID_ANY, wxT("Input File" ), wxDefaultPosition, wxDefaultSize, 0 ); staTxtIn  ->Wrap( -1 );
	wxStaticText* staTxtForm = new wxStaticText( this, wxID_ANY, wxT("Formula"    ), wxDefaultPosition, wxDefaultSize, 0 ); staTxtForm->Wrap( -1 );
	wxStaticText* staTxtName = new wxStaticText( this, wxID_ANY, wxT("Name"       ), wxDefaultPosition, wxDefaultSize, 0 ); staTxtName->Wrap( -1 );
	wxStaticText* staTxtSSyb = new wxStaticText( this, wxID_ANY, wxT("S. Syb."    ), wxDefaultPosition, wxDefaultSize, 0 ); staTxtSSyb->Wrap( -1 );
	wxStaticText* staTxtOut  = new wxStaticText( this, wxID_ANY, wxT("Output File"), wxDefaultPosition, wxDefaultSize, 0 ); staTxtOut ->Wrap( -1 );

	m_fileIn  = new wxFilePickerCtrl( this, wxID_ANY, wxEmptyString, wxT("Input Master Pattern" ), wxT("*.h5"), wxDefaultPosition, wxDefaultSize, wxFLP_DEFAULT_STYLE|wxFLP_SMALL );
	m_txtForm = new wxTextCtrl      ( this, wxID_ANY, wxEmptyString, wxDefaultPosition, wxDefaultSize, 0 );
	m_txtName = new wxTextCtrl      ( this, wxID_ANY, wxEmptyString, wxDefaultPosition, wxDefaultSize, 0 );
	m_txtSSyb = new wxTextCtrl      ( this, wxID_ANY, wxEmptyString, wxDefaultPosition, wxDefaultSize, 0 );
	m_fileOut = new wxFilePickerCtrl( this, wxID_ANY, wxEmptyString, wxT("Output Master Pattern"), wxT("*.sht"), wxDefaultPosition, wxDefaultSize, wxFLP_OVERWRITE_PROMPT|wxFLP_SAVE|wxFLP_SMALL|wxFLP_USE_TEXTCTRL );

	m_button = new wxButton( this, wxID_ANY, wxT("Convert"), wxDefaultPosition, wxDefaultSize, 0 );

	m_button->SetDefault();

	fgSizer->Add( staTxtIn  , 0, wxALL|wxALIGN_CENTER_VERTICAL              , 5 );
	fgSizer->Add( m_fileIn  , 0, wxALL|wxALIGN_CENTER_VERTICAL|wxEXPAND     , 5 );
	fgSizer->Add( staTxtForm, 0, wxALL|wxALIGN_CENTER_VERTICAL              , 5 );
	fgSizer->Add( m_txtForm , 0, wxALL|wxALIGN_CENTER_VERTICAL|wxEXPAND     , 5 );
	fgSizer->Add( staTxtName, 0, wxALL|wxALIGN_CENTER_VERTICAL              , 5 );
	fgSizer->Add( m_txtName , 0, wxALL|wxALIGN_CENTER_VERTICAL|wxEXPAND     , 5 );
	fgSizer->Add( staTxtSSyb, 0, wxALL|wxALIGN_CENTER_VERTICAL              , 5 );
	fgSizer->Add( m_txtSSyb , 0, wxALL|wxALIGN_CENTER_VERTICAL|wxEXPAND     , 5 );
	fgSizer->Add( staTxtOut , 0, wxALL|wxALIGN_CENTER_VERTICAL              , 5 );
	fgSizer->Add( m_fileOut , 0, wxALL|wxALIGN_CENTER_VERTICAL|wxEXPAND     , 5 );
	fgSizer->Add( 0, 0, 1, wxEXPAND);
	fgSizer->Add( m_button  , 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT, 5 );

	this->SetSizer( fgSizer );
	this->Layout();
	this->Centre( wxBOTH );

	// Connect Events
	m_fileIn ->Connect( wxEVT_COMMAND_FILEPICKER_CHANGED, wxFileDirPickerEventHandler( MpConvertDialog::OnInputChanged   ), NULL, this );
	m_txtForm->Connect( wxEVT_COMMAND_TEXT_UPDATED      , wxCommandEventHandler      ( MpConvertDialog::OnFormulaChanged ), NULL, this );
	m_fileOut->Connect( wxEVT_COMMAND_FILEPICKER_CHANGED, wxFileDirPickerEventHandler( MpConvertDialog::OnOutputChanged  ), NULL, this );
	m_button ->Connect( wxEVT_COMMAND_BUTTON_CLICKED    , wxCommandEventHandler      ( MpConvertDialog::OnConvert        ), NULL, this );
}

MpConvertDialog::~MpConvertDialog() {
	// Disconnect Events
	m_fileIn ->Disconnect( wxEVT_COMMAND_FILEPICKER_CHANGED, wxFileDirPickerEventHandler( MpConvertDialog::OnInputChanged   ), NULL, this );
	m_txtForm->Disconnect( wxEVT_COMMAND_TEXT_UPDATED      , wxCommandEventHandler      ( MpConvertDialog::OnFormulaChanged ), NULL, this );
	m_fileOut->Disconnect( wxEVT_COMMAND_FILEPICKER_CHANGED, wxFileDirPickerEventHandler( MpConvertDialog::OnOutputChanged  ), NULL, this );
	m_button ->Disconnect( wxEVT_COMMAND_BUTTON_CLICKED    , wxCommandEventHandler      ( MpConvertDialog::OnConvert        ), NULL, this );
}

#endif//_MP_CVT_DLG_H_
