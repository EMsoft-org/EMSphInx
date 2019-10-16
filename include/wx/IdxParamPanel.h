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

#ifndef _IDX_PARAM_PAN_H_
#define _IDX_PARAM_PAN_H_

#include <wx/artprov.h>
#include <wx/xrc/xmlres.h>
#include <wx/string.h>
#include <wx/stattext.h>
#include <wx/gdicmn.h>
#include <wx/font.h>
#include <wx/colour.h>
#include <wx/settings.h>
#include <wx/textctrl.h>
#include <wx/bitmap.h>
#include <wx/image.h>
#include <wx/icon.h>
#include <wx/button.h>
#include <wx/checkbox.h>
#include <wx/sizer.h>
#include <wx/statbox.h>
#include <wx/filepicker.h>
#include <wx/panel.h>
#include <wx/valnum.h>

#include "wx/ValidityPanel.h"

///////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
/// Class IdxParamPanel
///////////////////////////////////////////////////////////////////////////////
class IdxParamPanel : public ValidityPanel
{
	private:
		wxIntegerValidator      <int  > valBw;//bandwidth

	protected:
		wxTextCtrl      * m_txtBw   ;
		wxButton        * m_btnBwPrv;
		wxCheckBox      * m_chkNrm  ;
		wxCheckBox      * m_chkRef  ;
		wxFilePickerCtrl* m_fpData  ;
		wxFilePickerCtrl* m_fpVen   ;
		wxFilePickerCtrl* m_fpIpf   ;
		wxFilePickerCtrl* m_fpCi    ;

		// Virtual event handlers, overide them in your derived class
		virtual void BandwidthChanged( wxCommandEvent      & event ) { testValid();  }
		virtual void OnPreview       ( wxCommandEvent      & event ) { event.Skip(); }
		virtual void DataFileChanged ( wxFileDirPickerEvent& event ) { testValid();  }

	public:

		IdxParamPanel( wxWindow* parent, wxWindowID id = wxID_ANY, const wxPoint& pos = wxDefaultPosition, const wxSize& size = wxSize( 400,300 ), long style = wxTAB_TRAVERSAL, const wxString& name = wxEmptyString );
		~IdxParamPanel();

		//@brief : get bandwidth
		//@return: bandwidth
		int getBw() const {long v; m_txtBw->GetLineText(0).ToLong(&v); return v;}

		//@brief : get normalized flag
		//@return: normalized flag
		bool getNorm() const {return m_chkNrm->GetValue();}

		//@brief : get refinement flag
		//@return: refinement flag
		bool getRef() const {return m_chkRef->GetValue();}

		//@brief : get the data file
		//@return: data file
		wxString getDataFile() const {return m_fpData->GetPath();}

		//@brief : get the vendor file
		//@return: vendor file
		wxString getVendorFile() const {return m_fpVen->GetPath();}

		//@brief : get the ipf file
		//@return: ipf file
		wxString getIpfFile() const {return m_fpIpf->GetPath();}

		//@brief : get the CI file
		//@return: CI file
		wxString getCiFile() const {return m_fpCi->GetPath();}

		//require bandwidth + output file

		//@brief: sanity check the current state
		//@return: true if the values parsed from the panel are reasonable, false otherwise
		//@note  : checks for has a file, has detector sizes, and has an AHE value
		std::string validMsg() const {
			if(m_txtBw->GetLineText(0).IsEmpty()) return "bandwidth empty";
			if(getDataFile().IsEmpty()) return "data file empty";
			return "";
		}

};

///////////////////////////////////////////////////////////////////////////
// C++ code generated with wxFormBuilder (version Oct 26 2018)
// http://www.wxformbuilder.org/
//
// PLEASE DO *NOT* EDIT THIS FILE!
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////

IdxParamPanel::IdxParamPanel( wxWindow* parent, wxWindowID id, const wxPoint& pos, const wxSize& size, long style, const wxString& name ) : ValidityPanel( parent, id, pos, size, style, name ) {
	//split panel into 2 vertical boxes
	wxBoxSizer      * bIdxMisc  = new wxBoxSizer( wxVERTICAL );
	wxStaticBoxSizer* sbIdxPrm  = new wxStaticBoxSizer( new wxStaticBox( this, wxID_ANY, wxT("Indexing Parameters") ), wxVERTICAL );
	wxStaticBoxSizer* sbOutFile = new wxStaticBoxSizer( new wxStaticBox( this, wxID_ANY, wxT("Output Files"       ) ), wxVERTICAL );
	bIdxMisc->Add( sbIdxPrm , 1, wxEXPAND, 5 );
	bIdxMisc->Add( sbOutFile, 2, wxEXPAND, 5 );

	//top box is 6x2 grid
	wxFlexGridSizer* fgIdxPrm = new wxFlexGridSizer( 2, 6, 0, 0 );
	fgIdxPrm->AddGrowableCol( 0 ); fgIdxPrm->AddGrowableCol( 2 ); fgIdxPrm->AddGrowableCol( 5 );
	fgIdxPrm->AddGrowableRow( 0 ); fgIdxPrm->AddGrowableRow( 1 );
	fgIdxPrm->SetFlexibleDirection( wxHORIZONTAL );
	fgIdxPrm->SetNonFlexibleGrowMode( wxFLEX_GROWMODE_SPECIFIED );
	sbIdxPrm->Add( fgIdxPrm, 1, wxEXPAND, 5 );

	//bottom box is 4x2 grid
	wxFlexGridSizer* fgOutFile = new wxFlexGridSizer( 4, 2, 0, 0 );
	fgOutFile->AddGrowableCol( 1 );
	fgOutFile->AddGrowableRow( 0 ); fgOutFile->AddGrowableRow( 1 ); fgOutFile->AddGrowableRow( 2 ); fgOutFile->AddGrowableRow( 3 );
	fgOutFile->SetFlexibleDirection( wxHORIZONTAL );
	fgOutFile->SetNonFlexibleGrowMode( wxFLEX_GROWMODE_SPECIFIED );
	sbOutFile->Add( fgOutFile, 1, wxEXPAND, 5 );

	//build all labels
	wxStaticText* staTxtBw       = new wxStaticText( sbIdxPrm ->GetStaticBox(), wxID_ANY, wxT("Bandwidth"  ), wxDefaultPosition, wxDefaultSize, 0 ); staTxtBw      ->Wrap(-1);
	wxStaticText* staTxtDataFile = new wxStaticText( sbOutFile->GetStaticBox(), wxID_ANY, wxT("Data File"  ), wxDefaultPosition, wxDefaultSize, 0 ); staTxtDataFile->Wrap(-1);
	wxStaticText* staTxtVenFile  = new wxStaticText( sbOutFile->GetStaticBox(), wxID_ANY, wxT("Vendor File"), wxDefaultPosition, wxDefaultSize, 0 ); staTxtVenFile ->Wrap(-1);
	wxStaticText* staTxtIpf      = new wxStaticText( sbOutFile->GetStaticBox(), wxID_ANY, wxT("IPF Map"    ), wxDefaultPosition, wxDefaultSize, 0 ); staTxtIpf     ->Wrap(-1);
	wxStaticText* staTxtCi       = new wxStaticText( sbOutFile->GetStaticBox(), wxID_ANY, wxT("CI Map"     ), wxDefaultPosition, wxDefaultSize, 0 ); staTxtCi      ->Wrap(-1);

	//set text validator ranges
	valBw.SetRange(1, 350);

	//build elements for indexing parameters box
	m_txtBw    = new wxTextCtrl      ( sbIdxPrm ->GetStaticBox(), wxID_ANY, wxT(""          ), wxDefaultPosition, wxDefaultSize, wxTE_CENTER, valBw );
	m_btnBwPrv = new wxButton        ( sbIdxPrm ->GetStaticBox(), wxID_ANY, wxT("Preview..."), wxDefaultPosition, wxDefaultSize, 0                  );
	m_chkNrm   = new wxCheckBox      ( sbIdxPrm ->GetStaticBox(), wxID_ANY, wxT("Normalized"), wxDefaultPosition, wxDefaultSize, 0                  );
	m_chkRef   = new wxCheckBox      ( sbIdxPrm ->GetStaticBox(), wxID_ANY, wxT("Refinement"), wxDefaultPosition, wxDefaultSize, 0                  );
	m_btnBwPrv->Enable(false);

	//build elements for output file box
	m_fpData   = new wxFilePickerCtrl( sbOutFile->GetStaticBox(), wxID_ANY, wxEmptyString, wxT("Select a file"), wxT("*.h5"       ), wxDefaultPosition, wxDefaultSize, wxFLP_OVERWRITE_PROMPT|wxFLP_SAVE|wxFLP_SMALL|wxFLP_USE_TEXTCTRL );
	m_fpVen    = new wxFilePickerCtrl( sbOutFile->GetStaticBox(), wxID_ANY, wxEmptyString, wxT("Select a file"), wxT("*.ang;*.ctf"), wxDefaultPosition, wxDefaultSize, wxFLP_OVERWRITE_PROMPT|wxFLP_SAVE|wxFLP_SMALL|wxFLP_USE_TEXTCTRL );
	m_fpIpf    = new wxFilePickerCtrl( sbOutFile->GetStaticBox(), wxID_ANY, wxEmptyString, wxT("Select a file"), wxT("*.png"      ), wxDefaultPosition, wxDefaultSize,                        wxFLP_SAVE|wxFLP_SMALL|wxFLP_USE_TEXTCTRL );
	m_fpCi     = new wxFilePickerCtrl( sbOutFile->GetStaticBox(), wxID_ANY, wxEmptyString, wxT("Select a file"), wxT("*.png"      ), wxDefaultPosition, wxDefaultSize,                        wxFLP_SAVE|wxFLP_SMALL|wxFLP_USE_TEXTCTRL );

	//assemble indexing parameter grid
	fgIdxPrm ->Add( 0             , 0, 1    , wxEXPAND                                        , 5 );
	fgIdxPrm ->Add( staTxtBw      , 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_CENTER_HORIZONTAL, 5 );
	fgIdxPrm ->Add( 0             , 0, 1    , wxEXPAND                                        , 5 );
	fgIdxPrm ->Add( m_txtBw       , 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_CENTER_HORIZONTAL, 5 );
	fgIdxPrm ->Add( m_btnBwPrv    , 0, wxALL|wxALIGN_CENTER_VERTICAL                          , 5 );
	fgIdxPrm ->Add( 0             , 0, 1    , wxEXPAND                                        , 5 );
	fgIdxPrm ->Add( 0             , 0, 1    , wxEXPAND                                        , 5 );
	fgIdxPrm ->Add( m_chkNrm      , 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_CENTER_HORIZONTAL, 5 );
	fgIdxPrm ->Add( 0             , 0, 1    , wxEXPAND                                        , 5 );
	fgIdxPrm ->Add( m_chkRef      , 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_CENTER_HORIZONTAL, 5 );
	fgIdxPrm ->Add( 0             , 0, 1    , wxEXPAND                                        , 5 );
	fgIdxPrm ->Add( 0             , 0, 1    , wxEXPAND                                        , 5 );

	//assemble output file grid
	fgOutFile->Add( staTxtDataFile, 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_CENTER_HORIZONTAL, 5 );
	fgOutFile->Add( m_fpData      , 0, wxALL|wxALIGN_CENTER_VERTICAL|wxEXPAND                 , 5 );
	fgOutFile->Add( staTxtVenFile , 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_CENTER_HORIZONTAL, 5 );
	fgOutFile->Add( m_fpVen       , 0, wxALL|wxALIGN_CENTER_VERTICAL|wxEXPAND                 , 5 );
	fgOutFile->Add( staTxtIpf     , 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_CENTER_HORIZONTAL, 5 );
	fgOutFile->Add( m_fpIpf       , 0, wxALL|wxALIGN_CENTER_VERTICAL|wxEXPAND                 , 5 );
	fgOutFile->Add( staTxtCi      , 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_CENTER_HORIZONTAL, 5 );
	fgOutFile->Add( m_fpCi        , 0, wxALL|wxALIGN_CENTER_VERTICAL|wxEXPAND                 , 5 );

	this->SetSizer( bIdxMisc );
	this->Layout();

	// Connect Events
	m_txtBw   ->Connect( wxEVT_COMMAND_TEXT_UPDATED      , wxCommandEventHandler      ( IdxParamPanel::BandwidthChanged ), NULL, this );
	m_btnBwPrv->Connect( wxEVT_COMMAND_BUTTON_CLICKED    , wxCommandEventHandler      ( IdxParamPanel::OnPreview        ), NULL, this );
	m_fpData  ->Connect( wxEVT_COMMAND_FILEPICKER_CHANGED, wxFileDirPickerEventHandler( IdxParamPanel::DataFileChanged  ), NULL, this );
}

IdxParamPanel::~IdxParamPanel() {
	// Disconnect Events
	m_txtBw   ->Disconnect( wxEVT_COMMAND_TEXT_UPDATED      , wxCommandEventHandler      ( IdxParamPanel::BandwidthChanged ), NULL, this );
	m_btnBwPrv->Disconnect( wxEVT_COMMAND_BUTTON_CLICKED    , wxCommandEventHandler      ( IdxParamPanel::OnPreview        ), NULL, this );
	m_fpData  ->Disconnect( wxEVT_COMMAND_FILEPICKER_CHANGED, wxFileDirPickerEventHandler( IdxParamPanel::DataFileChanged  ), NULL, this );
}

#endif//_IDX_PARAM_PAN_H_
