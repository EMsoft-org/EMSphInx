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

#ifndef _PAT_LOAD_PAN_H_
#define _PAT_LOAD_PAN_H_

#include <wx/artprov.h>
#include <wx/xrc/xmlres.h>
#include <wx/string.h>
#include <wx/filepicker.h>
#include <wx/gdicmn.h>
#include <wx/font.h>
#include <wx/colour.h>
#include <wx/settings.h>
#include <wx/sizer.h>
#include <wx/statbox.h>
#include <wx/stattext.h>
#include <wx/textctrl.h>
#include <wx/bitmap.h>
#include <wx/image.h>
#include <wx/icon.h>
#include <wx/button.h>
#include <wx/checkbox.h>
#include <wx/panel.h>
#include <wx/valnum.h>
#include <wx/progdlg.h>

#include "wx/ValidityPanel.h"

///////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
/// Class PatternLoadPanelBase
///////////////////////////////////////////////////////////////////////////////
class PatternLoadPanelBase : public ValidityPanel {
		wxFilePickerCtrl* m_filePicker;
		wxTextCtrl      * m_txtW      ;
		wxTextCtrl      * m_txtH      ;
		wxTextCtrl      * m_txtBit    ;
		wxTextCtrl      * m_txtNum    ;
		wxTextCtrl      * m_txtPrvCnt ;
		wxButton        * m_btnPrv    ;
		wxTextCtrl      * m_mskRad    ;
		wxCheckBox      * m_chkBckg   ;
		wxTextCtrl      * m_txtNReg   ;

		wxIntegerValidator<int> valSize  ;
		wxIntegerValidator<int> valPrvNum;
		wxIntegerValidator<int> valCirc  ;
		wxIntegerValidator<int> valAhe   ;

	protected:

		// Virtual event handlers, overide them in your derived class
		virtual void FileChanged( wxFileDirPickerEvent& event ) { testValid(); event.Skip(); }
		virtual void WidthChanged( wxCommandEvent& event ) { testValid(); event.Skip(); }
		virtual void HeightChanged( wxCommandEvent& event ) { testValid(); event.Skip(); }
		virtual void NumPrvChanged( wxCommandEvent& event ) { testValid(); event.Skip(); }
		virtual void DoPreview( wxCommandEvent& event ) { event.Skip(); }
		virtual void CircChanged( wxCommandEvent& event ) { testValid(); event.Skip(); }
		virtual void AheChanged( wxCommandEvent& event ) { testValid(); event.Skip(); }


		void EnableWH(const bool enable) {m_txtW->SetEditable(enable); m_txtH->SetEditable(enable);}
		bool IsEnabledWH() const {return m_txtW->IsEnabled() && m_txtH->IsEnabled();}
		void EnablePrv(const bool enable) {m_btnPrv->Enable(enable);}
		void ClearInfo();

	public:

		PatternLoadPanelBase( wxWindow* parent, wxWindowID id = wxID_ANY, const wxPoint& pos = wxDefaultPosition, const wxSize& size = wxSize( 400,400 ), long style = wxTAB_TRAVERSAL, const wxString& name = wxEmptyString );
		~PatternLoadPanelBase();

		//@brief: sanity check the current state
		//@return: true if the values parsed from the panel are reasonable, false otherwise
		//@note  : checks for has a file, has detector sizes, and has an AHE value
		std::string validMsg() const {
			if(!wxFileExists(getFile())) return "pattern file doesn't exist";
			if(!hasW()) return "pattern width empty";
			if(!hasH()) return "pattern height empty";
			if(0 == getNum()) return "0 patterns";
			if(!hasNprv()) return "preview count empty";
			if(!hasCirc()) return "mask radius empty";
			if(!hasNreg()) return "AHE tiles empty";
			return "";
		}

		//@brief: get values from editable fields
		wxString getFile() const {return m_filePicker->GetPath();}
		long     getW   () const {long v; m_txtW     ->GetLineText(0).ToLong(&v); return v;}
		long     getH   () const {long v; m_txtH     ->GetLineText(0).ToLong(&v); return v;}
		long     getBit () const {long v; m_txtBit   ->GetLineText(0).ToLong(&v); return v;}
		long     getNum () const {long v; m_txtNum   ->GetLineText(0).ToLong(&v); return v;}
		long     getNprv() const {long v; m_txtPrvCnt->GetLineText(0).ToLong(&v); return v;}
		long     getCirc() const {long v; m_mskRad   ->GetLineText(0).ToLong(&v); return v;}
		bool     getBckg() const {return m_chkBckg->IsChecked();}
		long     getNreg() const {long v; m_txtNReg  ->GetLineText(0).ToLong(&v); return v;}

		//@brief: check if text fields are empty
		bool hasW   () const {return !m_txtW     ->GetLineText(0).IsEmpty();}
		bool hasH   () const {return !m_txtH     ->GetLineText(0).IsEmpty();}
		bool hasNprv() const {return !m_txtPrvCnt->GetLineText(0).IsEmpty();}
		bool hasCirc() const {return !m_mskRad   ->GetLineText(0).IsEmpty();}
		bool hasNreg() const {return !m_txtNReg  ->GetLineText(0).IsEmpty();}

		//@brief: set values of text fields
		void ClearFile() {wxFileName fn; m_filePicker->SetFileName(fn);}
		void setW   (int  w) {m_txtW     ->Clear(); m_txtW       ->operator<<(w);}
		void setH   (int  h) {m_txtH     ->Clear(); m_txtH       ->operator<<(h);}
		void setBit (int  b) {m_txtBit   ->Clear(); m_txtBit     ->operator<<(b);}
		void setNum (int  n);
		void setNprv(int  n) {m_txtPrvCnt->Clear(); m_txtPrvCnt  ->operator<<(n);}
		void setCirc(int  n) {m_mskRad   ->Clear(); m_mskRad     ->operator<<(n);}
		void setBckg(bool b) {                      m_chkBckg    ->SetValue  (b);}
		void setNreg(int  n) {m_txtNReg  ->Clear(); m_txtNReg    ->operator<<(n);}
};

#include <memory>
#include "modality/ebsd/nml.hpp"

class PatternLoadPanel : public PatternLoadPanelBase {
	std::string aux;//auxiliary string for pattern files
	size_t fileBytes;

	int imW, imH;
	std::shared_ptr< std::vector< std::vector<char> > > images;//preview images

	void DimChanged();

	protected:

		//make sure file is valid when changed
		void FileChanged( wxFileDirPickerEvent& event );
		void NumPrvChanged( wxFileDirPickerEvent& event ) {InavlidateImages(); event.Skip();}

		//update pattern number estimate if width/height aren't known from file
		void WidthChanged( wxCommandEvent& event ) {DimChanged();}
		void HeightChanged( wxCommandEvent& event ) {DimChanged();}

		//launch pattern preview on button click
		void DoPreview( wxCommandEvent& event );

	public:
		bool LoadImages();
		void InavlidateImages() {images->clear(), imW = -1, imH = -1;}
		bool ImagesValid() const {return !images->empty() && getNprv() == images->size() && getW() == imW && getH() == imH;}
		std::shared_ptr< std::vector< std::vector<char> > > GetImages() {return images;}
		PatternLoadPanel( wxWindow* parent, wxWindowID id = wxID_ANY) : PatternLoadPanelBase(parent, id) {images = std::make_shared< std::vector< std::vector<char> > >(); InavlidateImages();}

		bool validFile() const {return wxFileExists(getFile());}
		std::string getAux() const {return aux;}
		
};

///////////////////////////////////////////////////////////////////////////
// C++ code generated with wxFormBuilder (version Oct 26 2018)
// http://www.wxformbuilder.org/
//
// PLEASE DO *NOT* EDIT THIS FILE!
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////

void PatternLoadPanelBase::ClearInfo() {
	m_txtW->Clear();
	m_txtH->Clear();
	EnableWH(false);
	m_txtBit->ChangeValue("Unknown");
	m_txtNum->ChangeValue("Unknown");
	m_txtPrvCnt->ChangeValue("1");
	valPrvNum.SetMax(1);
}

void PatternLoadPanelBase::setNum (int n) {
	m_txtNum ->Clear();
	m_txtNum ->operator<<(n);
	valPrvNum.SetMax(n);
	m_txtPrvCnt->SetValidator(valPrvNum);
	if(getNprv() > n) {
		setNprv(std::max(n,1));
	} else {
		setNprv(std::min(n,100));
	}
}
		
PatternLoadPanelBase::PatternLoadPanelBase( wxWindow* parent, wxWindowID id, const wxPoint& pos, const wxSize& size, long style, const wxString& name ) : ValidityPanel( parent, id, pos, size, style, name ) {
	//split panel into 3 vertical boxes
	wxBoxSizer      * bPatInfo  = new wxBoxSizer( wxVERTICAL );
	wxStaticBoxSizer* sbPatFile = new wxStaticBoxSizer( new wxStaticBox( this, wxID_ANY, wxT("Pattern File"    ) ), wxVERTICAL );
	wxStaticBoxSizer* sbPatInfo = new wxStaticBoxSizer( new wxStaticBox( this, wxID_ANY, wxT("Pattern Info"    ) ), wxVERTICAL );
	wxStaticBoxSizer* sbImPrc   = new wxStaticBoxSizer( new wxStaticBox( this, wxID_ANY, wxT("Image Processing") ), wxVERTICAL );
	bPatInfo->Add( sbPatFile, 0, wxEXPAND, 5 );
	bPatInfo->Add( sbPatInfo, 1, wxEXPAND, 5 );
	bPatInfo->Add( sbImPrc  , 1, wxEXPAND, 5 );

	//middle box is 6x4 grid
	wxFlexGridSizer* fgPatInfo = new wxFlexGridSizer( 4, 6, 0, 0 );
	fgPatInfo->AddGrowableCol( 0 ); fgPatInfo->AddGrowableCol( 2 ); fgPatInfo->AddGrowableCol( 5 );
	fgPatInfo->AddGrowableRow( 0 ); fgPatInfo->AddGrowableRow( 1 ); fgPatInfo->AddGrowableRow( 2 ); fgPatInfo->AddGrowableRow( 3 );
	fgPatInfo->SetFlexibleDirection( wxHORIZONTAL );
	fgPatInfo->SetNonFlexibleGrowMode( wxFLEX_GROWMODE_SPECIFIED );
	sbPatInfo->Add( fgPatInfo, 1, wxEXPAND, 5 );

	//lower box is 6x2 grid
	wxFlexGridSizer* fgImPrc = new wxFlexGridSizer( 4, 6, 0, 0 );
	fgImPrc->AddGrowableCol( 0 ); fgImPrc->AddGrowableCol( 2 ); fgImPrc->AddGrowableCol( 5 );
	fgImPrc->AddGrowableRow( 0 ); fgImPrc->AddGrowableRow( 1 ); fgImPrc->AddGrowableRow( 2 ); fgImPrc->AddGrowableRow( 3 );
	fgImPrc->SetFlexibleDirection( wxHORIZONTAL );
	fgImPrc->SetNonFlexibleGrowMode( wxFLEX_GROWMODE_SPECIFIED );
	sbImPrc->Add( fgImPrc, 1, wxEXPAND, 5 );

	//build all labels
	wxStaticText* staTxtBinDetW = new wxStaticText( sbPatInfo->GetStaticBox(), wxID_ANY, wxT("Binned Detector Width" ), wxDefaultPosition, wxDefaultSize, 0            ); staTxtBinDetW->Wrap(-1);
	wxStaticText* staTxtPixW    = new wxStaticText( sbPatInfo->GetStaticBox(), wxID_ANY, wxT("pix"                   ), wxDefaultPosition, wxDefaultSize, wxALIGN_LEFT ); staTxtPixW   ->Wrap(-1);
	wxStaticText* staTxtH       = new wxStaticText( sbPatInfo->GetStaticBox(), wxID_ANY, wxT("Binned Detector Height"), wxDefaultPosition, wxDefaultSize, 0            ); staTxtH      ->Wrap(-1);
	wxStaticText* staTxtBitH    = new wxStaticText( sbPatInfo->GetStaticBox(), wxID_ANY, wxT("pix"                   ), wxDefaultPosition, wxDefaultSize, 0            ); staTxtBitH   ->Wrap(-1);
	wxStaticText* staTxtBitDpt  = new wxStaticText( sbPatInfo->GetStaticBox(), wxID_ANY, wxT("Bitdepth"              ), wxDefaultPosition, wxDefaultSize, 0            ); staTxtBitDpt ->Wrap(-1);
	wxStaticText* staTxtBits    = new wxStaticText( sbPatInfo->GetStaticBox(), wxID_ANY, wxT("bits"                  ), wxDefaultPosition, wxDefaultSize, wxALIGN_LEFT ); staTxtBits   ->Wrap(-1);
	wxStaticText* staTxtNum     = new wxStaticText( sbPatInfo->GetStaticBox(), wxID_ANY, wxT("Number"                ), wxDefaultPosition, wxDefaultSize, 0            ); staTxtNum    ->Wrap(-1);
	wxStaticText* staTxtPats    = new wxStaticText( sbPatInfo->GetStaticBox(), wxID_ANY, wxT("pats"                  ), wxDefaultPosition, wxDefaultSize, wxALIGN_LEFT ); staTxtPats   ->Wrap(-1);
	wxStaticText* staTxtPrvCnt  = new wxStaticText( sbImPrc  ->GetStaticBox(), wxID_ANY, wxT("Preview Count"         ), wxDefaultPosition, wxDefaultSize, 0            ); staTxtPrvCnt ->Wrap(-1);

	wxStaticText* staTxtCirc    = new wxStaticText( sbImPrc  ->GetStaticBox(), wxID_ANY, wxT("Circular Mask Radius"  ), wxDefaultPosition, wxDefaultSize, 0            ); staTxtCirc   ->Wrap(-1);
	wxStaticText* staTxtCircU   = new wxStaticText( sbImPrc  ->GetStaticBox(), wxID_ANY, wxT("pix"                   ), wxDefaultPosition, wxDefaultSize, 0            ); staTxtCircU  ->Wrap(-1);
	wxStaticText* staTxtNreg    = new wxStaticText( sbImPrc  ->GetStaticBox(), wxID_ANY, wxT("Adaptive Histogram Eq."), wxDefaultPosition, wxDefaultSize, 0            ); staTxtNreg   ->Wrap(-1);
	wxStaticText* staTxtNregU   = new wxStaticText( sbImPrc  ->GetStaticBox(), wxID_ANY, wxT("tiles"                 ), wxDefaultPosition, wxDefaultSize, 0            ); staTxtNregU  ->Wrap(-1);

	//set validator ranges
	valSize  .SetRange( 1, 8192);
	valPrvNum.SetRange( 1, 1   );
	valCirc  .SetRange(-1, 4096);
	valAhe   .SetRange( 0, 32  );

	//build elements for pattern file box and assemble
	wxString wildCards = wxT("EBSD Pattern Files (*.h5, *.upx, *.ebsp, or *.data)|*.h5;*.hdf;*.hdf5;*.up1;*.up2;*.ebsp;*.data");
	m_filePicker = new wxFilePickerCtrl( sbPatFile->GetStaticBox(), wxID_ANY, wxEmptyString, wxT("Select a file"), wildCards, wxDefaultPosition, wxDefaultSize, wxFLP_DEFAULT_STYLE|wxFLP_SMALL );
	sbPatFile->Add( m_filePicker, 1, wxALL|wxEXPAND, 5 );

	//build elements for pattern info box
	m_txtW   = new wxTextCtrl( sbPatInfo->GetStaticBox(), wxID_ANY, wxT(""       ), wxDefaultPosition, wxDefaultSize, wxTE_CENTER              , valSize );
	m_txtH   = new wxTextCtrl( sbPatInfo->GetStaticBox(), wxID_ANY, wxT(""       ), wxDefaultPosition, wxDefaultSize, wxTE_CENTER              , valSize );
	m_txtBit = new wxTextCtrl( sbPatInfo->GetStaticBox(), wxID_ANY, wxT("Unknown"), wxDefaultPosition, wxDefaultSize, wxTE_CENTER|wxTE_READONLY          );
	m_txtNum = new wxTextCtrl( sbPatInfo->GetStaticBox(), wxID_ANY, wxT("Unknown"), wxDefaultPosition, wxDefaultSize, wxTE_CENTER|wxTE_READONLY          );

	//build elements for image processing box
	m_txtPrvCnt = new wxTextCtrl( sbImPrc->GetStaticBox(), wxID_ANY, wxT("1"                  ), wxDefaultPosition, wxDefaultSize, wxTE_CENTER|wxTE_PROCESS_ENTER, valPrvNum );
	m_btnPrv    = new wxButton  ( sbImPrc->GetStaticBox(), wxID_ANY, wxT("Preview..."         ), wxDefaultPosition, wxDefaultSize, 0                                         );
	m_mskRad    = new wxTextCtrl( sbImPrc->GetStaticBox(), wxID_ANY, wxT("-1"                 ), wxDefaultPosition, wxDefaultSize, wxTE_CENTER                   , valCirc   );
	m_chkBckg   = new wxCheckBox( sbImPrc->GetStaticBox(), wxID_ANY, wxT("Gaussian Background"), wxDefaultPosition, wxDefaultSize, 0                                         );
	m_txtNReg   = new wxTextCtrl( sbImPrc->GetStaticBox(), wxID_ANY, wxT("0"                  ), wxDefaultPosition, wxDefaultSize, wxTE_CENTER                   , valAhe    );
	m_btnPrv->Enable(false);

	//assemble pattern info grid
	fgPatInfo->Add( 0            , 0, 1    , wxEXPAND                                        , 5 );
	fgPatInfo->Add( staTxtBinDetW, 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_CENTER_HORIZONTAL, 5 );
	fgPatInfo->Add( 0            , 0, 1    , wxEXPAND                                        , 5 );
	fgPatInfo->Add( m_txtW       , 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT            , 5 );
	fgPatInfo->Add( staTxtPixW   , 0, wxALL|wxALIGN_CENTER_VERTICAL                          , 5 );
	fgPatInfo->Add( 0            , 0, 1    , wxEXPAND                                        , 5 );

	fgPatInfo->Add( 0            , 0, 1    , wxEXPAND                                        , 5 );
	fgPatInfo->Add( staTxtH      , 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_CENTER_HORIZONTAL, 5 );
	fgPatInfo->Add( 0            , 0, 1    , wxEXPAND                                        , 5 );
	fgPatInfo->Add( m_txtH       , 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT            , 5 );
	fgPatInfo->Add( staTxtBitH   , 0, wxALL|wxALIGN_CENTER_VERTICAL                          , 5 );
	fgPatInfo->Add( 0            , 0, 1    , wxEXPAND                                        , 5 );

	fgPatInfo->Add( 0            , 0, 1    , wxEXPAND                                        , 5 );
	fgPatInfo->Add( staTxtBitDpt , 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_CENTER_HORIZONTAL, 5 );
	fgPatInfo->Add( 0            , 0, 1    , wxEXPAND                                        , 5 );
	fgPatInfo->Add( m_txtBit     , 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT            , 5 );
	fgPatInfo->Add( staTxtBits   , 0, wxALL|wxALIGN_CENTER_VERTICAL                          , 5 );
	fgPatInfo->Add( 0            , 0, 1    , wxEXPAND                                        , 5 );

	fgPatInfo->Add( 0            , 0, 1    , wxEXPAND                                        , 5 );
	fgPatInfo->Add( staTxtNum    , 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_CENTER_HORIZONTAL, 5 );
	fgPatInfo->Add( 0            , 0, 1    , wxEXPAND                                        , 5 );
	fgPatInfo->Add( m_txtNum     , 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT            , 5 );
	fgPatInfo->Add( staTxtPats   , 0, wxALL|wxALIGN_CENTER_VERTICAL                          , 5 );
	fgPatInfo->Add( 0            , 0, 1    , wxEXPAND                                        , 5 );

	//assemble image processing grid
	fgImPrc  ->Add( 0            , 0       , 1, wxEXPAND                                     , 5 );
	fgImPrc  ->Add( staTxtPrvCnt , 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_CENTER_HORIZONTAL, 5 );
	fgImPrc  ->Add( 0            , 0       , 1, wxEXPAND                                     , 5 );
	fgImPrc  ->Add( m_txtPrvCnt  , 0, wxALL|wxALIGN_CENTER_VERTICAL                          , 5 );
	fgImPrc  ->Add( m_btnPrv     , 0, wxALL|wxALIGN_CENTER_VERTICAL|wxEXPAND                 , 5 );
	fgImPrc  ->Add( 0            , 0       , 1, wxEXPAND                                     , 5 );

	fgImPrc  ->Add( 0            , 0       , 1, wxEXPAND                                     , 5 );
	fgImPrc  ->Add( staTxtCirc   , 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_CENTER_HORIZONTAL, 5 );
	fgImPrc  ->Add( 0            , 0       , 1, wxEXPAND                                     , 5 );
	fgImPrc  ->Add( m_mskRad     , 0, wxALL|wxALIGN_CENTER_VERTICAL                          , 5 );
	fgImPrc  ->Add( staTxtCircU  , 0, wxALL|wxALIGN_CENTER_VERTICAL|wxEXPAND                 , 5 );
	fgImPrc  ->Add( 0            , 0       , 1, wxEXPAND                                     , 5 );

	fgImPrc  ->Add( 0            , 0       , 1, wxEXPAND                                     , 5 );
	fgImPrc  ->Add( 0            , 0       , 1, wxEXPAND                                     , 5 );
	fgImPrc  ->Add( 0            , 0       , 1, wxEXPAND                                     , 5 );
	fgImPrc  ->Add( m_chkBckg    , 0, wxALL|wxALIGN_CENTER_VERTICAL                          , 5 );
	fgImPrc  ->Add( 0            , 0       , 1, wxEXPAND                                     , 5 );
	fgImPrc  ->Add( 0            , 0       , 1, wxEXPAND                                     , 5 );

	fgImPrc  ->Add( 0            , 0       , 1, wxEXPAND                                     , 5 );
	fgImPrc  ->Add( staTxtNreg   , 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_CENTER_HORIZONTAL, 5 );
	fgImPrc  ->Add( 0            , 0       , 1, wxEXPAND                                     , 5 );
	fgImPrc  ->Add( m_txtNReg    , 0, wxALL|wxALIGN_CENTER_VERTICAL                          , 5 );
	fgImPrc  ->Add( staTxtNregU  , 0, wxALL|wxALIGN_CENTER_VERTICAL|wxEXPAND                 , 5 );
	fgImPrc  ->Add( 0            , 0       , 1, wxEXPAND                                     , 5 );

	this->SetSizer( bPatInfo );
	this->Layout();

	// Connect Events
	m_filePicker->Connect( wxEVT_COMMAND_FILEPICKER_CHANGED, wxFileDirPickerEventHandler( PatternLoadPanelBase::FileChanged   ), NULL, this );
	m_txtW      ->Connect( wxEVT_COMMAND_TEXT_UPDATED      , wxCommandEventHandler      ( PatternLoadPanelBase::WidthChanged  ), NULL, this );
	m_txtH      ->Connect( wxEVT_COMMAND_TEXT_UPDATED      , wxCommandEventHandler      ( PatternLoadPanelBase::HeightChanged ), NULL, this );
	m_txtPrvCnt ->Connect( wxEVT_COMMAND_TEXT_UPDATED      , wxCommandEventHandler      ( PatternLoadPanelBase::NumPrvChanged ), NULL, this );
	m_txtPrvCnt ->Connect( wxEVT_COMMAND_TEXT_ENTER        , wxCommandEventHandler      ( PatternLoadPanelBase::DoPreview     ), NULL, this );
	m_btnPrv    ->Connect( wxEVT_COMMAND_BUTTON_CLICKED    , wxCommandEventHandler      ( PatternLoadPanelBase::DoPreview     ), NULL, this );
	m_mskRad    ->Connect( wxEVT_COMMAND_TEXT_UPDATED      , wxCommandEventHandler      ( PatternLoadPanelBase::CircChanged   ), NULL, this );
	m_txtNReg   ->Connect( wxEVT_COMMAND_TEXT_UPDATED      , wxCommandEventHandler      ( PatternLoadPanelBase::AheChanged    ), NULL, this );
}

PatternLoadPanelBase::~PatternLoadPanelBase() {
	// Disconnect Events
	m_filePicker->Disconnect( wxEVT_COMMAND_FILEPICKER_CHANGED, wxFileDirPickerEventHandler( PatternLoadPanelBase::FileChanged   ), NULL, this );
	m_txtW      ->Disconnect( wxEVT_COMMAND_TEXT_UPDATED      , wxCommandEventHandler      ( PatternLoadPanelBase::WidthChanged  ), NULL, this );
	m_txtH      ->Disconnect( wxEVT_COMMAND_TEXT_UPDATED      , wxCommandEventHandler      ( PatternLoadPanelBase::HeightChanged ), NULL, this );
	m_txtPrvCnt ->Disconnect( wxEVT_COMMAND_TEXT_UPDATED      , wxCommandEventHandler      ( PatternLoadPanelBase::NumPrvChanged ), NULL, this );
	m_txtPrvCnt ->Disconnect( wxEVT_COMMAND_TEXT_ENTER        , wxCommandEventHandler      ( PatternLoadPanelBase::DoPreview     ), NULL, this );
	m_btnPrv    ->Disconnect( wxEVT_COMMAND_BUTTON_CLICKED    , wxCommandEventHandler      ( PatternLoadPanelBase::DoPreview     ), NULL, this );
	m_mskRad    ->Disconnect( wxEVT_COMMAND_TEXT_UPDATED      , wxCommandEventHandler      ( PatternLoadPanelBase::CircChanged   ), NULL, this );
	m_txtNReg   ->Disconnect( wxEVT_COMMAND_TEXT_UPDATED      , wxCommandEventHandler      ( PatternLoadPanelBase::AheChanged    ), NULL, this );

}

////////////////////////////////////////////////////////////////////

#include <thread>
#include <atomic>

#include <wx/msgdlg.h>
#include <wx/choicdlg.h>
#include "PatternPreviewPanel.h"
#include "SinglePanelDialog.h"

void PatternLoadPanel::DimChanged() {
	if(32 == getBit() && IsEnabledWH()) {//a *.data file is currently loaded
		size_t bytes = 4 * getW() * getH();//patterns per byte
		if(!(hasW() && hasH())) bytes = 0;
		size_t numPat = 0;
		if(bytes > 0) {
			numPat = fileBytes / bytes;
		}
		setNum(numPat);
		InavlidateImages();
		testValid();
	}
}

#include "modality/ebsd/pattern.hpp"

void PatternLoadPanel::FileChanged( wxFileDirPickerEvent& event ) {
	EnablePrv(false);
	InavlidateImages();
	ClearInfo();
	if(wxFileExists(getFile())) {
		//get aux path for hdf files
		if(H5::H5File::isHdf5(getFile().ToStdString().c_str())) {
			std::vector<std::string> auxPaths = emsphinx::ebsd::PatternFile::SearchH5(getFile().ToStdString());
			std::sort(auxPaths.begin(), auxPaths.end());

			//now we have all viable datasets, make sure we have at least 1 choice
			if(auxPaths.empty()) {
				ClearFile();
				wxMessageDialog msgDlg(this, "No suitable datasets found", "Error");
				msgDlg.ShowModal();
				return;
			}

			//select dataset if multiple are found
			int idx = 0;//first choice
			if(auxPaths.size() > 1) {//there are multiple choices
				wxArrayString choices;
				for(const std::string& str : auxPaths) choices.Add(str);
				wxSingleChoiceDialog dSetSelectDlg(this, wxT("Select Pattern Dataset"), "caption", choices);//2 empty strings are defaultDir and defaultFile
				if(dSetSelectDlg.ShowModal() != wxID_OK) return;//user didn't pick a dataset
				idx = dSetSelectDlg.GetSelection();//get the selected number
			}
			aux = auxPaths[idx];
		} else {
			aux = "";
		}

		//get dimensions
		int w, h;
		uint64_t num;
		emsphinx::ImageSource::Bits bits;
		emsphinx::ebsd::PatternFile::GetFileDims(getFile().ToStdString(), w, h, bits, num, aux);

		//most pattern files have everything we need
		if(w > 0 && h > 0 && num > 0) {
			switch(bits) {
				case emsphinx::ImageSource::Bits::U8 : setBit( 8); break;
				case emsphinx::ImageSource::Bits::U16: setBit(16); break;
				case emsphinx::ImageSource::Bits::F32: setBit(32); break;
				case emsphinx::ImageSource::Bits::UNK:
					ClearFile();
					wxMessageDialog msgDlg(this, "Unknown Pattern Bitdepth", "Error");
					msgDlg.ShowModal();
					return;
			}
			setW(w); setH(h); setNum(num);
			EnableWH(false); EnablePrv(true);
		} else if(-1 == w && -1 == h && emsphinx::ImageSource::Bits::F32 == bits) {//*.data files have unknown size/count
			fileBytes = num;
			EnableWH(true); setBit(32); setNum(0);
		} else {
			ClearFile();
			wxMessageDialog msgDlg(this, "Unknown Pattern Bitdepth", "Error");
			msgDlg.ShowModal();
			return;
		}
		event.Skip();//keep passing
	}
	testValid();
}

//launch pattern preview on button click
bool PatternLoadPanel::LoadImages() {
	if(ImagesValid()) return true;//already up to date

	//build the pattern file
	std::shared_ptr<emsphinx::ebsd::PatternFile> pat;
	wxProgressDialog dlg("Loading Patterns", "Loading Preview Patterns", 1, this, wxPD_CAN_ABORT | wxPD_APP_MODAL | wxPD_AUTO_HIDE | wxPD_ELAPSED_TIME | wxPD_REMAINING_TIME);
	dlg.Show();
	if(!dlg.Pulse()) return false;//we don't know how long setting up the input file will take (but it should be pretty fast for most types)
	try {
		//if this takes too long we should put it in a thread and call Pulse periodically to keep program responsive
		pat = emsphinx::ebsd::PatternFile::Read(getFile().ToStdString(), aux, getW(), getH());
		if(NULL == pat.get()) throw std::runtime_error("failed to read patterns (nullptr)\n");
	} catch (std::exception& e) {
		wxMessageDialog msgDlg(this, e.what(), "Error Reading Patterns");
		msgDlg.ShowModal();
		ClearInfo();
		ClearFile();
		return false;
	}

	//initiailize loading
	size_t numLoad = getNprv();
	images->clear();
	images->reserve(numLoad);
	dlg.SetRange(pat->numPat());

	//start a thread actually load the patterns
	std::atomic<size_t> prg;//counter for current progress
	prg.store(0);
	std::atomic_flag flg;//keep working flag
	flg.test_and_set();
	std::thread thd = std::thread([pat,this,numLoad,&prg,&flg](){
		size_t idxNext = 0;
		const double spacing = double(pat->numPat()) / numLoad;
		std::vector<char> buff(pat->imBytes());
		std::vector<char> buff8(pat->imBytes(), 0);//fill with 0 for emsphinx::ImageSource::Bits::UNK
		for(size_t i = 0; i < pat->numPat(); i++) {
			if(!flg.test_and_set()) {//was a cancel requested?
				flg.clear();//reclear flag
				break;//stop
			}
			prg.store(i);//store progress
			pat->extract(buff.data(), 1);//read a pattern
			if(i == idxNext) {
				switch(pat->pixelType()) {
					case emsphinx::ImageSource::Bits::UNK: break;
					case emsphinx::ImageSource::Bits::U8 : buff8.swap(buff); break;
					case emsphinx::ImageSource::Bits::U16: {
						uint16_t* ptr = (uint16_t*)buff.data();
						std::pair<uint16_t*, uint16_t*> minMax = std::minmax_element(ptr, ptr + pat->numPix());//compute minimum and maximum quality
						const double scale = double(255) / (*minMax.second - *minMax.first);
						std::transform(ptr, ptr + pat->numPix(), (uint8_t*)buff8.data(), [&](const uint16_t& v){return (uint8_t)std::round(scale * (v - (*minMax.first)));});
					} break;
					case emsphinx::ImageSource::Bits::F32: {
						float* ptr = (float*)buff.data();
						std::pair<float*, float*> minMax = std::minmax_element(ptr, ptr + pat->numPix());//compute minimum and maximum quality
						const float scale = 255.0f / (*minMax.second - *minMax.first);
						std::transform(ptr, ptr + pat->numPix(), (uint8_t*)buff8.data(), [&](const float& v){return (uint8_t)std::round(scale * (v - (*minMax.first)));});
					} break;
				}
				images->push_back(buff8);
				idxNext = (size_t) std::round(spacing * images->size());
			}
		}

		//flag completion
		prg.store(pat->numPat());
	});

	//wait for loading to finish
	size_t curPrg = prg.load();
	while(curPrg != pat->numPat()) {
		std::this_thread::sleep_for(std::chrono::milliseconds(30));//~30 fps update
		curPrg = prg.load();//get current progress
		if(!dlg.Update(curPrg)) {
			flg.clear();//tell worker to stop
			thd.join();
			return false;
		}
	}
	thd.join();

	//update
	imW = pat->width ();
	imH = pat->height();
	return true;
}

void PatternLoadPanel::DoPreview( wxCommandEvent& event ) {
	if(!LoadImages()) return;
	SinglePanelDialog<PatternPreviewPanel> dlg(this);
	dlg.getPanel()->SetImages(images, imW, imH);
	dlg.getPanel()->SetCirc(getCirc());
	dlg.getPanel()->SetBckg(getBckg());
	dlg.getPanel()->SetNreg(getNreg());
	dlg.getPanel()->ForceRedraw();
	if( wxID_OK == dlg.ShowModal() ) {
		setCirc(dlg.getPanel()->GetCirc());
		setBckg(dlg.getPanel()->GetBckg());
		setNreg(dlg.getPanel()->GetNreg());
	}
}

#endif//_PAT_LOAD_PAN_H_
