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

#ifndef _SCAN_DIMS_PAN_H_
#define _SCAN_DIMS_PAN_H_

#include <wx/artprov.h>
#include <wx/xrc/xmlres.h>
#include <wx/string.h>
#include <wx/stattext.h>
#include <wx/gdicmn.h>
#include <wx/font.h>
#include <wx/colour.h>
#include <wx/settings.h>
#include <wx/textctrl.h>
#include <wx/sizer.h>
#include <wx/statbox.h>
#include <wx/choice.h>
#include <wx/bitmap.h>
#include <wx/image.h>
#include <wx/icon.h>
#include <wx/button.h>
#include <wx/panel.h>
#include <wx/valnum.h>
#include <wx/checkbox.h>

#include "wx/ValidityPanel.h"

#include "modality/ebsd/pattern.hpp"
#include "idx/roi.h"

///////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
/// Class ScanDimsPanel
///////////////////////////////////////////////////////////////////////////////
class ScanDimsPanel : public ValidityPanel
{
	private:
		int m_numPat = -1;
		std::shared_ptr< std::vector<float> > m_iq, m_ci, m_iqCalc;//existing maps (can be NULL)

		emsphinx::RoiSelection m_roi;//roi
		std::vector<char> m_mask;//current roi mask

		wxIntegerValidator      <int  > valPix;//scan size
		wxFloatingPointValidator<float> valUm ;//pixel size

		void updateChoices();

	protected:
		wxTextCtrl* m_txtW    ;
		wxTextCtrl* m_txtH    ;
		wxTextCtrl* m_txtX    ;
		wxTextCtrl* m_txtY    ;
		wxChoice  * m_cmbRoiIm;
		wxButton  * m_btnRoi  ;
		wxTextCtrl* m_txtCvg  ;
		wxButton  * m_btnClear;
		wxCheckBox* m_chkDims ;

		// Virtual event handlers, overide them in your derived class
		virtual void ScanDimsChanged( wxCommandEvent& event ) { testValid(); }
		virtual void ScanStepChanged( wxCommandEvent& event ) { testValid(); }
		virtual void DoROI          ( wxCommandEvent& event );
		virtual void ClearROI       ( wxCommandEvent& event ) { m_txtCvg->SetValue("No ROI"); m_btnClear->Enable(false); m_roi.clear(); }
		virtual void OnCheckbox     ( wxCommandEvent& event ) { testValid(); }

	public:

		ScanDimsPanel( wxWindow* parent, wxWindowID id = wxID_ANY, const wxPoint& pos = wxDefaultPosition, const wxSize& size = wxSize( 350,300 ), long style = wxTAB_TRAVERSAL, const wxString& name = wxEmptyString );
		~ScanDimsPanel();

		//@brief : get the scan width
		//@return: scan width in pixels
		int getW() const {long v; m_txtW->GetLineText(0).ToLong(&v); return v;}

		//@brief : get the scan height
		//@return: scan height in pixels
		int getH() const {long v; m_txtH->GetLineText(0).ToLong(&v); return v;}

		//@brief : get the pixel width
		//@return: pixel width in microns
		double getX() const {double v; m_txtX->GetLineText(0).ToDouble(&v); return v;}

		//@brief : get the pixel width
		//@return: pixel width in microns
		double getY() const {double v; m_txtY->GetLineText(0).ToDouble(&v); return v;}

		//@brief  : set the scan width
		//@param w: scan width in pixels
		void setW(const int w) {m_txtW->Clear(); m_txtW->operator<<(w);}

		//@brief  : set the scan height
		//@param h: scan height in pixels
		void setH(const int h) {m_txtH->Clear(); m_txtH->operator<<(h);}

		//@brief  : set the pixel width
		//@param y: pixel width in microns
		void setX(const float x) {m_txtX->Clear(); m_txtX->operator<<(x);}

		//@brief  : set the pixel width
		//@param z: pixel width in microns
		void setY(const float y) {m_txtY->Clear(); m_txtY->operator<<(y);}

		//@brief    : set existing maps (for preview) and # patterns
		//@param iq : iq map
		//@param ci : ci map
		//@param num: number of patterns
		void setMaps(std::shared_ptr< std::vector<float> > iq, std::shared_ptr< std::vector<float> > ci, int num) {m_iq = iq; m_ci = ci; m_numPat = num; updateChoices();}

		//@brief    : set self computed image quality map
		//@param iq : iq map computed directly from patterns
		void setCalcIq(std::shared_ptr< std::vector<float> > iq) {m_iqCalc = iq; updateChoices();}

		//@brief : get the selected map (or empty pointer for no selection)
		//@return: currently selected image
		wxImage getMap() const;

		//@brief : get the ROI
		//@return: roi
		emsphinx::RoiSelection getRoi() const {return m_roi;}

		//@brief: sanity check the current state
		//@return: true if the values parsed from the panel are reasonable, false otherwise
		//@note  : checks for has a file, has detector sizes, and has an AHE value
		std::string validMsg() const {
			m_btnRoi->Enable(false);
			if(m_txtW->GetLineText(0).IsEmpty()) return "scan width empty";
			if(m_txtH->GetLineText(0).IsEmpty()) return "scan height empty";
			if(m_txtX->GetLineText(0).IsEmpty()) return "x step empty";
			if(m_txtY->GetLineText(0).IsEmpty()) return "y step empty";
			int numPix = getW() * getH();
			std::ostringstream ss;
			ss << m_numPat << " patterns for " << numPix << " pixels";
			if(numPix > m_numPat) return "not enough patterns (" + ss.str() + ')';
			if(numPix < m_numPat && !m_chkDims->IsChecked()) return "too many patterns (" + ss.str() + ')';
			m_btnRoi->Enable(true);
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

#include <wx/msgdlg.h>

#include "SinglePanelDialog.h"
#include "RoiSelectionPanel.h"

void ScanDimsPanel::updateChoices() {
	//accumluate choices
	int numChoice = 0;
	wxString choices[3];
	if(NULL != m_iqCalc.get()) choices[numChoice++] = wxT("Computed IQ");//prefer computed iq over everything (most consistent)
	if(NULL == m_iq.get() ? false : !m_iq->empty()) choices[numChoice++] = wxT("Existing IQ");//prefer existing iq over existing ci (nicer image)
	if(NULL == m_ci.get() ? false : !m_ci->empty()) choices[numChoice++] = wxT("Existing CI");

	//now build box and set choice 
	m_cmbRoiIm->Clear();
	for(int i = 0; i < numChoice; i++) m_cmbRoiIm->Append(choices[i]);
	m_cmbRoiIm->SetSelection(0 == numChoice ? wxNOT_FOUND : 0);
	// m_btnRoi->Enable(numChoice > 0);
}

//@brief : get the selected map (or empty pointer for no selection)
//@return: currently selected image
wxImage ScanDimsPanel::getMap() const {
	//get float image
	std::shared_ptr< std::vector<float> > ptr;
	wxString sel = m_cmbRoiIm->GetString(m_cmbRoiIm->GetSelection());
	if     (wxString("Existing CI") == sel) ptr = m_ci;
	else if(wxString("Existing IQ") == sel) ptr = m_iq;
	else if(wxString("Computed IQ") == sel) ptr = m_iqCalc;

	//convert to 8 bit (this needs to be moved to another file)
	std::vector<char> buff(m_numPat);
	if(NULL != ptr.get())  {
		std::pair<float*, float*> minMax = std::minmax_element(ptr->data(), ptr->data() + buff.size());//compute minimum and maximum quality
		const float scale = 255.0f / (*minMax.second - *minMax.first);
		std::transform(ptr->data(), ptr->data() + buff.size(), (uint8_t*)buff.data(), [&](const float& v){return (uint8_t)std::round(scale * (v - (*minMax.first)));});
	}

	//wrap as image (likewise)
	size_t w = getW();
	size_t h = getH();
	wxImage im(getW(), getH());
	unsigned char* alpha = (unsigned char*)malloc(buff.size());//SetAlpha requires this...
	std::fill(alpha, alpha + buff.size(), 0xFF);
	im.SetAlpha(alpha, false);//false -> image takes ownership and will free
	for(size_t j = 0; j < h; j++) {
		for(size_t i = 0; i < w; i++) {
			char c = buff[j * w + i];
			im.SetRGB(i, j, c, c, c);
		}
	}
	return im;
}

void ScanDimsPanel::DoROI( wxCommandEvent& event ) {
	if(0 == m_cmbRoiIm->GetCount()) {
		wxMessageDialog msgDlg(this, "No suitable maps found, try selecting a pattern file with more data (*.h5 scan or *.upX/*.ebsp with associated *.ang/*.ctf) or computing image quailty during pattern load", "Error");
		msgDlg.ShowModal();
		return;
	}

	SinglePanelDialog<RoiSelectionPanel> dlg(this, wxID_ANY, "Select ROI");
	wxImage im = getMap();
	dlg.getPanel()->setImage(im);
	dlg.getPanel()->setRoi(m_roi);
	if(wxID_OK == dlg.ShowModal()) {//this currently crashes...
		m_roi = dlg.getPanel()->getRoi();
		std::vector<char> mask = m_roi.buildMask(getW(), getH());
		size_t num0 = std::count_if(mask.begin(), mask.end(), [](const char&c){return c == 0;});
		double frac = 1.0 - double(num0) / (mask.size());
		m_txtCvg->ChangeValue(wxString::Format(wxT("%0.1f"), frac * 100));
		m_btnClear->Enable(true);
	}
}

ScanDimsPanel::ScanDimsPanel( wxWindow* parent, wxWindowID id, const wxPoint& pos, const wxSize& size, long style, const wxString& name ) : ValidityPanel( parent, id, pos, size, style, name ) {
	//split panel into 2 vertical boxes
	wxBoxSizer      * bScnDm   = new wxBoxSizer( wxVERTICAL );
	wxStaticBoxSizer* sbScnDim = new wxStaticBoxSizer( new wxStaticBox( this, wxID_ANY, wxT("Scan Dimensions"   ) ), wxVERTICAL );
	wxStaticBoxSizer* sbRoi    = new wxStaticBoxSizer( new wxStaticBox( this, wxID_ANY, wxT("Region of Interest") ), wxVERTICAL );
	bScnDm->Add( sbScnDim, 2, wxEXPAND, 5 );
	bScnDm->Add( sbRoi   , 1, wxEXPAND, 5 );

	//top pox is 6x4 grid
	wxFlexGridSizer* fgScnDim = new wxFlexGridSizer( 4, 6, 0, 0 );
	fgScnDim->AddGrowableCol( 0 ); fgScnDim->AddGrowableCol( 2 ); fgScnDim->AddGrowableCol( 5 );
	fgScnDim->AddGrowableRow( 0 ); fgScnDim->AddGrowableRow( 1 ); fgScnDim->AddGrowableRow( 2 ); fgScnDim->AddGrowableRow( 3 );
	fgScnDim->SetFlexibleDirection( wxHORIZONTAL );
	fgScnDim->SetNonFlexibleGrowMode( wxFLEX_GROWMODE_SPECIFIED );
	sbScnDim->Add( fgScnDim, 1, wxEXPAND, 5 );

	//bottom box is 6x2 grid
	wxFlexGridSizer* fgRoi = new wxFlexGridSizer( 2, 6, 0, 0 );
	fgRoi->AddGrowableCol( 0 ); fgRoi->AddGrowableCol( 2 ); fgRoi->AddGrowableCol( 5 );
	fgRoi->AddGrowableRow( 0 ); fgRoi->AddGrowableRow( 1 );
	fgRoi->SetFlexibleDirection( wxHORIZONTAL );
	fgRoi->SetNonFlexibleGrowMode( wxFLEX_GROWMODE_SPECIFIED );
	sbRoi->Add( fgRoi, 1, wxEXPAND, 5 );

	//sub grid for # and clear ubbton
	wxFlexGridSizer* fgpctClr = new wxFlexGridSizer( 1, 2, 0, 0 );
	fgpctClr->AddGrowableCol( 0 ); fgpctClr->AddGrowableCol( 1 );
	fgpctClr->AddGrowableRow( 0 );
	fgpctClr->SetFlexibleDirection( wxBOTH );
	fgpctClr->SetNonFlexibleGrowMode( wxFLEX_GROWMODE_SPECIFIED );

	//build all labels
	wxStaticText* staTxtScnW     = new wxStaticText( sbScnDim->GetStaticBox(), wxID_ANY, wxT("Scan Width" ), wxDefaultPosition, wxDefaultSize, 0            ); staTxtScnW    ->Wrap(-1);
	wxStaticText* staTxtPixW     = new wxStaticText( sbScnDim->GetStaticBox(), wxID_ANY, wxT("pix"        ), wxDefaultPosition, wxDefaultSize, wxALIGN_LEFT ); staTxtPixW    ->Wrap(-1);
	wxStaticText* staTxtScnH     = new wxStaticText( sbScnDim->GetStaticBox(), wxID_ANY, wxT("Scan Height"), wxDefaultPosition, wxDefaultSize, 0            ); staTxtScnH    ->Wrap(-1);
	wxStaticText* staTxtPixH     = new wxStaticText( sbScnDim->GetStaticBox(), wxID_ANY, wxT("pix"        ), wxDefaultPosition, wxDefaultSize, 0            ); staTxtPixH    ->Wrap(-1);
	wxStaticText* staTxtStpX     = new wxStaticText( sbScnDim->GetStaticBox(), wxID_ANY, wxT("X Step"     ), wxDefaultPosition, wxDefaultSize, 0            ); staTxtStpX    ->Wrap(-1);
	wxStaticText* staTxtUmX      = new wxStaticText( sbScnDim->GetStaticBox(), wxID_ANY, wxT("um"         ), wxDefaultPosition, wxDefaultSize, wxALIGN_LEFT ); staTxtUmX     ->Wrap(-1);
	wxStaticText* staTxtScnStepY = new wxStaticText( sbScnDim->GetStaticBox(), wxID_ANY, wxT("Y Step"     ), wxDefaultPosition, wxDefaultSize, 0            ); staTxtScnStepY->Wrap(-1);
	wxStaticText* staTxtUmY      = new wxStaticText( sbScnDim->GetStaticBox(), wxID_ANY, wxT("um"         ), wxDefaultPosition, wxDefaultSize, wxALIGN_LEFT ); staTxtUmY     ->Wrap(-1);
	wxStaticText* staTxtIm       = new wxStaticText( sbRoi   ->GetStaticBox(), wxID_ANY, wxT("Image"      ), wxDefaultPosition, wxDefaultSize, 0            ); staTxtIm      ->Wrap(-1);
	wxStaticText* staTxtCvg      = new wxStaticText( sbRoi   ->GetStaticBox(), wxID_ANY, wxT("Coverage"   ), wxDefaultPosition, wxDefaultSize, 0            ); staTxtCvg     ->Wrap(-1);
	wxStaticText* staTxtCvgU     = new wxStaticText( sbRoi   ->GetStaticBox(), wxID_ANY, wxT("%"          ), wxDefaultPosition, wxDefaultSize, 0            ); staTxtCvgU    ->Wrap(-1);

	//set text validator ranges/precisions
	valPix.SetRange( 1    , 5000 );
	valUm .SetRange( 0.01f, 100  );
	valUm.SetPrecision(2);//max 6 for float

	//build elements for scan dimensions box
	m_txtW     = new wxTextCtrl( sbScnDim->GetStaticBox(), wxID_ANY, wxT(""                     ), wxDefaultPosition, wxDefaultSize, wxTE_CENTER             , valPix         );
	m_txtH     = new wxTextCtrl( sbScnDim->GetStaticBox(), wxID_ANY, wxT(""                     ), wxDefaultPosition, wxDefaultSize, wxTE_CENTER             , valPix         );
	m_txtX     = new wxTextCtrl( sbScnDim->GetStaticBox(), wxID_ANY, wxT(""                     ), wxDefaultPosition, wxDefaultSize, wxTE_CENTER             , valUm          );
	m_txtY     = new wxTextCtrl( sbScnDim->GetStaticBox(), wxID_ANY, wxT(""                     ), wxDefaultPosition, wxDefaultSize, wxTE_CENTER             , valUm          );
	m_chkDims  = new wxCheckBox( sbScnDim->GetStaticBox(), wxID_ANY, wxT("Ignore Extra Patterns"), wxDefaultPosition, wxDefaultSize, 0                                        );
	sbScnDim->Add( m_chkDims, 0, wxALIGN_RIGHT, 5 );

	//build elements for ROI box
	wxString m_cmbRoiImChoices[] = { wxT("Compute IQ"), wxT("Existing IQ"), wxT("Existing CI") };
	int m_cmbRoiImNChoices = sizeof( m_cmbRoiImChoices ) / sizeof( wxString );
	m_cmbRoiIm = new wxChoice  ( sbRoi   ->GetStaticBox(), wxID_ANY,                       wxDefaultPosition, wxDefaultSize, 0                 , m_cmbRoiImChoices, 0 );
	m_btnRoi   = new wxButton  ( sbRoi   ->GetStaticBox(), wxID_ANY, wxT("Select ROI..."), wxDefaultPosition, wxDefaultSize,                         0                );
	m_txtCvg   = new wxTextCtrl( sbRoi   ->GetStaticBox(), wxID_ANY, wxT("No ROI"       ), wxDefaultPosition, wxDefaultSize, wxTE_CENTER|wxTE_READONLY                );
	m_btnClear = new wxButton  ( sbRoi   ->GetStaticBox(), wxID_ANY, wxT("Clear ROI"    ), wxDefaultPosition, wxDefaultSize,                         0                );
	m_cmbRoiIm->SetSelection( -1 );
	
	//assemble scan dims grid
	fgScnDim->Add( 0             , 0, 1    , wxEXPAND                                        , 5 );
	fgScnDim->Add( staTxtScnW    , 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_CENTER_HORIZONTAL, 5 );
	fgScnDim->Add( 0             , 0, 1    , wxEXPAND                                        , 5 );
	fgScnDim->Add( m_txtW        , 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT            , 5 );
	fgScnDim->Add( staTxtPixW    , 0, wxALL|wxALIGN_CENTER_VERTICAL                          , 5 );
	fgScnDim->Add( 0             , 0, 1    , wxEXPAND                                        , 5 );
	fgScnDim->Add( 0             , 0, 1    , wxEXPAND                                        , 5 );
	fgScnDim->Add( staTxtScnH    , 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_CENTER_HORIZONTAL, 5 );
	fgScnDim->Add( 0             , 0, 1    , wxEXPAND                                        , 5 );
	fgScnDim->Add( m_txtH        , 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT            , 5 );
	fgScnDim->Add( staTxtPixH    , 0, wxALL|wxALIGN_CENTER_VERTICAL                          , 5 );
	fgScnDim->Add( 0             , 0, 1    , wxEXPAND                                        , 5 );
	fgScnDim->Add( 0             , 0, 1    , wxEXPAND                                        , 5 );
	fgScnDim->Add( staTxtStpX    , 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_CENTER_HORIZONTAL, 5 );
	fgScnDim->Add( 0             , 0, 1    , wxEXPAND                                        , 5 );
	fgScnDim->Add( m_txtX        , 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT            , 5 );
	fgScnDim->Add( staTxtUmX     , 0, wxALL|wxALIGN_CENTER_VERTICAL                          , 5 );
	fgScnDim->Add( 0             , 0, 1    , wxEXPAND                                        , 5 );
	fgScnDim->Add( 0             , 0, 1    , wxEXPAND                                        , 5 );
	fgScnDim->Add( staTxtScnStepY, 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_CENTER_HORIZONTAL, 5 );
	fgScnDim->Add( 0             , 0, 1    , wxEXPAND                                        , 5 );
	fgScnDim->Add( m_txtY        , 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT            , 5 );
	fgScnDim->Add( staTxtUmY     , 0, wxALL|wxALIGN_CENTER_VERTICAL                          , 5 );
	fgScnDim->Add( 0             , 0, 1    , wxEXPAND                                        , 5 );

	//assemble %/clear roi subgrid
	fgpctClr->Add( staTxtCvgU    , 0, wxALL|wxALIGN_CENTER_VERTICAL                          , 5 );
	fgpctClr->Add( m_btnClear    , 0, wxALL|wxALIGN_CENTER_VERTICAL                          , 5 );

	//assemble roi grid
	fgRoi   ->Add( 0             , 0, 1    , wxEXPAND                                        , 5 );
	fgRoi   ->Add( staTxtIm      , 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_CENTER_HORIZONTAL, 5 );
	fgRoi   ->Add( 0             , 0, 1    , wxEXPAND                                        , 5 );
	fgRoi   ->Add( m_cmbRoiIm    , 0, wxALL|wxALIGN_CENTER_VERTICAL|wxEXPAND                 , 5 );
	fgRoi   ->Add( m_btnRoi      , 0, wxALL|wxALIGN_CENTER_VERTICAL|wxEXPAND                 , 5 );
	fgRoi   ->Add( 0             , 0, 1    , wxEXPAND                                        , 5 );
	fgRoi   ->Add( 0             , 0, 1    , wxEXPAND                                        , 5 );
	fgRoi   ->Add( staTxtCvg     , 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_CENTER_HORIZONTAL, 5 );
	fgRoi   ->Add( 0             , 0, 1    , wxEXPAND                                        , 5 );
	fgRoi   ->Add( m_txtCvg      , 0, wxALL|wxALIGN_CENTER_VERTICAL|wxEXPAND                 , 5 );
	fgRoi   ->Add( fgpctClr      , 1,                               wxEXPAND                 , 5 );
	fgRoi   ->Add( 0             , 0, 1    , wxEXPAND                                        , 5 );

	this->SetSizer( bScnDm );
	this->Layout();
	// m_btnRoi->Enable(false);//we don't have an image source yet
	m_btnClear->Enable(false);//we don't have an ROI yet

	// Connect Events
	m_txtW    ->Connect( wxEVT_COMMAND_TEXT_UPDATED    , wxCommandEventHandler( ScanDimsPanel::ScanDimsChanged ), NULL, this );
	m_txtH    ->Connect( wxEVT_COMMAND_TEXT_UPDATED    , wxCommandEventHandler( ScanDimsPanel::ScanDimsChanged ), NULL, this );
	m_txtX    ->Connect( wxEVT_COMMAND_TEXT_UPDATED    , wxCommandEventHandler( ScanDimsPanel::ScanStepChanged ), NULL, this );
	m_txtY    ->Connect( wxEVT_COMMAND_TEXT_UPDATED    , wxCommandEventHandler( ScanDimsPanel::ScanStepChanged ), NULL, this );
	m_btnRoi  ->Connect( wxEVT_COMMAND_BUTTON_CLICKED  , wxCommandEventHandler( ScanDimsPanel::DoROI           ), NULL, this );
	m_btnClear->Connect( wxEVT_COMMAND_BUTTON_CLICKED  , wxCommandEventHandler( ScanDimsPanel::ClearROI        ), NULL, this );
	m_chkDims ->Connect( wxEVT_COMMAND_CHECKBOX_CLICKED, wxCommandEventHandler( ScanDimsPanel::OnCheckbox      ), NULL, this );
}

ScanDimsPanel::~ScanDimsPanel() {
	// Disconnect Events
	m_txtW    ->Disconnect( wxEVT_COMMAND_TEXT_UPDATED    , wxCommandEventHandler( ScanDimsPanel::ScanDimsChanged ), NULL, this );
	m_txtH    ->Disconnect( wxEVT_COMMAND_TEXT_UPDATED    , wxCommandEventHandler( ScanDimsPanel::ScanDimsChanged ), NULL, this );
	m_txtX    ->Disconnect( wxEVT_COMMAND_TEXT_UPDATED    , wxCommandEventHandler( ScanDimsPanel::ScanStepChanged ), NULL, this );
	m_txtY    ->Disconnect( wxEVT_COMMAND_TEXT_UPDATED    , wxCommandEventHandler( ScanDimsPanel::ScanStepChanged ), NULL, this );
	m_btnRoi  ->Disconnect( wxEVT_COMMAND_BUTTON_CLICKED  , wxCommandEventHandler( ScanDimsPanel::DoROI           ), NULL, this );
	m_btnClear->Disconnect( wxEVT_COMMAND_BUTTON_CLICKED  , wxCommandEventHandler( ScanDimsPanel::ClearROI        ), NULL, this );
	m_chkDims ->Disconnect( wxEVT_COMMAND_CHECKBOX_CLICKED, wxCommandEventHandler( ScanDimsPanel::OnCheckbox      ), NULL, this );
}

#endif//_SCAN_DIMS_PAN_H_
