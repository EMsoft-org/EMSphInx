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

#ifndef _MP_FILT_PAN_H_
#define _MP_FILT_PAN_H_

#include <wx/string.h>
#include <wx/checkbox.h>
#include <wx/settings.h>
#include <wx/textctrl.h>
#include <wx/stattext.h>
#include <wx/choice.h>
#include <wx/sizer.h>
#include <wx/panel.h>
#include <wx/valnum.h>
#include <wx/button.h>
#include <wx/tglbtn.h>

///////////////////////////////////////////////////////////////////////////

#include "PeriodicTablePanel.h"

///////////////////////////////////////////////////////////////////////////////
/// Class MasterFileFilterPanel
///////////////////////////////////////////////////////////////////////////////
class MasterFileFilterPanel : public wxPanel
{
	private:
		wxFloatingPointValidator<float> valKv ;
		wxFloatingPointValidator<float> valTlt;
		wxIntegerValidator      <int  > valSg ;

	protected:
		wxCheckBox        * m_chkKv ;
		wxTextCtrl        * m_kvMin ;
		wxTextCtrl        * m_kvMax ;
		wxCheckBox        * m_chkTlt;
		wxTextCtrl        * m_tltMin;
		wxTextCtrl        * m_tltMax;
		wxCheckBox        * m_chkSg ;
		wxTextCtrl        * m_sgMin ;
		wxTextCtrl        * m_sgMax ;
		wxChoice          * m_laue  ;
		wxCheckBox        * m_chkEl ;
		wxToggleButton    * m_btnExt;
		wxButton          * m_btnClr;
		PeriodicTablePanel* m_prdTbl;

		void ChkKvChanged (wxCommandEvent& event) { if(!event.IsChecked()) clearKv (); }//m_chkKv  clicked
		void ChkTltChanged(wxCommandEvent& event) { if(!event.IsChecked()) clearTlt(); }//m_chkTlt clicked
		void ChkSgChanged (wxCommandEvent& event) { if(!event.IsChecked()) clearSg (); }//m_chkSg  clicked
		void ChkElChanged (wxCommandEvent& event) { if(!event.IsChecked()) clearEl (); }//m_chkSg  clicked
		void ClearPressed (wxCommandEvent& event) { clearEl (); }

		void ChoiceChanged  ( wxCommandEvent& event );
		void ElementsChanged( wxCommandEvent& event ) {m_chkEl->SetValue(ElementMask() != m_prdTbl->getMask());}

		//tick the filter use checkbox when text fields are changed
		void TxtKvMinChanged ( wxCommandEvent& event ) { m_chkKv ->SetValue(true); }
		void TxtKvMaxChanged ( wxCommandEvent& event ) { m_chkKv ->SetValue(true); }
		void TxtTltMinChanged( wxCommandEvent& event ) { m_chkTlt->SetValue(true); }
		void TxtTltMaxChanged( wxCommandEvent& event ) { m_chkTlt->SetValue(true); }
		void TxtSgMinChanged ( wxCommandEvent& event ) { m_chkSg ->SetValue(true); }
		void TxtSgMaxChanged ( wxCommandEvent& event ) { m_chkSg ->SetValue(true); }

		//@brief: clear the current acclerating voltage filter
		void clearKv();

		//@brief: clear the current tilt filter
		void clearTlt();

		//@brief: clear the current space group filter
		void clearSg();

		//@brief: clear the current elements filter
		void clearEl();

	public:
		//@brief    : set filter bounds
		//@param kv : upper and lower bounds for kv filtering
		//@param tlt: upper and lower bounds for tilt filtering
		//@param el : bitmask of elements
		//@param ext: is element bitmask exact or inclusive
		//@param sg : upper and lower bounds for space group filtering
		void setBounds(std::pair<float, float> kv, std::pair<float, float> tlt, ElementMask el, bool ext, std::pair<int  , int  > sg);

		//@brief    : get filter bounds
		//@param kv : location to write upper and lower bounds for kv filtering
		//@param tlt: location to write upper and lower bounds for tilt filtering
		//@param el : location to write bitmask of elements
		//@param ext: is element bitmask exact or inclusive
		//@param sg : location to write upper and lower bounds for space group filtering
		void getBounds(std::pair<float, float>& kv, std::pair<float, float>& tlt, ElementMask& el, bool& ext, std::pair<int  , int  >& sg) const;

		MasterFileFilterPanel( wxWindow* parent, wxWindowID id = wxID_ANY, const wxPoint& pos = wxDefaultPosition, const wxSize& size = wxSize( 695,300 ), long style = wxTAB_TRAVERSAL, const wxString& name = wxEmptyString );
		~MasterFileFilterPanel();

};


///////////////////////////////////////////////////////////////////////////
// C++ code generated with wxFormBuilder (version Oct 26 2018)
// http://www.wxformbuilder.org/
//
// PLEASE DO *NOT* EDIT THIS FILE!
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////

//@brief: clear the current acclerating voltage filter
void MasterFileFilterPanel::clearKv() {
	m_chkKv ->SetValue(false);
	m_kvMin->ChangeValue( "0.0");
	m_kvMax->ChangeValue("50.0");
}

//@brief: clear the current tilt filter
void MasterFileFilterPanel::clearTlt() {
	m_chkTlt->SetValue(false);
	m_tltMin->ChangeValue( "0.0");
	m_tltMax->ChangeValue("90.0");
}

//@brief: clear the current space group filter
void MasterFileFilterPanel::clearSg() {
	m_chkSg ->SetValue(false);
	m_sgMin ->ChangeValue(  "1");
	m_sgMax ->ChangeValue("230");
	m_laue  ->SetSelection(wxNOT_FOUND);
}

//@brief: clear the current elements filter
void MasterFileFilterPanel::clearEl() {
	m_chkEl ->SetValue(false);
	m_prdTbl->setMask(ElementMask());
}

//@brief    : set filter bounds
//@param kv : upper and lower bounds for kv filtering
//@param tlt: upper and lower bounds for tilt filtering
//@param el : bitmask of elements
//@param ext: is element bitmask exact or inclusive
//@param sg : upper and lower bounds for space group filtering
void MasterFileFilterPanel::setBounds(std::pair<float, float> kv, std::pair<float, float> tlt, ElementMask el, bool ext, std::pair<int  , int  > sg) {
	if(kv .first == kv .first) {//not nan
		m_chkKv ->SetValue(true);
		m_kvMin ->Clear(); m_kvMin ->operator<<(kv .first );
		m_kvMax ->Clear(); m_kvMax ->operator<<(kv .second);
	}

	if(tlt.first == tlt.first) {//not nan
		m_chkTlt->SetValue(true);
		m_tltMin->Clear(); m_tltMin->operator<<(tlt.first );
		m_tltMax->Clear(); m_tltMax->operator<<(tlt.second);
	}

	if(sg.first > 1 || sg.second < 230) {
		m_chkSg ->SetValue(true);
		m_sgMin ->Clear(); m_sgMin ->operator<<(sg .first );
		m_sgMax ->Clear(); m_sgMax ->operator<<(sg .second);

		//check if this range corresponds to a laue group
		if(  1 == sg.first &&    2 == sg.second) m_laue->SetSelection( 0);
		if(  3 == sg.first &&   15 == sg.second) m_laue->SetSelection( 1);
		if( 16 == sg.first &&   74 == sg.second) m_laue->SetSelection( 2);
		if( 75 == sg.first &&   88 == sg.second) m_laue->SetSelection( 3);
		if( 89 == sg.first &&  142 == sg.second) m_laue->SetSelection( 4);
		if(143 == sg.first &&  148 == sg.second) m_laue->SetSelection( 5);
		if(149 == sg.first &&  167 == sg.second) m_laue->SetSelection( 6);
		if(168 == sg.first &&  176 == sg.second) m_laue->SetSelection( 7);
		if(177 == sg.first &&  194 == sg.second) m_laue->SetSelection( 8);
		if(195 == sg.first &&  206 == sg.second) m_laue->SetSelection( 9);
		if(207 == sg.first &&  230 == sg.second) m_laue->SetSelection(10);
	}

	m_chkEl->SetValue(ElementMask() != el);
	m_prdTbl->setMask(el);
	m_btnExt->SetValue(ext);
}

//@brief    : get filter bounds
//@param kv : location to write upper and lower bounds for kv filtering
//@param tlt: location to write upper and lower bounds for tilt filtering
//@param el : location to write bitmask of elements
//@param ext: is element bitmask exact or inclusive
//@param sg : location to write upper and lower bounds for space group filtering
void MasterFileFilterPanel::getBounds(std::pair<float, float>& kv, std::pair<float, float>& tlt, ElementMask& el, bool& ext, std::pair<int  , int  >& sg) const {
	if(m_chkKv ->GetValue()) {
		double vMin, vMax;
		m_kvMin ->GetLineText(0).ToDouble(&vMin);
		m_kvMax ->GetLineText(0).ToDouble(&vMax);
		kv .first  = (int)vMin;
		kv .second = (int)vMax;
	} else {
		kv.first = kv.second = NAN;
	}

	if(m_chkTlt->GetValue()) {
		double vMin, vMax;
		m_tltMin->GetLineText(0).ToDouble(&vMin);
		m_tltMax->GetLineText(0).ToDouble(&vMax);
		tlt.first  = (int)vMin;
		tlt.second = (int)vMax;
	} else {
		tlt.first = tlt.second = NAN;
	}

	if(m_chkSg ->GetValue()) {
		long vMin, vMax;
		m_sgMin ->GetLineText(0).ToLong(&vMin);
		m_sgMax ->GetLineText(0).ToLong(&vMax);
		sg .first  = (int)vMin;
		sg .second = (int)vMax;
	} else {
		sg.first = 1; sg.second = 230;
	}

	el = m_prdTbl->getMask();
	ext = m_btnExt->GetValue();
}

void MasterFileFilterPanel::ChoiceChanged( wxCommandEvent& event ) {
	int vMin = 1, vMax = 230;
	switch(event.GetInt()) {
		case  0: vMin =   1, vMax =   2; break;
		case  1: vMin =   3, vMax =  15; break;
		case  2: vMin =  16, vMax =  74; break;
		case  3: vMin =  75, vMax =  88; break;
		case  4: vMin =  89, vMax = 142; break;
		case  5: vMin = 143, vMax = 148; break;
		case  6: vMin = 149, vMax = 167; break;
		case  7: vMin = 168, vMax = 176; break;
		case  8: vMin = 177, vMax = 194; break;
		case  9: vMin = 195, vMax = 206; break;
		case 10: vMin = 207, vMax = 230; break;
	}
	if(vMin > 1 || vMax < 230) {
		m_chkSg ->SetValue(true);
		m_sgMin ->Clear(); m_sgMin ->operator<<(vMin);
		m_sgMax ->Clear(); m_sgMax ->operator<<(vMax);
	}
}

MasterFileFilterPanel::MasterFileFilterPanel( wxWindow* parent, wxWindowID id, const wxPoint& pos, const wxSize& size, long style, const wxString& name ) : wxPanel( parent, id, pos, size, style, name ) {
	//build sizers
	wxBoxSizer     * bSizer  = new wxBoxSizer( wxVERTICAL );
	wxFlexGridSizer* fgSizer = new wxFlexGridSizer( 4, 5, 0, 0 );
	fgSizer->AddGrowableCol( 0 ); fgSizer->AddGrowableCol( 1 ); fgSizer->AddGrowableCol( 2 ); fgSizer->AddGrowableCol( 3 ); fgSizer->AddGrowableCol( 4 );
	fgSizer->AddGrowableRow( 0 ); fgSizer->AddGrowableRow( 1 ); fgSizer->AddGrowableRow( 2 ); fgSizer->AddGrowableRow( 3 );
	fgSizer->SetFlexibleDirection( wxHORIZONTAL );
	fgSizer->SetNonFlexibleGrowMode( wxFLEX_GROWMODE_SPECIFIED );

	//build labels
	wxStaticText* dash1 = new wxStaticText( this, wxID_ANY, wxT("-"), wxDefaultPosition, wxDefaultSize, 0 );
	wxStaticText* dash2 = new wxStaticText( this, wxID_ANY, wxT("-"), wxDefaultPosition, wxDefaultSize, 0 );
	wxStaticText* dash3 = new wxStaticText( this, wxID_ANY, wxT("-"), wxDefaultPosition, wxDefaultSize, 0 );

	//build validators
	valKv .SetRange(0.0f, 50.0f);
	valTlt.SetRange(0.0f, 90.0f);
	valSg .SetRange(1   , 230  );
	valKv .SetPrecision(1);
	valTlt.SetPrecision(1);

	//build elements
	wxString choices[] = { 
		wxT("Triclinic (-1) [1, 2]"),
		wxT("Monoclinic (2/m) [3, 15]"),
		wxT("Orthorhombic (mmm) [16, 74]"),
		wxT("Tetragonal Low (4/m) [75, 88]"),
		wxT("Tetragonal High (4/mmm) [89, 142]"),
		wxT("Triclinic Low (-3) [143, 148]"),
		wxT("Triclinic High (3/m) [149, 167]"),
		wxT("Hexagonal Low (6/m) [168, 176]"),
		wxT("Hexagaonl High (6/mmm) [177, 194]"),
		wxT("Cubic Low (m-3) [195, 206]"),
		wxT("Cubic High (m-3m) [207, 230]")
	};
	m_chkKv  = new wxCheckBox        ( this, wxID_ANY, wxT("kV"         ), wxDefaultPosition, wxDefaultSize, 0               );
	m_kvMin  = new wxTextCtrl        ( this, wxID_ANY, wxT("0.0"        ), wxDefaultPosition, wxDefaultSize, 0 , valKv       );
	m_kvMax  = new wxTextCtrl        ( this, wxID_ANY, wxT("50.0"       ), wxDefaultPosition, wxDefaultSize, 0 , valKv       );
	m_chkTlt = new wxCheckBox        ( this, wxID_ANY, wxT("Tilt"       ), wxDefaultPosition, wxDefaultSize, 0               );
	m_tltMin = new wxTextCtrl        ( this, wxID_ANY, wxT("0.0"        ), wxDefaultPosition, wxDefaultSize, 0 , valTlt      );
	m_tltMax = new wxTextCtrl        ( this, wxID_ANY, wxT("90.0"       ), wxDefaultPosition, wxDefaultSize, 0 , valTlt      );
	m_chkSg  = new wxCheckBox        ( this, wxID_ANY, wxT("Space Group"), wxDefaultPosition, wxDefaultSize, 0               );
	m_sgMin  = new wxTextCtrl        ( this, wxID_ANY, wxT("1"          ), wxDefaultPosition, wxDefaultSize, 0 , valSg       );
	m_sgMax  = new wxTextCtrl        ( this, wxID_ANY, wxT("230"        ), wxDefaultPosition, wxDefaultSize, 0 , valSg       );
	m_laue   = new wxChoice          ( this, wxID_ANY,                     wxDefaultPosition, wxDefaultSize, 11, choices, 0  );
	m_prdTbl = new PeriodicTablePanel( this, wxID_ANY,                     wxDefaultPosition, wxDefaultSize, wxTAB_TRAVERSAL );
	m_laue->SetSelection( -1 );
	m_chkEl  = new wxCheckBox        ( this, wxID_ANY, wxT("Elements"   ), wxDefaultPosition, wxDefaultSize, 0               );
	m_btnExt = new wxToggleButton    ( this, wxID_ANY, wxT("Exact Match"), wxDefaultPosition, wxDefaultSize, 0               );
	m_btnClr = new wxButton          ( this, wxID_ANY, wxT("Clear"      ), wxDefaultPosition, wxDefaultSize, 0               );

	//assemble grid
	fgSizer->Add( m_chkKv , 0, wxALL|wxALIGN_CENTER_VERTICAL                          , 5 );
	fgSizer->Add( m_kvMin , 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_CENTER_HORIZONTAL, 5 );
	fgSizer->Add( dash1   , 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_CENTER_HORIZONTAL, 5 );
	fgSizer->Add( m_kvMax , 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_CENTER_HORIZONTAL, 5 );
	fgSizer->AddStretchSpacer();

	fgSizer->Add( m_chkTlt, 0, wxALL|wxALIGN_CENTER_VERTICAL                          , 5 );
	fgSizer->Add( m_tltMin, 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_CENTER_HORIZONTAL, 5 );
	fgSizer->Add( dash2   , 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_CENTER_HORIZONTAL, 5 );
	fgSizer->Add( m_tltMax, 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_CENTER_HORIZONTAL, 5 );
	fgSizer->AddStretchSpacer();

	fgSizer->Add( m_chkSg , 0, wxALL|wxALIGN_CENTER_VERTICAL                          , 5 );
	fgSizer->Add( m_sgMin , 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_CENTER_HORIZONTAL, 5 );
	fgSizer->Add( dash3   , 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_CENTER_HORIZONTAL, 5 );
	fgSizer->Add( m_sgMax , 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_CENTER_HORIZONTAL, 5 );
	fgSizer->Add( m_laue  , 0, wxALL, 5 );

	fgSizer->Add( m_chkEl , 0, wxALL|wxALIGN_CENTER_VERTICAL                          , 5 );
	fgSizer->Add( m_btnExt, 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_CENTER_HORIZONTAL, 5 );
	fgSizer->AddStretchSpacer();
	fgSizer->Add( m_btnClr, 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_CENTER_HORIZONTAL, 5 );
	fgSizer->AddStretchSpacer();

	//final assembly of parts
	bSizer->Add( fgSizer, 0, wxALIGN_CENTER_HORIZONTAL, 5 );
	// bSizer->Add( m_prdTbl, 1, wxSHAPED|wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_HORIZONTAL );
	bSizer->Add( m_prdTbl, 1, wxEXPAND );
	this->SetMinSize(bSizer->GetMinSize());
	this->SetSizer( bSizer );
	this->Layout();

	m_laue->Connect( wxEVT_COMMAND_CHOICE_SELECTED, wxCommandEventHandler( MasterFileFilterPanel::ChoiceChanged ), NULL, this );
	m_chkKv ->Connect( wxEVT_CHECKBOX              , wxCommandEventHandler( MasterFileFilterPanel::ChkKvChanged     ), NULL, this );
	m_chkTlt->Connect( wxEVT_CHECKBOX              , wxCommandEventHandler( MasterFileFilterPanel::ChkTltChanged    ), NULL, this );
	m_chkSg ->Connect( wxEVT_CHECKBOX              , wxCommandEventHandler( MasterFileFilterPanel::ChkSgChanged     ), NULL, this );
	m_chkEl ->Connect( wxEVT_CHECKBOX              , wxCommandEventHandler( MasterFileFilterPanel::ChkElChanged     ), NULL, this );
	m_kvMin ->Connect( wxEVT_COMMAND_TEXT_UPDATED  , wxCommandEventHandler( MasterFileFilterPanel::TxtKvMinChanged  ), NULL, this );
	m_kvMax ->Connect( wxEVT_COMMAND_TEXT_UPDATED  , wxCommandEventHandler( MasterFileFilterPanel::TxtKvMaxChanged  ), NULL, this );
	m_tltMin->Connect( wxEVT_COMMAND_TEXT_UPDATED  , wxCommandEventHandler( MasterFileFilterPanel::TxtTltMinChanged ), NULL, this );
	m_tltMax->Connect( wxEVT_COMMAND_TEXT_UPDATED  , wxCommandEventHandler( MasterFileFilterPanel::TxtTltMaxChanged ), NULL, this );
	m_sgMin ->Connect( wxEVT_COMMAND_TEXT_UPDATED  , wxCommandEventHandler( MasterFileFilterPanel::TxtSgMinChanged  ), NULL, this );
	m_sgMax ->Connect( wxEVT_COMMAND_TEXT_UPDATED  , wxCommandEventHandler( MasterFileFilterPanel::TxtSgMaxChanged  ), NULL, this );
	m_btnClr->Connect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( MasterFileFilterPanel::ClearPressed     ), NULL, this );
	m_prdTbl->connectToggleEvent(wxCommandEventHandler( MasterFileFilterPanel::ElementsChanged ),  this);

}

MasterFileFilterPanel::~MasterFileFilterPanel() {
	m_laue->Disconnect( wxEVT_COMMAND_CHOICE_SELECTED, wxCommandEventHandler( MasterFileFilterPanel::ChoiceChanged ), NULL, this );
	m_chkKv ->Disconnect( wxEVT_CHECKBOX              , wxCommandEventHandler( MasterFileFilterPanel::ChkKvChanged     ), NULL, this );
	m_chkTlt->Disconnect( wxEVT_CHECKBOX              , wxCommandEventHandler( MasterFileFilterPanel::ChkTltChanged    ), NULL, this );
	m_chkSg ->Disconnect( wxEVT_CHECKBOX              , wxCommandEventHandler( MasterFileFilterPanel::ChkSgChanged     ), NULL, this );
	m_chkEl ->Disconnect( wxEVT_CHECKBOX              , wxCommandEventHandler( MasterFileFilterPanel::ChkElChanged     ), NULL, this );
	m_kvMin ->Disconnect( wxEVT_COMMAND_TEXT_UPDATED  , wxCommandEventHandler( MasterFileFilterPanel::TxtKvMinChanged  ), NULL, this );
	m_kvMax ->Disconnect( wxEVT_COMMAND_TEXT_UPDATED  , wxCommandEventHandler( MasterFileFilterPanel::TxtKvMaxChanged  ), NULL, this );
	m_tltMin->Disconnect( wxEVT_COMMAND_TEXT_UPDATED  , wxCommandEventHandler( MasterFileFilterPanel::TxtTltMinChanged ), NULL, this );
	m_tltMax->Disconnect( wxEVT_COMMAND_TEXT_UPDATED  , wxCommandEventHandler( MasterFileFilterPanel::TxtTltMaxChanged ), NULL, this );
	m_sgMin ->Disconnect( wxEVT_COMMAND_TEXT_UPDATED  , wxCommandEventHandler( MasterFileFilterPanel::TxtSgMinChanged  ), NULL, this );
	m_sgMax ->Disconnect( wxEVT_COMMAND_TEXT_UPDATED  , wxCommandEventHandler( MasterFileFilterPanel::TxtSgMaxChanged  ), NULL, this );
	m_btnClr->Disconnect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( MasterFileFilterPanel::ClearPressed     ), NULL, this );
	m_prdTbl->disconnectToggleEvent(wxCommandEventHandler( MasterFileFilterPanel::ElementsChanged ),  this);
}

#endif//_MP_FILT_PAN_H_
