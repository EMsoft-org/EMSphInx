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

#ifndef _PAT_CEN_PAN_H_
#define _PAT_CEN_PAN_H_

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

#include "wx/ValidityPanel.h"

///////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
/// Class PatternCenterPanel
///////////////////////////////////////////////////////////////////////////////
class PatternCenterPanel : public ValidityPanel
{
	private:
		size_t m_lstBin = 1;//last binning value
		size_t m_binW, m_binH;//binned detector size (must be set from external knowledge)

		wxIntegerValidator      <int  > valBin;//bining
		wxFloatingPointValidator<float> valDlt;//pixel size
		wxFloatingPointValidator<float> valWdt;//detector width
		wxFloatingPointValidator<float> valPc ;//pcx/pcy
		wxFloatingPointValidator<float> valL  ;//L
		wxFloatingPointValidator<float> xyStar;//x*/y*
		wxFloatingPointValidator<float> zStar ;//z*
		wxFloatingPointValidator<float> valTlt;//tilt

	protected:
		wxTextCtrl* m_txtBin   ;
		wxTextCtrl* m_txtUnDetW;
		wxTextCtrl* m_txtPxSz  ;
		wxTextCtrl* m_txtDetW  ;
		wxChoice  * m_chcVendor;
		wxTextCtrl* m_txtPcx   ;
		wxTextCtrl* m_txtXstar ;
		wxTextCtrl* m_txtPcy   ;
		wxTextCtrl* m_txtYstar ;
		wxTextCtrl* m_txtL     ;
		wxTextCtrl* m_txtZstar ;
		wxTextCtrl* m_txtTlt   ;
		wxButton  * m_btnFit   ;

		//@brief    : enable/disable the pattern center part of the panel
		//@param enb: enable state
		void EnablePatCen(bool enb);

		//@brief: compute the emsoft pattern center from x*, y*, z* if possible
		void computeEMsoft();

		//@brief: compute the vendor pattern center from emsoft pattern center if possible
		void computeVendor();

		// Virtual event handlers, overide them in your derived class
		void BinningChanged( wxCommandEvent& event ) { binChged  ();    testValid();}
		void PixSzChanged  ( wxCommandEvent& event ) { pixSzChged();    testValid();}
		void DetWidChanged ( wxCommandEvent& event ) { detWChged ();    testValid();}
		void VendorChanged ( wxCommandEvent& event ) { computeVendor(); testValid();}
		void PcxChanged    ( wxCommandEvent& event ) { computeVendor(); testValid();}
		void PcyChanged    ( wxCommandEvent& event ) { computeVendor(); testValid();}
		void LChanged      ( wxCommandEvent& event ) { computeVendor(); testValid();}
		void XstarChanged  ( wxCommandEvent& event ) { computeEMsoft(); testValid();}
		void YstarChanged  ( wxCommandEvent& event ) { computeEMsoft(); testValid();}
		void ZstarChanged  ( wxCommandEvent& event ) { computeEMsoft(); testValid();}
		void TltChanged    ( wxCommandEvent& event ) {                  testValid();}

		void DoFit         ( wxCommandEvent& event ) { event.Skip(); }

		void binChged ();
		void pixSzChged();//update detector size from pixel size
		void detWChged ();//update pixel size from pixel size

		//check if text fields are empty
		bool hasBin () const {return !m_txtBin   ->GetLineText(0).IsEmpty();}
		bool hasUbDW() const {return !m_txtUnDetW->GetLineText(0).IsEmpty();}
		bool hasPxSz() const {return !m_txtPxSz  ->GetLineText(0).IsEmpty();}
		bool hasDetW() const {return !m_txtDetW  ->GetLineText(0).IsEmpty();}
		bool hasPcx () const {return !m_txtPcx   ->GetLineText(0).IsEmpty();}
		bool hasPcy () const {return !m_txtPcy   ->GetLineText(0).IsEmpty();}
		bool hasL   () const {return !m_txtL     ->GetLineText(0).IsEmpty();}
		bool hasXst () const {return !m_txtXstar ->GetLineText(0).IsEmpty();}
		bool hasYst () const {return !m_txtYstar ->GetLineText(0).IsEmpty();}
		bool hasZst () const {return !m_txtZstar ->GetLineText(0).IsEmpty();}
		bool hasTlt () const {return !m_txtTlt   ->GetLineText(0).IsEmpty();}

		//get numeric values of fields
		long   getBin () const {long   v; m_txtBin   ->GetLineText(0).ToLong  (&v); return v;}
		long   getUbDW() const {long   v; m_txtUnDetW->GetLineText(0).ToLong  (&v); return v;}
		double getPxSz() const {double v; m_txtPxSz  ->GetLineText(0).ToDouble(&v); return v;}
		double getDetW() const {double v; m_txtDetW  ->GetLineText(0).ToDouble(&v); return v;}
		double getPcx () const {double v; m_txtPcx   ->GetLineText(0).ToDouble(&v); return v;}
		double getPcy () const {double v; m_txtPcy   ->GetLineText(0).ToDouble(&v); return v;}
		double getL   () const {double v; m_txtL     ->GetLineText(0).ToDouble(&v); return v;}
		double getXst () const {double v; m_txtXstar ->GetLineText(0).ToDouble(&v); return v;}
		double getYst () const {double v; m_txtYstar ->GetLineText(0).ToDouble(&v); return v;}
		double getZst () const {double v; m_txtZstar ->GetLineText(0).ToDouble(&v); return v;}
		double getTlt () const {double v; m_txtTlt   ->GetLineText(0).ToDouble(&v); return v;}

		//set numeric values of fields
		void setBin (long   v) const {m_txtBin   ->Clear(); m_txtBin   ->operator<<(v);}
		void setUbDW(long   v) const {m_txtUnDetW->Clear(); m_txtUnDetW->operator<<(v);}
		void setPxSz(double v) const {m_txtPxSz  ->Clear(); m_txtPxSz  ->operator<<(v);}
		void setDetW(double v) const {m_txtDetW  ->Clear(); m_txtDetW  ->operator<<(v);}
		void setPcx (double v) const {m_txtPcx   ->Clear(); m_txtPcx   ->operator<<(v);}
		void setPcy (double v) const {m_txtPcy   ->Clear(); m_txtPcy   ->operator<<(v);}
		void setL   (double v) const {m_txtL     ->Clear(); m_txtL     ->operator<<(v);}
		void setXst (double v) const {m_txtXstar ->Clear(); m_txtXstar ->operator<<(v);}
		void setYst (double v) const {m_txtYstar ->Clear(); m_txtYstar ->operator<<(v);}
		void setZst (double v) const {m_txtZstar ->Clear(); m_txtZstar ->operator<<(v);}
		void setTlt (double v) const {m_txtTlt   ->Clear(); m_txtTlt   ->operator<<(v);}

	public:
		//@brief: sanity check the current state
		//@return: true if the values parsed from the panel are reasonable, false otherwise
		//@note  : checks for has a file, has detector sizes, and has an AHE value
		std::string validMsg() const {
			if(!hasBin()) return "binning empty";
			if(!hasPxSz()) return "binned pixel size empty";
			if(!hasDetW()) return "detector width empty";
			if(!hasPcx()) return "pcx empty";
			if(!hasPcy()) return "pcy empty";
			if(!hasL()) return "L empty";
			if(!hasTlt()) return "detector tilt empty";
			return "";
		}

		PatternCenterPanel( wxWindow* parent, wxWindowID id = wxID_ANY, const wxPoint& pos = wxDefaultPosition, const wxSize& size = wxSize( 400,400 ), long style = wxTAB_TRAVERSAL, const wxString& name = wxEmptyString );
		~PatternCenterPanel();

		//@brief  : set the binned detector width in pixels (hidden field)
		//@param w: binned detector width in pixels
		//@param w: binned detector height in pixels
		void setBinnedPix(const size_t w, const size_t h) {m_binW = w; m_binH = h; setBin(1);}//should also call same function as BinningChanged

		//@brief : get the pixel size
		//@return: pixel size in microns
		double getDelta() const {return getPxSz() * getBin();}

		//@brief    : get the pixel size
		//@param dlt: pixel size in microns
		void setDelta(double dlt) {setPxSz(dlt); pixSzChged();}

		//@brief    : update the pattern center
		//@param x  : xstar (or pcx)
		//@param y  : ystar (or pcy)
		//@param z  : zstar (or L  )
		//@param ven: vendor (must be "EMsoft", "EDAX", "Oxford", or "Bruker")
		//@note     : if ven == EMsoft the current x/y/z* will be updated, otherwise EMsoft will be updated
		void setPatternCenter(const double x, const double y, const double z, std::string ven);

		//@brief: clear pattern center info
		void clear();

		//@brief  : get EMsoft pattern center
		//@param x: location to write pcx
		//@param y: location to write pcy
		//@param z: location to write L
		void getPatternCenter(double& x, double& y, double& z) {x = getPcx(); y = getPcy(); z = getL();}

		//@brief  : set the detector tilt
		//@param t: detector tilt in degrees
		void setDetTlt(const double t) {setTlt(t);}

		//@brief : get the detector tilt
		//@return: detector tilt in degrees
		double getDetTlt() const {return getTlt();}

};

///////////////////////////////////////////////////////////////////////////
// C++ code generated with wxFormBuilder (version Oct 26 2018)
// http://www.wxformbuilder.org/
//
// PLEASE DO *NOT* EDIT THIS FILE!
///////////////////////////////////////////////////////////////////////////

void PatternCenterPanel::clear() {
	EnablePatCen(false);
	m_txtBin   ->Clear();
	m_txtUnDetW->Clear();
	m_txtPxSz  ->Clear();
	m_txtDetW  ->Clear();
	m_txtPcx   ->Clear();
	m_txtPcy   ->Clear();
	m_txtL     ->Clear();
	m_txtXstar ->Clear();
	m_txtYstar ->Clear();
	m_txtZstar ->Clear();
	m_txtTlt   ->Clear();
	m_chcVendor->SetSelection(0);
}

//@brief    : enable/disable the pattern center part of the panel
//@param enb: enable state
void PatternCenterPanel::EnablePatCen(bool enb) {
	m_chcVendor->Enable(enb);
	m_txtPcx   ->Enable(enb);
	m_txtXstar ->Enable(enb);
	m_txtPcy   ->Enable(enb);
	m_txtYstar ->Enable(enb);
	m_txtL     ->Enable(enb);
	m_txtZstar ->Enable(enb);
	m_txtTlt   ->Enable(enb);
	// m_btnFit   ->Enable(enb);
}

//@brief: compute the emsoft pattern center from x*, y*, z* if possible
void PatternCenterPanel::computeEMsoft() {
	if(hasXst() && hasYst() && hasZst() && hasPxSz() && hasBin()) {
		double pctr[3] = {
			getXst(),
			getYst(),
			getZst(),
		};
		double delta = getDelta();
		switch(m_chcVendor->GetSelection()) {
			case 0://bruker
				pctr[0] = pctr[0] * m_binW - 0.5 * m_binW;
				pctr[1] = 0.5 * m_binH - pctr[1] * m_binH;
				pctr[2] = pctr[2] * m_binH * delta;
			break;

			case 1://edax
				pctr[0] = pctr[0] * m_binW - 0.5 * m_binW;
				pctr[1] = pctr[1] * m_binW - 0.5 * m_binH;
				pctr[2] = pctr[2] * m_binW * delta;
			break;

			case 2://oxford
				pctr[0] = pctr[0] * m_binW - 0.5 * m_binW;
				pctr[1] = pctr[1] * m_binH - 0.5 * m_binH;
				pctr[2] = pctr[2] * m_binW * delta;
			break;

			default: return;//nothing selected
		}
		m_txtPcx->ChangeValue(wxString::Format(wxT("%0.3f"), pctr[0]));
		m_txtPcy->ChangeValue(wxString::Format(wxT("%0.3f"), pctr[1]));
		m_txtL  ->ChangeValue(wxString::Format(wxT("%0.2f"), pctr[2]));
	}
}

//@brief: compute the vendor pattern center from emsoft pattern center if possible
void PatternCenterPanel::computeVendor() {
	if(hasPcx() && hasPcy() && hasL() && hasPxSz() && hasBin()) {
		double pctr[3] = {
			getPcx(),
			getPcy(),
			getL  (),
		};
		double delta = getDelta();
		switch(m_chcVendor->GetSelection()) {
			case 0://bruker
				pctr[0] = (pctr[0] + 0.5 * m_binW) / m_binW;
				pctr[1] = (0.5 * m_binH - pctr[1]) / m_binH;
				pctr[2] = pctr[2] / (delta * m_binH);
			break;

			case 1://edax
				pctr[0] = (pctr[0] + 0.5 * m_binW) / m_binW;
				pctr[1] = (pctr[1] + 0.5 * m_binH) / m_binW;
				pctr[2] = pctr[2] / (delta * m_binW);
			break;

			case 2://oxford
				pctr[0] = (pctr[0] + 0.5 * m_binW) / m_binW;
				pctr[1] = (pctr[1] + 0.5 * m_binH) / m_binH;
				pctr[2] = pctr[2] / (delta * m_binW);
			break;

			default: return;//nothing selected
		}
		m_txtXstar->ChangeValue(wxString::Format(wxT("%0.6f"), pctr[0]));
		m_txtYstar->ChangeValue(wxString::Format(wxT("%0.6f"), pctr[1]));
		m_txtZstar->ChangeValue(wxString::Format(wxT("%0.6f"), pctr[2]));
	}
}

void PatternCenterPanel::binChged() {
	//compute unbinned size in pixels
	if(hasBin()) {
		long bin = getBin();
		long unbinPix = bin * m_binW;
		m_txtUnDetW->Clear();
		m_txtUnDetW->operator<<(unbinPix);
		if(hasPxSz()) detWChged();//if we had a pixel size, update it to keep the detector width fixed
		m_lstBin = bin;
	}
}

//update detector size from pixel size
void PatternCenterPanel::pixSzChged() {
	if(hasPxSz()) {
		EnablePatCen(true);
		m_txtDetW->ChangeValue(wxString::Format(wxT("%0.4f"), getPxSz() * (m_binW * getBin()) / 1000));//change value so we dont' issue an event
		computeEMsoft();
	} else {
		EnablePatCen(false);
	}
}

//update pixel size from pixel size
void PatternCenterPanel::detWChged () {
	if(hasDetW()) {
		EnablePatCen(true);
		m_txtPxSz->ChangeValue(wxString::Format(wxT("%0.2f"), getDetW() * 1000.0 / (m_binW * getBin())));
		computeEMsoft();
	} else {
		EnablePatCen(false);
	}
}

//@brief    : update the pattern center
//@param x  : xstar (or pcx)
//@param y  : ystar (or pcy)
//@param z  : zstar (or L  )
//@param ven: vendor (must be "EMsoft", "EDAX", "Oxford", or "Bruker")
//@note     : if ven == EMsoft the current x/y/z* will be updated, otherwise EMsoft will be updated
void PatternCenterPanel::setPatternCenter(const double x, const double y, const double z, std::string ven) {
	std::transform(ven.begin(), ven.end(), ven.begin(), [](unsigned char c){ return std::tolower(c); });
	bool vendor = true;
	clear();
	if("bruker" == ven) {
		m_chcVendor->SetSelection(0);
	} else if("edax" == ven) {
		m_chcVendor->SetSelection(1);
	} else if("oxford" == ven) {
		m_chcVendor->SetSelection(2);
	} else if("emsoft" == ven) {
		m_txtPcx->ChangeValue(wxString::Format(wxT("%0.3f"), x));
		m_txtPcy->ChangeValue(wxString::Format(wxT("%0.3f"), y));
		m_txtL  ->ChangeValue(wxString::Format(wxT("%0.2f"), z));
		vendor = false;
		computeVendor();
	} else {
		throw std::runtime_error("unknown vendor " + ven);
	}
	if(vendor) {
		m_txtXstar->ChangeValue(wxString::Format(wxT("%0.6f"), x));
		m_txtYstar->ChangeValue(wxString::Format(wxT("%0.6f"), y));
		m_txtZstar->ChangeValue(wxString::Format(wxT("%0.6f"), z));
		computeEMsoft();
	}
}

///////////////////////////////////////////////////////////////////////////

PatternCenterPanel::PatternCenterPanel( wxWindow* parent, wxWindowID id, const wxPoint& pos, const wxSize& size, long style, const wxString& name ) : ValidityPanel( parent, id, pos, size, style, name ) {
	//split panel into 2 vertical boxes
	wxBoxSizer      * bPatCent  = new wxBoxSizer( wxVERTICAL );
	wxStaticBoxSizer* sbPixSize = new wxStaticBoxSizer( new wxStaticBox( this, wxID_ANY, wxT("Pixel Size"    ) ), wxVERTICAL );
	wxStaticBoxSizer* sbPatCen  = new wxStaticBoxSizer( new wxStaticBox( this, wxID_ANY, wxT("Pattern Center") ), wxVERTICAL );
	bPatCent->Add( sbPixSize, 3, wxEXPAND, 5 );
	bPatCent->Add( sbPatCen , 4, wxEXPAND, 5 );

	//top box is 3x4 grid
	wxFlexGridSizer* fgPixSz = new wxFlexGridSizer( 4, 3, 0, 0 );
	fgPixSz->AddGrowableCol( 0 ); fgPixSz->AddGrowableCol( 1 ); fgPixSz->AddGrowableCol( 2 );
	fgPixSz->AddGrowableRow( 0 ); fgPixSz->AddGrowableRow( 1 ); fgPixSz->AddGrowableRow( 2 ); fgPixSz->AddGrowableRow( 3 );
	fgPixSz->SetFlexibleDirection( wxHORIZONTAL );
	fgPixSz->SetNonFlexibleGrowMode( wxFLEX_GROWMODE_SPECIFIED );
	sbPixSize->Add( fgPixSz, 1, wxEXPAND, 5 );

	//bottom box is 8x5 grid
	wxFlexGridSizer* fgPtCen = new wxFlexGridSizer( 5, 8, 0, 0 );
	fgPtCen->AddGrowableCol( 0 ); fgPtCen->AddGrowableCol( 4 ); fgPtCen->AddGrowableCol( 7 );
	fgPtCen->AddGrowableRow( 0 ); fgPtCen->AddGrowableRow( 1 ); fgPtCen->AddGrowableRow( 2 ); fgPtCen->AddGrowableRow( 3 ); fgPtCen->AddGrowableRow( 4 );
	fgPtCen->SetFlexibleDirection( wxHORIZONTAL );
	fgPtCen->SetNonFlexibleGrowMode(  wxFLEX_GROWMODE_SPECIFIED );
	sbPatCen->Add( fgPtCen, 1, wxEXPAND, 5 );

	//build all labels
	wxStaticText* staTxtBin   = new wxStaticText( sbPixSize->GetStaticBox(), wxID_ANY, wxT("Binning"                ), wxDefaultPosition, wxDefaultSize, 0            ); staTxtBin  ->Wrap(-1);
	wxStaticText* staTxtX     = new wxStaticText( sbPixSize->GetStaticBox(), wxID_ANY, wxT("x"                      ), wxDefaultPosition, wxDefaultSize, 0            ); staTxtX    ->Wrap(-1);
	wxStaticText* staTxtUnbW  = new wxStaticText( sbPixSize->GetStaticBox(), wxID_ANY, wxT("Unbinned Detector Width"), wxDefaultPosition, wxDefaultSize, 0            ); staTxtUnbW ->Wrap(-1);
	wxStaticText* staTxtPix   = new wxStaticText( sbPixSize->GetStaticBox(), wxID_ANY, wxT("pix"                    ), wxDefaultPosition, wxDefaultSize, wxALIGN_LEFT ); staTxtPix  ->Wrap(-1);
	wxStaticText* staTxtPxSz  = new wxStaticText( sbPixSize->GetStaticBox(), wxID_ANY, wxT("Binned Pixel Size"      ), wxDefaultPosition, wxDefaultSize, 0            ); staTxtPxSz ->Wrap(-1);
	wxStaticText* staTxtUm    = new wxStaticText( sbPixSize->GetStaticBox(), wxID_ANY, wxT("um"                     ), wxDefaultPosition, wxDefaultSize, wxALIGN_LEFT ); staTxtUm   ->Wrap(-1);
	wxStaticText* staTxtDetW  = new wxStaticText( sbPixSize->GetStaticBox(), wxID_ANY, wxT("Detector Width"         ), wxDefaultPosition, wxDefaultSize, 0            ); staTxtDetW ->Wrap(-1);
	wxStaticText* staTxtMm    = new wxStaticText( sbPixSize->GetStaticBox(), wxID_ANY, wxT("mm"                     ), wxDefaultPosition, wxDefaultSize, wxALIGN_LEFT ); staTxtMm   ->Wrap(-1);
	wxStaticText* staTxtEmSft = new wxStaticText( sbPatCen ->GetStaticBox(), wxID_ANY, wxT("EMsoft"                 ), wxDefaultPosition, wxDefaultSize, 0            ); staTxtEmSft->Wrap(-1);
	wxStaticText* staTxtPcx   = new wxStaticText( sbPatCen ->GetStaticBox(), wxID_ANY, wxT("pcx"                    ), wxDefaultPosition, wxDefaultSize, 0            ); staTxtPcx  ->Wrap(-1);
	wxStaticText* staTxtPcxU  = new wxStaticText( sbPatCen ->GetStaticBox(), wxID_ANY, wxT("pix"                    ), wxDefaultPosition, wxDefaultSize, 0            ); staTxtPcxU ->Wrap(-1);
	wxStaticText* staTxtXstar = new wxStaticText( sbPatCen ->GetStaticBox(), wxID_ANY, wxT("x*"                     ), wxDefaultPosition, wxDefaultSize, 0            ); staTxtXstar->Wrap(-1);
	wxStaticText* staTxtPcy   = new wxStaticText( sbPatCen ->GetStaticBox(), wxID_ANY, wxT("pcy"                    ), wxDefaultPosition, wxDefaultSize, 0            ); staTxtPcy  ->Wrap(-1);
	wxStaticText* staTxtPcyU  = new wxStaticText( sbPatCen ->GetStaticBox(), wxID_ANY, wxT("pix"                    ), wxDefaultPosition, wxDefaultSize, 0            ); staTxtPcyU ->Wrap(-1);
	wxStaticText* staTxtYstar = new wxStaticText( sbPatCen ->GetStaticBox(), wxID_ANY, wxT("y*"                     ), wxDefaultPosition, wxDefaultSize, 0            ); staTxtYstar->Wrap(-1);
	wxStaticText* staTxtL     = new wxStaticText( sbPatCen ->GetStaticBox(), wxID_ANY, wxT("L"                      ), wxDefaultPosition, wxDefaultSize, 0            ); staTxtL    ->Wrap(-1);
	wxStaticText* staTxtLU    = new wxStaticText( sbPatCen ->GetStaticBox(), wxID_ANY, wxT("um"                     ), wxDefaultPosition, wxDefaultSize, 0            ); staTxtLU   ->Wrap(-1);
	wxStaticText* staTxtZstar = new wxStaticText( sbPatCen ->GetStaticBox(), wxID_ANY, wxT("z*"                     ), wxDefaultPosition, wxDefaultSize, 0            ); staTxtZstar->Wrap(-1);
	wxStaticText* staTxtTlt   = new wxStaticText( sbPatCen ->GetStaticBox(), wxID_ANY, wxT("Detector Tilt"          ), wxDefaultPosition, wxDefaultSize, 0            ); staTxtTlt  ->Wrap(-1);
	wxStaticText* staTxtDeg   = new wxStaticText( sbPatCen ->GetStaticBox(), wxID_ANY, wxT("deg"                    ), wxDefaultPosition, wxDefaultSize, wxALIGN_LEFT ); staTxtDeg  ->Wrap(-1);

	//set text validator ranges
	valBin.SetRange( 1    , 16    );
	valDlt.SetRange( 1    , 1200  );//this is [5,75] um * 16x binning (actualy [1,75] to make text entry easier)
	valWdt.SetRange( 1    , 100   );//[1,100] mm
	valPc .SetRange(-1024 , 1024  );//this is x/y* from -1,1 for a 1k detector
	valL  .SetRange( 1    , 50000 );//[1,50] mm (actualy 1um to 50mm to make text entry easier)
	xyStar.SetRange(-1    , 1     );//should be ~0.5 and ~0.7 for pcx and pcy respectively
	zStar .SetRange( 0.01f, 2     );//minimum is 1 mm for 100 mm detector
	valTlt.SetRange(-90   , 90    );

	//specify validator precisions, max 6 for float
	valDlt.SetPrecision(2);
	valWdt.SetPrecision(4);
	valPc .SetPrecision(3);
	valL  .SetPrecision(2);
	xyStar.SetPrecision(6);
	zStar .SetPrecision(6);
	valTlt.SetPrecision(1);

	//build elements for pixel size box
	wxString m_chcVendorChoices[] = { wxT("Bruker"), wxT("EDAX"), wxT("Oxford") };
	int m_chcVendorNChoices = sizeof( m_chcVendorChoices ) / sizeof( wxString );
	m_txtBin    = new wxTextCtrl( sbPixSize->GetStaticBox(), wxID_ANY, wxT(""      ), wxDefaultPosition, wxDefaultSize, wxTE_CENTER              , valBin );
	m_txtUnDetW = new wxTextCtrl( sbPixSize->GetStaticBox(), wxID_ANY, wxT(""      ), wxDefaultPosition, wxDefaultSize, wxTE_CENTER|wxTE_READONLY         );
	m_txtPxSz   = new wxTextCtrl( sbPixSize->GetStaticBox(), wxID_ANY, wxT(""      ), wxDefaultPosition, wxDefaultSize, wxTE_CENTER              , valDlt );
	m_txtDetW   = new wxTextCtrl( sbPixSize->GetStaticBox(), wxID_ANY, wxT(""      ), wxDefaultPosition, wxDefaultSize, wxTE_CENTER              , valWdt );
	
	//build elements for pattern center box
	m_txtPcx    = new wxTextCtrl( sbPatCen ->GetStaticBox(), wxID_ANY, wxT(""      ), wxDefaultPosition, wxDefaultSize, wxTE_CENTER              , valPc  );
	m_txtXstar  = new wxTextCtrl( sbPatCen ->GetStaticBox(), wxID_ANY, wxT(""      ), wxDefaultPosition, wxDefaultSize, wxTE_CENTER              , xyStar );
	m_txtPcy    = new wxTextCtrl( sbPatCen ->GetStaticBox(), wxID_ANY, wxT(""      ), wxDefaultPosition, wxDefaultSize, wxTE_CENTER              , valPc  );
	m_txtYstar  = new wxTextCtrl( sbPatCen ->GetStaticBox(), wxID_ANY, wxT(""      ), wxDefaultPosition, wxDefaultSize, wxTE_CENTER              , xyStar );
	m_txtL      = new wxTextCtrl( sbPatCen ->GetStaticBox(), wxID_ANY, wxT(""      ), wxDefaultPosition, wxDefaultSize, wxTE_CENTER              , valL   );
	m_txtZstar  = new wxTextCtrl( sbPatCen ->GetStaticBox(), wxID_ANY, wxT(""      ), wxDefaultPosition, wxDefaultSize, wxTE_CENTER              , zStar  );
	m_txtTlt    = new wxTextCtrl( sbPatCen ->GetStaticBox(), wxID_ANY, wxT(""      ), wxDefaultPosition, wxDefaultSize, wxTE_CENTER              , valTlt );
	m_chcVendor = new wxChoice  ( sbPatCen ->GetStaticBox(), wxID_ANY,                wxDefaultPosition, wxDefaultSize, m_chcVendorNChoices, m_chcVendorChoices, 0 );
	m_btnFit    = new wxButton  ( sbPatCen ->GetStaticBox(), wxID_ANY, wxT("Fit..."), wxDefaultPosition, wxDefaultSize, 0 );
	m_chcVendor->SetSelection( -1 );
	m_btnFit->Enable(false);//need to write fit routine

	//assemble pixel size grid
	fgPixSz->Add( staTxtBin  , 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_CENTER_HORIZONTAL, 5 );
	fgPixSz->Add( m_txtBin   , 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT            , 5 );
	fgPixSz->Add( staTxtX    , 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_LEFT             , 5 );

	fgPixSz->Add( staTxtUnbW , 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_CENTER_HORIZONTAL, 5 );
	fgPixSz->Add( m_txtUnDetW, 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT            , 5 );
	fgPixSz->Add( staTxtPix  , 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_LEFT             , 5 );

	fgPixSz->Add( staTxtPxSz , 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_CENTER_HORIZONTAL, 5 );
	fgPixSz->Add( m_txtPxSz  , 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT            , 5 );
	fgPixSz->Add( staTxtUm   , 0, wxALL|wxALIGN_CENTER_VERTICAL                          , 5 );

	fgPixSz->Add( staTxtDetW , 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_CENTER_HORIZONTAL, 5 );
	fgPixSz->Add( m_txtDetW  , 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT            , 5 );
	fgPixSz->Add( staTxtMm   , 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_LEFT             , 5 );

	//assemble pattern center grid
	fgPtCen->Add( 0          , 0, 1    , wxEXPAND                                        , 5 );
	fgPtCen->Add( 0          , 0, 1    , wxEXPAND                                        , 5 );
	fgPtCen->Add( staTxtEmSft, 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_CENTER_HORIZONTAL, 5 );
	fgPtCen->Add( 0          , 0, 1    , wxEXPAND                                        , 5 );
	fgPtCen->Add( 0          , 0, 1    , wxEXPAND                                        , 5 );
	fgPtCen->Add( 0          , 0, 1    , wxEXPAND                                        , 5 );
	fgPtCen->Add( m_chcVendor, 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_CENTER_HORIZONTAL, 5 );
	fgPtCen->Add( 0          , 0, 1    , wxEXPAND                                        , 5 );

	fgPtCen->Add( 0          , 0, 1    , wxEXPAND                                        , 5 );
	fgPtCen->Add( staTxtPcx  , 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_CENTER_HORIZONTAL, 5 );
	fgPtCen->Add( m_txtPcx   , 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT            , 5 );
	fgPtCen->Add( staTxtPcxU , 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_LEFT             , 5 );
	fgPtCen->Add( 0          , 0, 1    , wxEXPAND                                        , 5 );
	fgPtCen->Add( staTxtXstar, 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_CENTER_HORIZONTAL, 5 );
	fgPtCen->Add( m_txtXstar , 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT            , 5 );
	fgPtCen->Add( 0          , 0, 1    , wxEXPAND                                        , 5 );

	fgPtCen->Add( 0          , 0, 1    , wxEXPAND                                        , 5 );
	fgPtCen->Add( staTxtPcy  , 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_CENTER_HORIZONTAL, 5 );
	fgPtCen->Add( m_txtPcy   , 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT            , 5 );
	fgPtCen->Add( staTxtPcyU , 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_LEFT             , 5 );
	fgPtCen->Add( 0          , 0, 1    , wxEXPAND                                        , 5 );
	fgPtCen->Add( staTxtYstar, 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_CENTER_HORIZONTAL, 5 );
	fgPtCen->Add( m_txtYstar , 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT            , 5 );
	fgPtCen->Add( 0          , 0, 1    , wxEXPAND                                        , 5 );

	fgPtCen->Add( 0          , 0, 1    , wxEXPAND                                        , 5 );
	fgPtCen->Add( staTxtL    , 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_CENTER_HORIZONTAL, 5 );
	fgPtCen->Add( m_txtL     , 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT            , 5 );
	fgPtCen->Add( staTxtLU   , 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_LEFT             , 5 );
	fgPtCen->Add( 0          , 0, 1    , wxEXPAND                                        , 5 );
	fgPtCen->Add( staTxtZstar, 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_CENTER_HORIZONTAL, 5 );
	fgPtCen->Add( m_txtZstar , 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT            , 5 );
	fgPtCen->Add( 0          , 0, 1    , wxEXPAND                                        , 5 );

	fgPtCen->Add( 0          , 0, 1    , wxEXPAND                                        , 5 );
	fgPtCen->Add( staTxtTlt  , 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_CENTER_HORIZONTAL, 5 );
	fgPtCen->Add( m_txtTlt   , 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT            , 5 );
	fgPtCen->Add( staTxtDeg  , 0, wxALL|wxALIGN_CENTER_VERTICAL                          , 5 );
	fgPtCen->Add( 0          , 0, 1    , wxEXPAND                                        , 5 );
	fgPtCen->Add( 0          , 0, 1    , wxEXPAND                                        , 5 );
	fgPtCen->Add( m_btnFit   , 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT            , 5 );
	fgPtCen->Add( 0          , 0, 1    , wxEXPAND                                        , 5 );

	this->SetSizer( bPatCent );
	this->Layout();
	EnablePatCen(false);

	// Connect Events
	m_txtBin   ->Connect( wxEVT_COMMAND_TEXT_UPDATED   , wxCommandEventHandler( PatternCenterPanel::BinningChanged ), NULL, this );
	m_txtPxSz  ->Connect( wxEVT_COMMAND_TEXT_UPDATED   , wxCommandEventHandler( PatternCenterPanel::PixSzChanged   ), NULL, this );
	m_txtDetW  ->Connect( wxEVT_COMMAND_TEXT_UPDATED   , wxCommandEventHandler( PatternCenterPanel::DetWidChanged  ), NULL, this );
	m_chcVendor->Connect( wxEVT_COMMAND_CHOICE_SELECTED, wxCommandEventHandler( PatternCenterPanel::VendorChanged  ), NULL, this );
	m_txtPcx   ->Connect( wxEVT_COMMAND_TEXT_UPDATED   , wxCommandEventHandler( PatternCenterPanel::PcxChanged     ), NULL, this );
	m_txtXstar ->Connect( wxEVT_COMMAND_TEXT_UPDATED   , wxCommandEventHandler( PatternCenterPanel::XstarChanged   ), NULL, this );
	m_txtPcy   ->Connect( wxEVT_COMMAND_TEXT_UPDATED   , wxCommandEventHandler( PatternCenterPanel::PcyChanged     ), NULL, this );
	m_txtYstar ->Connect( wxEVT_COMMAND_TEXT_UPDATED   , wxCommandEventHandler( PatternCenterPanel::YstarChanged   ), NULL, this );
	m_txtL     ->Connect( wxEVT_COMMAND_TEXT_UPDATED   , wxCommandEventHandler( PatternCenterPanel::LChanged       ), NULL, this );
	m_txtZstar ->Connect( wxEVT_COMMAND_TEXT_UPDATED   , wxCommandEventHandler( PatternCenterPanel::ZstarChanged   ), NULL, this );
	m_txtTlt   ->Connect( wxEVT_COMMAND_TEXT_UPDATED   , wxCommandEventHandler( PatternCenterPanel::TltChanged     ), NULL, this );
	m_btnFit   ->Connect( wxEVT_COMMAND_BUTTON_CLICKED , wxCommandEventHandler( PatternCenterPanel::DoFit          ), NULL, this );
}

PatternCenterPanel::~PatternCenterPanel() {
	// Disconnect Events
	m_txtBin   ->Disconnect( wxEVT_COMMAND_TEXT_UPDATED   , wxCommandEventHandler( PatternCenterPanel::BinningChanged ), NULL, this );
	m_txtPxSz  ->Disconnect( wxEVT_COMMAND_TEXT_UPDATED   , wxCommandEventHandler( PatternCenterPanel::PixSzChanged   ), NULL, this );
	m_txtDetW  ->Disconnect( wxEVT_COMMAND_TEXT_UPDATED   , wxCommandEventHandler( PatternCenterPanel::DetWidChanged  ), NULL, this );
	m_chcVendor->Disconnect( wxEVT_COMMAND_CHOICE_SELECTED, wxCommandEventHandler( PatternCenterPanel::VendorChanged  ), NULL, this );
	m_txtPcx   ->Disconnect( wxEVT_COMMAND_TEXT_UPDATED   , wxCommandEventHandler( PatternCenterPanel::PcxChanged     ), NULL, this );
	m_txtXstar ->Disconnect( wxEVT_COMMAND_TEXT_UPDATED   , wxCommandEventHandler( PatternCenterPanel::XstarChanged   ), NULL, this );
	m_txtPcy   ->Disconnect( wxEVT_COMMAND_TEXT_UPDATED   , wxCommandEventHandler( PatternCenterPanel::PcyChanged     ), NULL, this );
	m_txtYstar ->Disconnect( wxEVT_COMMAND_TEXT_UPDATED   , wxCommandEventHandler( PatternCenterPanel::YstarChanged   ), NULL, this );
	m_txtL     ->Disconnect( wxEVT_COMMAND_TEXT_UPDATED   , wxCommandEventHandler( PatternCenterPanel::LChanged       ), NULL, this );
	m_txtZstar ->Disconnect( wxEVT_COMMAND_TEXT_UPDATED   , wxCommandEventHandler( PatternCenterPanel::ZstarChanged   ), NULL, this );
	m_txtTlt   ->Disconnect( wxEVT_COMMAND_TEXT_UPDATED   , wxCommandEventHandler( PatternCenterPanel::TltChanged     ), NULL, this );
	m_btnFit   ->Disconnect( wxEVT_COMMAND_BUTTON_CLICKED , wxCommandEventHandler( PatternCenterPanel::DoFit          ), NULL, this );
}

#endif//_PAT_CEN_PAN_H_
