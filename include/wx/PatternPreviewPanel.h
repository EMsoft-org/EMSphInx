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

#ifndef _PAT_PREVIEW_H_
#define _PAT_PREVIEW_H_

#include <wx/statbmp.h>
#include <wx/sizer.h>
#include <wx/scrolbar.h>
#include <wx/stattext.h>
#include <wx/spinctrl.h>
#include <wx/frame.h>
#include <wx/checkbox.h>

#include "ImagePanel.h"

#include "modality/ebsd/imprc.hpp"

//@brief: frame for pattern preview
class PatternPreviewPanel : public wxPanel {
	std::shared_ptr< std::vector< std::vector<char> > > m_images ;
	emsphinx::ebsd::PatternProcessor<float>             m_proc   ;//image processor
	size_t                                              m_idxCur ;//current image index
	int                                                 m_circCur;//current circular radius
	bool                                                m_bckgCur;//current background flag
	size_t                                              m_nRegCur;//current nregions

	//manually catch arrow keys (they don't seem to be captured on osx)
	void keyLeft ();
	void keyRight();
	void keyUp   ();
	void keyDown ();
	void keySpace();

	//@brief     : update the the displayed images
	//@param idx : new index to display (must be < m_images->size())
	//@param circ: circular mask radius
	//@param bckg: new background subtraction flag
	//@param nReg: new nRegions to display
	void updateImages(size_t idx, int circ, bool bckg, size_t nReg);

	protected:
		wxImagePanel * m_panelRaw ;
		wxImagePanel * m_panelPrc ;
		wxScrollBar  * m_scrollBar;
		wxSpinCtrl   * m_spinCtlR ;
		wxCheckBox   * m_bckgChk  ;
		wxSpinCtrl   * m_spinCtlNr;

		// Virtual event handlers, overide them in your derived class
		virtual void scrlPat ( wxScrollEvent & event) { updateImages(GetIdx(), m_circCur, m_bckgCur, m_nRegCur); }//@brief: triggered by scroll bar changing [change displayed pattern]
		virtual void circPat ( wxSpinEvent   & event) { updateImages(m_idxCur, GetCirc(), m_bckgCur, m_nRegCur); }//@brief: triggered by spinner changing [change mask radius]
		virtual void bckgPat ( wxCommandEvent& event) { updateImages(m_idxCur, m_circCur, GetBckg(), m_nRegCur); }//@brief: triggered by check box changing [change background subtraction]
		virtual void procPat ( wxSpinEvent   & event) { updateImages(m_idxCur, m_circCur, m_bckgCur, GetNreg()); }//@brief: triggered by spinner changing [change nregions]
		// virtual void procChk ( wxCommandEvent& event);//triggered by check box changing [changing circ mask]
		        void keyPress( wxKeyEvent    & event);
		
		DECLARE_EVENT_TABLE()

	public:

		PatternPreviewPanel( wxWindow* parent, wxWindowID id = wxID_ANY, const wxPoint& pos = wxDefaultPosition, const wxSize& size = wxSize( 600,400 ), long style = wxDEFAULT_FRAME_STYLE|wxTAB_TRAVERSAL );
		~PatternPreviewPanel();

		long GetIdx () const {return m_scrollBar->GetThumbPosition();}
		long GetCirc() const {return m_spinCtlR ->GetValue        ();}
		bool GetBckg() const {return m_bckgChk  ->GetValue        ();}
		long GetNreg() const {return m_spinCtlNr->GetValue        ();}
		void SetIdx (const long v) {m_scrollBar->SetThumbPosition(v); updateImages(v       , m_circCur, m_bckgCur, m_nRegCur);}
		void SetCirc(const long v) {m_spinCtlR ->SetValue        (v); updateImages(m_idxCur, v        , m_bckgCur, m_nRegCur);}
		void SetBckg(const bool v) {m_bckgChk  ->SetValue        (v); updateImages(m_idxCur, m_circCur, v        , m_nRegCur);}
		void SetNreg(const long v) {m_spinCtlNr->SetValue        (v); updateImages(m_idxCur, m_circCur, m_bckgCur, v        );}

		void ForceRedraw() {m_idxCur = -1; updateImages(GetIdx(), m_circCur, m_bckgCur, m_nRegCur);}

		void SetImages(std::shared_ptr< std::vector< std::vector<char> > > im, const size_t w, const size_t h);

};

//there are a bunch of scroll events, so this is easier
BEGIN_EVENT_TABLE(PatternPreviewPanel, wxPanel)
	EVT_SCROLL    (PatternPreviewPanel::scrlPat )
END_EVENT_TABLE()

///////////////////////////////////////////////////////////////////////////
//               PatternPreviewPanel logic Implementation                //
///////////////////////////////////////////////////////////////////////////

//@brief     : update the the displayed images
//@param idx : new index to display (must be < m_images->size())
//@param circ: circular mask radius
//@param bckg: new background subtraction flag
//@param nReg: new nRegions to display
void PatternPreviewPanel::updateImages(size_t idx, int circ, bool bckg, size_t nReg) {
	//handle insufficent images
	if(NULL == m_images.get() ? true : idx >= m_images->size()) {
		wxImage empty;
		m_panelRaw->setImage(empty);
		m_panelPrc->setImage(empty);
		return;
	}
	wxSize sz = m_panelRaw->GetSize();

	//start by updating the image processing paramters (builds a new mask)
	const bool imProcChanged = m_circCur != circ || m_bckgCur != bckg|| m_nRegCur != nReg;
	if(imProcChanged) {
		m_proc.setSize(sz.GetWidth(), sz.GetHeight(), circ, bckg, nReg);//update image processing + mask
		if(m_circCur != circ) {//the mask radius changed, update alpha channels
			char const * ptr = m_proc.getMask();
			const size_t nPix = sz.GetWidth() * sz.GetHeight();
			std::transform(ptr, ptr + nPix, m_panelRaw->GetAlpha(), [](const char& c){return 0 == c ? 0x00 : 0xFF;});
			std::transform(ptr, ptr + nPix, m_panelPrc->GetAlpha(), [](const char& c){return 0 == c ? 0x00 : 0xFF;});
			m_panelRaw->markStale();
			m_panelPrc->markStale();
			m_panelRaw->Refresh();
		}
		m_bckgCur = bckg;
		m_nRegCur = nReg;
		m_circCur = circ;
	}

	//first update the raw image if needed
	static bool first = true;
	if(m_idxCur != idx || first) {//the input image changed
		//loop over image updating
		char const * pBuff = m_images->operator[](idx).data();
		for(int j = 0; j < sz.GetHeight(); j++) {
			for(int i = 0; i < sz.GetWidth(); i++) {
				m_panelRaw->SetGray(i, j, *pBuff);
				++pBuff;
			}
		}
		m_panelRaw->Refresh();
		m_idxCur = idx;
		first = true;//force update on processed image even if nReg is the same
	}

	//next update the processed image if needed
	if(imProcChanged || first) {//either the input image changed, this is the first pass, or the image processing changed
		std::vector<char> ahe = m_images->operator[](idx);//copy input image
		m_proc.process((uint8_t*)ahe.data());//do AHE

		//update processed image
		char const * pBuff = ahe.data();
		for(int j = 0; j < sz.GetHeight(); j++) {
			for(int i = 0; i < sz.GetWidth(); i++) {
				m_panelPrc->SetGray(i, j, *pBuff);
				++pBuff;
			}
		}
		m_panelPrc->Refresh();
	}
	first = false;
}

void PatternPreviewPanel::SetImages(std::shared_ptr< std::vector< std::vector<char> > > img, const size_t w, const size_t h) {
	//update image sizes
	wxImage im(w, h);
	m_panelRaw->setImage(im);
	m_panelPrc->setImage(im);

	//update the spin control if needed
	const size_t minDim = std::min(w, h);
	const double diag = std::hypot(double(w) / 2, double(h) / 2);
	m_spinCtlR ->SetRange(-1, (int)std::ceil( diag ) );//don't allow too many regions (at least 3x3 window for each histogram)
	m_spinCtlNr->SetRange( 0, std::min<size_t>(32, minDim / 3));//don't allow too many regions (at least 3x3 window for each histogram)

	//update the scroll bar
	int idx = 0;
	m_images = img;

	int pos  = idx;
	int thmSz = m_images->size() > 100 ? m_images->size() / 100 : 1;
	int rng   = m_images->size();
	int pgSz  = std::min(thmSz * 10, rng);
	m_scrollBar->SetScrollbar(idx, thmSz, rng, pgSz);
	m_scrollBar->Enable(rng > 1);

	//get the minimum size of the panel (from contrl bar) and update minimum height to reflect images
	const wxSize& sz = this->GetSizer()->GetMinSize();
	double numIm = double(sz.GetWidth()) / (2 * w);//this is how many times we could fit 2 full sized images size by size
	int imH = std::round(numIm * h) + 10;//this is how tall the images should be for the widhts to nicely fill the frame
	this->GetSizer()->SetMinSize(wxSize(sz.GetWidth(), sz.GetHeight() + imH));//make images fill the frame
	this->SetMinSize(wxSize(sz.GetWidth(), sz.GetHeight() + imH));//make images fill the frame
	this->GetSizer()->Fit(GetParent());//force the parent window to grow to accomodate our new size

	//update displayed image
	m_circCur = -2;//force mask recalculation
	updateImages(idx, GetCirc(), GetBckg(), GetNreg());
}

///////////////////////////////////////////////////////////////////////////
//                 PatternPreviewPanel wxImplementation                  //
///////////////////////////////////////////////////////////////////////////

//@brief: decrement the scroll bar and update images
void PatternPreviewPanel::keyLeft () {
	if( !m_scrollBar->IsEnabled()) return;
	if( !m_scrollBar ->HasFocus() ) {m_panelRaw->SetFocus(); m_scrollBar ->SetFocus();}//scorr bar doesn't take focus from spin ctrl on osx
	m_scrollBar ->SetThumbPosition(m_scrollBar->GetThumbPosition() - 1);
	wxScrollEvent evt;
	scrlPat(evt);
}

//@brief: increment the scroll bar and update images
void PatternPreviewPanel::keyRight() {
	if( !m_scrollBar->IsEnabled()) return;
	if( !m_scrollBar ->HasFocus() ) {m_panelRaw->SetFocus(); m_scrollBar ->SetFocus();}//scorr bar doesn't take focus from spin ctrl on osx
	m_scrollBar ->SetThumbPosition(m_scrollBar->GetThumbPosition() + 1);
	wxScrollEvent evt;
	scrlPat(evt);
}

//@brief: increment the spinner and update images
void PatternPreviewPanel::keyUp   () {
	wxSpinEvent evt;
	if( !m_spinCtlNr->HasFocus() ) {
		if(m_spinCtlR->HasFocus()) {//dont steal focus from radius spin control
			m_spinCtlR->SetValue(m_spinCtlR->GetValue() + 1);
			circPat(evt);
			return;
		}
		m_spinCtlNr->SetFocus();
	}
	m_spinCtlNr->SetValue(m_spinCtlNr->GetValue() + 1);
	procPat(evt);
}

//@brief: decrement the spinner and update images
void PatternPreviewPanel::keyDown () {
	wxSpinEvent evt;
	if( !m_spinCtlNr->HasFocus() ) {
		if(m_spinCtlR->HasFocus()) {//dont steal focus from radius spin control
			m_spinCtlR->SetValue(m_spinCtlR->GetValue() - 1);
			circPat(evt);
			return;
		}
		m_spinCtlNr->SetFocus();
	}
	m_spinCtlNr->SetValue(m_spinCtlNr->GetValue() - 1);
	procPat(evt);
}

//@brief: toggle check box on space bar
void PatternPreviewPanel::keySpace () {
	m_bckgChk->SetValue(!m_bckgChk->GetValue());
	wxCommandEvent evt;
	bckgPat(evt);
}

//@brief: handle arrows globally + auto switch focus
void PatternPreviewPanel::keyPress( wxKeyEvent    & event) { 
	switch(event.GetKeyCode()) {
		case WXK_LEFT : keyLeft (); break;
		case WXK_RIGHT: keyRight(); break;
		case WXK_UP   : keyUp   (); break;
		case WXK_DOWN : keyDown (); break;
		case WXK_SPACE: keySpace(); break;
		default: event.Skip();
	}
}

PatternPreviewPanel::PatternPreviewPanel( wxWindow* parent, wxWindowID id, const wxPoint& pos, const wxSize& size, long style ) : wxPanel( parent, id, pos, size, style ) {
	this->SetSizeHints( wxDefaultSize, wxDefaultSize );

	//build sizers
	wxBoxSizer* vSizer     = new wxBoxSizer( wxVERTICAL   );//main sizer is vertical boxes
	wxBoxSizer* hSizerIm   = new wxBoxSizer( wxHORIZONTAL );//horizontal box for 2 images
	wxBoxSizer* hSizerCtrl = new wxBoxSizer( wxHORIZONTAL );//horizontal box for spin control

	//build elements
	m_panelRaw            = new wxImagePanel  ( this, wxID_ANY      , wxDefaultPosition            , wxDefaultSize                                                  );//raw image
	m_panelPrc            = new wxImagePanel  ( this, wxID_ANY      , wxDefaultPosition            , wxDefaultSize                                                  );//raw image
	m_scrollBar           = new wxScrollBar   ( this, wxID_ANY      , wxDefaultPosition            , wxDefaultSize    , wxSB_HORIZONTAL                             );//horizontal scroll bar
	wxStaticText* txtR    = new wxStaticText  ( this, wxID_ANY      , wxT("Mask Radius")           , wxDefaultPosition, wxDefaultSize  , 0                          );//spinner label
	m_spinCtlR            = new wxSpinCtrl    ( this, wxID_ANY      , wxEmptyString                , wxDefaultPosition, wxDefaultSize  , wxSP_ARROW_KEYS, -1, 0 , -1);
	m_bckgChk             = new wxCheckBox    ( this, wxID_ANY      , wxT("Gaussian Background")   , wxDefaultPosition, wxDefaultSize                               );//check box
	wxStaticText* txtNReg = new wxStaticText  ( this, wxID_ANY      , wxT("Histogram Equalization"), wxDefaultPosition, wxDefaultSize  , 0                          );//spinner label
	m_spinCtlNr           = new wxSpinCtrl    ( this, wxID_ANY      , wxEmptyString                , wxDefaultPosition, wxDefaultSize  , wxSP_ARROW_KEYS,  0, 32,  4);//spinner (initial value 4)

	//build image box
	hSizerIm->Add( m_panelRaw, 1, wxALL|wxEXPAND, 5 );
	hSizerIm->Add( m_panelPrc, 1, wxALL|wxEXPAND, 5 );

	//build controls box
	hSizerCtrl->AddStretchSpacer();
	hSizerCtrl->Add( txtR       , 0, wxALL|wxALIGN_CENTER_VERTICAL, 5 );
	hSizerCtrl->Add( m_spinCtlR , 0, wxALL|wxALIGN_CENTER_VERTICAL, 5 );
	hSizerCtrl->AddStretchSpacer();
	hSizerCtrl->Add( m_bckgChk  , 0, wxALL|wxALIGN_CENTER_VERTICAL, 5 );
	hSizerCtrl->AddStretchSpacer();
	hSizerCtrl->Add( txtNReg    , 0, wxALL|wxALIGN_CENTER_VERTICAL, 5 );
	hSizerCtrl->Add( m_spinCtlNr, 0, wxALL|wxALIGN_CENTER_VERTICAL, 5 );
	hSizerCtrl->AddStretchSpacer();

	//assemble into overall frame
	vSizer->Add( hSizerIm   , 1, wxEXPAND      , 5 );
	vSizer->Add( m_scrollBar, 0, wxALL|wxEXPAND, 5 );
	vSizer->Add( hSizerCtrl , 0, wxEXPAND      , 5 );

	//assign sizer to frame
	this->SetSizer( vSizer );
	this->Layout();
	this->Centre( wxBOTH );

	m_spinCtlR ->Connect( wxEVT_COMMAND_SPINCTRL_UPDATED, wxSpinEventHandler   ( PatternPreviewPanel::circPat ), NULL, this );
	m_bckgChk  ->Connect( wxEVT_CHECKBOX                , wxCommandEventHandler( PatternPreviewPanel::bckgPat ), NULL, this );
	m_spinCtlNr->Connect( wxEVT_COMMAND_SPINCTRL_UPDATED, wxSpinEventHandler   ( PatternPreviewPanel::procPat ), NULL, this );
	this       ->Connect( wxEVT_CHAR_HOOK               , wxKeyEventHandler    ( PatternPreviewPanel::keyPress), NULL, this );
}

PatternPreviewPanel::~PatternPreviewPanel() {
	m_spinCtlR ->Disconnect( wxEVT_COMMAND_SPINCTRL_UPDATED, wxSpinEventHandler   ( PatternPreviewPanel::circPat ), NULL, this );
	m_bckgChk  ->Disconnect( wxEVT_CHECKBOX                , wxCommandEventHandler( PatternPreviewPanel::bckgPat ), NULL, this );
	m_spinCtlNr->Disconnect( wxEVT_COMMAND_SPINCTRL_UPDATED, wxSpinEventHandler   ( PatternPreviewPanel::procPat ), NULL, this );
	this       ->Disconnect( wxEVT_CHAR_HOOK               , wxKeyEventHandler    ( PatternPreviewPanel::keyPress), NULL, this );
}

#endif//_PAT_PREVIEW_H_