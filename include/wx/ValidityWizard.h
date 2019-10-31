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

#ifndef _VALID_WIZ_H_
#define _VALID_WIZ_H_

#include <wx/artprov.h>
#include <wx/xrc/xmlres.h>
#include <wx/panel.h>
#include <wx/gdicmn.h>
#include <wx/font.h>
#include <wx/colour.h>
#include <wx/settings.h>
#include <wx/string.h>
#include <wx/simplebook.h>
#include <wx/statline.h>
#include <wx/bitmap.h>
#include <wx/image.h>
#include <wx/icon.h>
#include <wx/button.h>
#include <wx/sizer.h>
#include <wx/statusbr.h>
#include <wx/frame.h>

#include "wx/ValidityPanel.h"

///////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
/// Class ValidityWizard
///////////////////////////////////////////////////////////////////////////////
class ValidityWizard : public wxFrame
{
	private:
		bool m_finished = false;

	protected:
		wxSimplebook* m_book;
		wxSimplebook* m_ctrlBook;
		wxPanel* m_ctrlPanel;
		wxButton* m_btnPrev;
		wxButton* m_btnNext;
		wxButton* m_btnCancel;
		wxStatusBar* m_statusBar;

		// Virtual event handlers, overide them in your derived class
		virtual void PageChanged( wxBookCtrlEvent& event );
		virtual void PageChanging( wxBookCtrlEvent& event ) {event.Skip();}// if(!((ValidityPanel*)m_book->GetPage(event.GetOldSelection()))->isValid()) event.Veto(); }

		virtual void OnBack( wxCommandEvent& event );
		virtual void OnNext( wxCommandEvent& event );

		virtual void OnCancel( wxCommandEvent& event ) { Close(); }
		virtual void ValidityChanged( ValidityChangedEvent& event );

		void AddPage(ValidityPanel* panel);

	public:


		ValidityWizard( wxWindow* parent, wxWindowID id = wxID_ANY, const wxString& title = wxT("My Wizard"), const wxPoint& pos = wxDefaultPosition);

		~ValidityWizard();

		bool isFinished() {return m_finished;}

};


///////////////////////////////////////////////////////////////////////////
// C++ code generated with wxFormBuilder (version Aug 31 2019)
// http://www.wxformbuilder.org/
//
// PLEASE DO *NOT* EDIT THIS FILE!
///////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////

void ValidityWizard::PageChanged( wxBookCtrlEvent& event ) {
	m_btnPrev->Enable(0 != event.GetSelection());//disable the previous button on first page
	m_btnNext->SetLabel( 1 + event.GetSelection() == m_book->GetPageCount() ? "Finish" : "Next >");//update next button on last page
	ValidityPanel* pan = (ValidityPanel*)m_book->GetPage(event.GetSelection());
	m_btnNext->Enable(pan->isValid());//
	m_statusBar->SetStatusText(pan->validMsg());
}

void ValidityWizard::OnBack( wxCommandEvent& event ) {
	m_book->SetSelection(m_book->GetSelection() - 1);
}

void ValidityWizard::OnNext( wxCommandEvent& event ) {
	int idx = m_book->GetSelection() + 1;//get index of next page
	int num = m_book->GetPageCount();
	if(idx == num) {//handle the last page specially
		m_finished = true;
		Close();//finish button
	} else {
		if(idx + 1 == num) m_btnNext->SetLabel("Finish");//update label on last page
		m_book->SetSelection(idx);//move to current page

	}
}

void ValidityWizard::ValidityChanged( ValidityChangedEvent& event ) {
	if(event.GetEventObject() == (wxObject*)m_book->GetPage(m_book->GetSelection())) {//make sure the current page emmited the event
		m_btnNext->Enable(event.GetValid());
		m_statusBar->SetStatusText(event.GetMsg());
	}
}

void ValidityWizard::AddPage(ValidityPanel* panel) {
	if(0 == m_book->GetPageCount()) m_statusBar->SetStatusText(panel->validMsg());
	m_book->AddPage(panel, wxT(""), false );
	wxSize minSize = m_book->GetPage(0)->GetBestSize();
	for(size_t i = 1; i < m_book->GetPageCount(); i++) {
		wxSize pgSize = m_book->GetPage(i)->GetBestSize();
		minSize.x = std::max(minSize.x, pgSize.x);
		minSize.y = std::max(minSize.y, pgSize.y);
	}
	m_book->SetPageSize(minSize);
	panel->Connect( VALID_CHANGED, ValidityChangedEventHandler( ValidityWizard::ValidityChanged ), NULL, this );
}

ValidityWizard::ValidityWizard( wxWindow* parent, wxWindowID id, const wxString& title, const wxPoint& pos) : wxFrame( parent, id, title, pos, wxSize( 325,350 ), wxCAPTION|wxCLOSE_BOX|wxSYSTEM_MENU|wxTAB_TRAVERSAL )
{
	if(NULL != parent) parent->Enable(false);//approximate modal behavoir

	this->SetSizeHints( wxDefaultSize, wxDefaultSize );

	wxBoxSizer* bSizer;
	bSizer = new wxBoxSizer( wxVERTICAL );

	m_book = new wxSimplebook( this, wxID_ANY, wxDefaultPosition, wxDefaultSize, 0 );

	bSizer->Add( m_book, 1, wxEXPAND, 5 );

	m_ctrlBook = new wxSimplebook( this, wxID_ANY, wxDefaultPosition, wxDefaultSize, 0 );
	m_ctrlPanel = new wxPanel( m_ctrlBook, wxID_ANY, wxDefaultPosition, wxDefaultSize, wxTAB_TRAVERSAL );
	wxBoxSizer* vSizer;
	vSizer = new wxBoxSizer( wxVERTICAL );

	wxStaticLine* staticLine = new wxStaticLine( m_ctrlPanel, wxID_ANY, wxDefaultPosition, wxDefaultSize, wxLI_HORIZONTAL );
	vSizer->Add( staticLine, 0, wxEXPAND|wxALL, 5 );

	wxBoxSizer* hSizer;
	hSizer = new wxBoxSizer( wxHORIZONTAL );


	hSizer->Add( 0, 0, 1, wxEXPAND, 5 );

	m_btnPrev = new wxButton( m_ctrlPanel, wxID_ANY, wxT("< Back"), wxDefaultPosition, wxDefaultSize, 0 );
	hSizer->Add( m_btnPrev, 0, wxALL, 5 );

	m_btnNext = new wxButton( m_ctrlPanel, wxID_ANY, wxT("Next >"), wxDefaultPosition, wxDefaultSize, 0 );
	hSizer->Add( m_btnNext, 0, wxALL, 5 );

	m_btnCancel = new wxButton( m_ctrlPanel, wxID_ANY, wxT("Cancel"), wxDefaultPosition, wxDefaultSize, 0 );
	hSizer->Add( m_btnCancel, 0, wxALL, 5 );


	vSizer->Add( hSizer, 1, wxEXPAND, 5 );


	m_ctrlPanel->SetSizer( vSizer );
	m_ctrlPanel->Layout();
	vSizer->Fit( m_ctrlPanel );
	m_ctrlBook->AddPage( m_ctrlPanel, wxT("a page"), false );

	bSizer->Add( m_ctrlBook, 0, wxEXPAND, 5 );


	this->SetSizer( bSizer );
	this->Layout();
	m_statusBar = this->CreateStatusBar( 1, wxSTB_ELLIPSIZE_END|wxSTB_SHOW_TIPS, wxID_ANY );

	this->Centre( wxBOTH );

	//disable both buttons and make next button defaultx
	m_btnPrev->Disable();
	m_btnNext->Disable();
	m_btnNext->SetDefault();

	// Connect Events
	m_book->Connect( wxEVT_COMMAND_BOOKCTRL_PAGE_CHANGED, wxBookCtrlEventHandler( ValidityWizard::PageChanged ), NULL, this );
	m_book->Connect( wxEVT_COMMAND_BOOKCTRL_PAGE_CHANGING, wxBookCtrlEventHandler( ValidityWizard::PageChanging ), NULL, this );
	m_btnPrev->Connect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( ValidityWizard::OnBack ), NULL, this );
	m_btnNext->Connect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( ValidityWizard::OnNext ), NULL, this );
	m_btnCancel->Connect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( ValidityWizard::OnCancel ), NULL, this );
}

ValidityWizard::~ValidityWizard()
{
	if(NULL != m_parent) m_parent->Enable(true);//approximate modal behavoir

	// Disconnect Events
	m_book->Disconnect( wxEVT_COMMAND_BOOKCTRL_PAGE_CHANGED, wxBookCtrlEventHandler( ValidityWizard::PageChanged ), NULL, this );
	m_book->Disconnect( wxEVT_COMMAND_BOOKCTRL_PAGE_CHANGING, wxBookCtrlEventHandler( ValidityWizard::PageChanging ), NULL, this );
	m_btnPrev->Disconnect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( ValidityWizard::OnBack ), NULL, this );
	m_btnNext->Disconnect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( ValidityWizard::OnNext ), NULL, this );
	m_btnCancel->Disconnect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( ValidityWizard::OnCancel ), NULL, this );


	for(size_t i = 0; i < m_book->GetPageCount(); i++) {
		m_book->GetPage(i)->Disconnect( VALID_CHANGED, ValidityChangedEventHandler( ValidityWizard::ValidityChanged ), NULL, this );
	}
	m_book->DeleteAllPages();
}

#endif//_VALID_WIZ_H_
