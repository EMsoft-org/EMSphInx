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

#ifndef _WIDSOM_PROMPT_H_
#define _WIDSOM_PROMPT_H_

#include <wx/string.h>
#include <wx/stattext.h>
#include <wx/textctrl.h>
#include <wx/button.h>
#include <wx/sizer.h>
#include <wx/dialog.h>

class WisdomPrompt : public wxDialog {
	private:
		bool m_valid = true;

	protected:
		wxTextCtrl* m_txtBw;
		wxButton  * m_btn  ;

		//@brief: sanity check text box on button click
		virtual void OnCompute( wxCommandEvent& event );

	public:

		WisdomPrompt( wxWindow* parent, wxWindowID id = wxID_ANY, const wxString& title = wxT("FFT Wisdom Builder"), const wxPoint& pos = wxDefaultPosition, const wxSize& size = wxSize( 375,185 ), long style = wxDEFAULT_DIALOG_STYLE|wxSYSTEM_MENU );
		~WisdomPrompt();

		//@brief : parse bandwidths
		//@return: list of selected bandwidths
		std::vector<size_t> getBandwidths();


};

//@brief: sanity check text box on button click
void WisdomPrompt::OnCompute( wxCommandEvent& event ) {
	try {
		std::vector<size_t> bw = getBandwidths();
		EndModal(wxID_OK);
	} catch (std::runtime_error& e) {
		wxMessageDialog msgDlg(this, e.what(), "Error Parsing Bandwidths");	
		msgDlg.ShowModal();	
	}
}

WisdomPrompt::WisdomPrompt( wxWindow* parent, wxWindowID id, const wxString& title, const wxPoint& pos, const wxSize& size, long style ) : wxDialog( parent, id, title, pos, size, style ) {

	//build sizers
	this->SetSizeHints( wxDefaultSize, wxDefaultSize );
	wxBoxSizer* vSizer = new wxBoxSizer( wxVERTICAL );
	wxBoxSizer* hSizer = new wxBoxSizer( wxHORIZONTAL );

	//build/configure elements
	wxStaticText* staTxt = new wxStaticText( this, wxID_ANY, wxT("Enter bandwidth(s) to compute FFT wisdom for. Multiple bandwidths can be seperated with a comma and ranges can be specified with a dash. For example “53, 63-68, 74, 88” will compute FFT plans for bandwidths 53, 63, 64, 65, 66, 67, 68, 74, and 88."), wxDefaultPosition, wxDefaultSize, 0 );
	m_txtBw = new wxTextCtrl( this, wxID_ANY, wxT("53, 63, 68, 74, 88, 95, 113, 122, 158, 172, 188, 203, 221, 263"), wxDefaultPosition, wxDefaultSize, wxTE_CENTER|wxTE_MULTILINE );
	m_btn = new wxButton( this, wxID_ANY, wxT("Compute"), wxDefaultPosition, wxDefaultSize, 0 );//any instead of wxID_OK to make overriding close behavoir easier
	staTxt->Wrap( 350 );
	m_btn->SetDefault();

	//assemble
	hSizer->Add( 5, 0, 0, wxEXPAND, 5 );
	hSizer->Add( m_txtBw, 1, wxALL, 5 );
	hSizer->Add( 5, 0, 0, wxEXPAND, 5 );
	hSizer->Add( m_btn, 0, wxALL|wxALIGN_CENTER_VERTICAL, 5 );
	hSizer->Add( 5, 0, 0, wxEXPAND, 5 );
	vSizer->Add( staTxt, 0, wxALL|wxALIGN_CENTER_HORIZONTAL, 5 );
	vSizer->Add( hSizer, 0, wxEXPAND, 5 );

	//layout
	this->SetSizer( vSizer );
	this->Layout();
	this->Centre( wxBOTH );

	//connect Events
	m_btn->Connect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( WisdomPrompt::OnCompute ), NULL, this );
}

WisdomPrompt::~WisdomPrompt() {
	//disconnect Events
	m_btn->Disconnect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( WisdomPrompt::OnCompute ), NULL, this );
}

#include <sstream>
#include <algorithm>

//@brief : parse bandwidths
//@return: list of selected bandwidths
std::vector<size_t> WisdomPrompt::getBandwidths() {
	//create empty vector and wrap text as stringstream for parsing
	std::vector<size_t> ret(0);//initially no values
	std::istringstream ss(m_txtBw->GetValue().ToStdString());
	if(m_txtBw->GetValue().IsEmpty()) throw std::runtime_error("empty bandwidth string");
	ss >> std::skipws;//ignore whitespace

	//parse first value
	char c;
	size_t v;
	size_t idx = 0;
	if(!(ss >> v)) throw std::runtime_error("failed to parse any numbers from bandwidth string");
	ret.push_back(v);
	idx = ss.tellg();

	//keep parsing while there is string left
	while(ss >> c) {//get next seperator (should be ',' or '-')
		bool range = '-' == c;//check if this is a range (bandwidths are non-negative so we don't have to worry about ambiguity)
		idx = ss.tellg();//upate position
		if(!(ss >> v)) throw std::runtime_error("failed to parse number from `" + ss.str().substr(idx) + '\'');//try to extract the numeric value
		if(range) {//v is the upper bound of a range
			if(v <= ret.back()) throw std::runtime_error("upper bound of range <= lower bound");
			for(size_t i = ret.back() + 1; i <= v; i++) ret.push_back(i);//add all numbers in range
		} else {//v is a single number
			ret.push_back(v);
		}
	}
	std::sort(ret.begin(), ret.end());//sort bandwidths
	ret.erase( std::unique( ret.begin(), ret.end() ), ret.end() );//remove duplicates
	if(ret.front() < 16 || ret.back() > 384) throw std::runtime_error("bandwidths should be in [16,384]");

	return ret;
}

#endif//_WIDSOM_PROMPT_H_


