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

#ifndef _PER_TAB_PAN_H_
#define _PER_TAB_PAN_H_

#include <bitset>
#include <string>

//@brief: helper class to store a bitmask of elements
//@note : bits are 1 indexed to e.g. He is at .test(2) not .test(1)
struct ElementMask : public std::bitset<128> {//we only need [1,118]
	//@brief : convert to string list of elements
	//@return: string representation, e.g. "H, He, Li"
	std::string str() const;

	//@brief  : get the symbol for an element
	//@param z: atomic number
	//@return : e.g. "Li" for 3
	static std::string Symbol(const size_t z);
};

#include <wx/image.h>
#include <wx/button.h>
#include <wx/tglbtn.h>
#include <wx/sizer.h>
#include <wx/panel.h>

//@brief: panel to display a periodic table of toggle buttons
class PeriodicTablePanel : public wxPanel {
	std::vector<wxToggleButton*> m_btns;//a button for each element in order (0 indexed)

	public:
		//@brief : get a mask of which elements are toggled
		//@return: bitmask of set elements
		ElementMask getMask() const;

		//@brief   : set which elements are toggled
		//@param el: bitmask of toggle states
		void setMask(ElementMask el);

		PeriodicTablePanel( wxWindow* parent, wxWindowID id = wxID_ANY, const wxPoint& pos = wxDefaultPosition, const wxSize& size = wxDefaultSize, long style = wxTAB_TRAVERSAL, const wxString& name = wxEmptyString );
		~PeriodicTablePanel();

};

///////////////////////////////////////////////////////////////////////////

//@brief : get a mask of which elements are toggled
//@return: bitmask of set elements
ElementMask PeriodicTablePanel::getMask() const {
	ElementMask msk;
	if(118 != m_btns.size()) throw std::logic_error("unexpected number of elements in periodic table panel");
	for(size_t i = 0; i < 118; i++) {
		if(m_btns[i]->GetValue()) msk.set(i+1);
	}
	return msk;
}

//@brief   : set which elements are toggled
//@param el: bitmask of toggle states
void PeriodicTablePanel::setMask(ElementMask el) {
	if(118 != m_btns.size()) throw std::logic_error("unexpected number of elements in periodic table panel");
	for(size_t i = 0; i < 118; i++) m_btns[i]->SetValue(el.test(i+1));
}

PeriodicTablePanel::PeriodicTablePanel( wxWindow* parent, wxWindowID id, const wxPoint& pos, const wxSize& size, long style, const wxString& name ) : wxPanel( parent, id, pos, size, style, name ) {
	wxGridSizer* gSizer  = new wxGridSizer( 10, 18, 0, 0 );
	this  ->SetMinSize(wxSize(18 * 40, 10 * 50));
	gSizer->SetMinSize(wxSize(18 * 40, 10 * 50));

	//element to label
	std::function<std::string(size_t, bool)> labelFunc = [](const size_t z, const bool mu){
		std::ostringstream ss;
		if(mu)
			ss << "<big><b>" << ElementMask::Symbol(z) << "</b></big>\n<small>" << z << "</small>";
		else
			ss               << ElementMask::Symbol(z) <<           "\n"        << z              ;//msw doesn't support multiline markup
		return ss.str();
	};

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
	const bool markup = false;
#else
	const bool markup = true ;
#endif

	//build all buttons up front
	for(size_t i = 1; i <= 118; i++) {
		std::ostringstream ss;
		m_btns.push_back( new wxToggleButton( this, wxID_ANY, labelFunc(i, false), wxDefaultPosition, wxDefaultSize, wxBU_EXACTFIT) );
		if(markup) m_btns.back()->SetLabelMarkup(labelFunc(i, markup));
	}

	//build special disabled buttons for rare earth breaks
	wxToggleButton* reBtn[2] = {
		new wxToggleButton( this, wxID_ANY, labelFunc(57, false), wxDefaultPosition, wxDefaultSize, wxBU_EXACTFIT),
		new wxToggleButton( this, wxID_ANY, labelFunc(89, false), wxDefaultPosition, wxDefaultSize, wxBU_EXACTFIT)
	};
	if(markup) reBtn[0]->SetLabelMarkup(labelFunc(57, markup));
	if(markup) reBtn[1]->SetLabelMarkup(labelFunc(89, markup));
	reBtn[0]->Disable();
	reBtn[1]->Disable();

	//add rows without rare earth breaks
	const int brd = 2;//border
	for(size_t i = 0; i < 54; i++) {
		gSizer->Add( m_btns[i], 1, wxALL|wxEXPAND, brd );
		if(i == 0 ) for(size_t i = 0; i < 16; i++) gSizer->AddStretchSpacer();
		if(i == 3 ) for(size_t i = 0; i < 10; i++) gSizer->AddStretchSpacer();
		if(i == 11) for(size_t i = 0; i < 10; i++) gSizer->AddStretchSpacer();
	}

	//add rows with rare earth breaks
	for(size_t j = 0; j < 2; j++) {
		size_t start = 54 + j * 32;
		for(size_t i = 0; i <  2; i++) gSizer->Add( m_btns[start+i   ], 1, wxALL|wxEXPAND, brd );//before break
			                           gSizer->Add( reBtn [      j   ], 1, wxALL|wxEXPAND, brd );//break
		for(size_t i = 3; i < 18; i++) gSizer->Add( m_btns[start+i+14], 1, wxALL|wxEXPAND, brd );//after break
	}

	//add rare earth rows
	for(size_t i = 0; i < 18; i++) gSizer->AddStretchSpacer();//empty row
	for(size_t j = 0; j < 2; j++) {
		size_t start = 54 + j * 32;
		for(size_t i = 0; i < 2; i++) gSizer->AddStretchSpacer();
		for(size_t i = 2; i < 17; i++) gSizer->Add( m_btns[start+i], 1, wxALL|wxEXPAND, brd );//after break
		gSizer->AddStretchSpacer();
	}

	this->SetSizer( gSizer );
	gSizer->Fit(this);
	this->Layout();
}

PeriodicTablePanel::~PeriodicTablePanel() {
}

////////////////////////////////////////////////////////////////////////
//                            ElementMask                             //
////////////////////////////////////////////////////////////////////////

#include <vector>

//@brief  : get the symbol for an element
//@param z: atomic number
//@return : e.g. "Li" for 3
std::string ElementMask::Symbol(const size_t z) {
	static const std::vector<std::string> syms = {
		"H" , "He", "Li", "Be", "B" , "C" , "N" , "O" , "F" , "Ne", "Na", "Mg", "Al", "Si", "P" , "S" ,
		"Cl", "Ar", "K" , "Ca", "Sc", "Ti", "V" , "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge",
		"As", "Se", "Br", "Kr", "Rb", "Sr", "Y" , "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
		"In", "Sn", "Sb", "Te", "I" , "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd",
		"Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W" , "Re", "Os", "Ir", "Pt", "Au", "Hg",
		"Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U" , "Np", "Pu", "Am", "Cm",
		"Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn",
		"Nh", "Fl", "Mc", "Lv", "Ts", "Og", 
	};
	return syms[z-1];
}

//@brief : convert to string list of elements
//@return: string representation, e.g. "H, He, Li"
std::string ElementMask::str() const {
	std::string s;
	for(size_t i = 1; i <= 118; i++) {
		if(test(i)) {
			if(!s.empty()) s += ',';
			s += Symbol(i);
		}
	}
	return s;
}

#endif//_PER_TAB_PAN_H_
