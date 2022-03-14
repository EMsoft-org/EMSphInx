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
#include <wx/valnum.h>

///////////////////////////////////////////////////////////////////////////////
/// Class MpConvertDialog
///////////////////////////////////////////////////////////////////////////////
class MpConvertDialog : public wxDialog {
	private:
		wxIntegerValidator<int> valBw;//bandwidth validator

	protected:
		wxFilePickerCtrl* m_fileIn ;
		wxTextCtrl      * m_txtForm;
		wxTextCtrl      * m_txtName;
		wxTextCtrl      * m_txtSSyb;
		wxTextCtrl      * m_txtBw  ;
		wxFilePickerCtrl* m_fileOut;
		wxButton        * m_button ;

		// Virtual event handlers, overide them in your derived class
		virtual void OnInputChanged  ( wxFileDirPickerEvent& event );
		virtual void OnFormulaChanged( wxCommandEvent      & event ) {m_button->Enable(isValid());}
		virtual void OnBwChanged     ( wxCommandEvent      & event ) {m_button->Enable(isValid());}
		virtual void OnOutputChanged ( wxFileDirPickerEvent& event ) {m_button->Enable(isValid());}
		virtual void OnConvert       ( wxCommandEvent      & event ) {EndModal(wxID_OK);}

		bool isValid() const {
			return wxFileExists(m_fileIn->GetPath()) &&//we need an input file
			       !m_txtForm->GetValue().IsEmpty()  &&//we need a non-empty formula
			       !m_txtBw  ->GetValue().IsEmpty()  &&//we need a bandwidth
			       !m_fileOut->GetPath ().IsEmpty()   ;//we need an output file
		}

	public:

		std::string getInput  () const {return m_fileIn ->GetPath ().ToStdString();}
		std::string getFormula() const {return m_txtForm->GetValue().ToStdString();}
		std::string getName   () const {return m_txtName->GetValue().ToStdString();}
		std::string getSSyb   () const {return m_txtSSyb->GetValue().ToStdString();}
		size_t      getBw     () const {long v; m_txtBw->GetLineText(0).ToLong(&v); return (size_t)v;}
		std::string getOutput () const {return m_fileOut->GetPath ().ToStdString();}

		MpConvertDialog( wxWindow* parent, wxWindowID id = wxID_ANY, const wxString& title = wxT("Convert EMsoft Master Pattern"), const wxPoint& pos = wxDefaultPosition, const wxSize& size = wxSize( 300,350 ), long style = wxDEFAULT_DIALOG_STYLE );
		~MpConvertDialog();

};

#include <wx/filename.h>
#include <wx/msgdlg.h>
#include <wx/progdlg.h>
#include "idx/master.hpp"

MpConvertDialog::MpConvertDialog( wxWindow* parent, wxWindowID id, const wxString& title, const wxPoint& pos, const wxSize& size, long style ) : wxDialog( parent, id, title, pos, size, style ) {
	this->SetSizeHints( wxDefaultSize, wxDefaultSize );

	valBw.SetRange(64, 512);
	/*
	Unwanted behavior here! If this is active and the text entry field is blank, numbers cannot be typed in.
	This is because any single-digit value (0-9) is below the minimum value allowed and will therefore be rejected.
	While typing in a number will not work, copying and pasting a value will work - including for values outside of the allowed range.
	Not sure how this reacts to the above, so I simply removed it from line 126 instead. I still don't recommend going below 64 or above ~1000.
	*/

	wxFlexGridSizer* fgSizer = new wxFlexGridSizer( 7, 2, 0, 0 );
	fgSizer->AddGrowableCol( 1 );
	fgSizer->AddGrowableRow( 0 ); fgSizer->AddGrowableRow( 1 );
	fgSizer->AddGrowableRow( 2 ); fgSizer->AddGrowableRow( 3 );
	fgSizer->AddGrowableRow( 4 ); fgSizer->AddGrowableRow( 5 );
	fgSizer->AddGrowableRow( 6 );
	fgSizer->SetFlexibleDirection( wxHORIZONTAL );
	fgSizer->SetNonFlexibleGrowMode( wxFLEX_GROWMODE_ALL );

	wxStaticText* staTxtIn   = new wxStaticText( this, wxID_ANY, wxT("Input File" ), wxDefaultPosition, wxDefaultSize, 0 ); staTxtIn  ->Wrap( -1 );
	wxStaticText* staTxtForm = new wxStaticText( this, wxID_ANY, wxT("Formula"    ), wxDefaultPosition, wxDefaultSize, 0 ); staTxtForm->Wrap( -1 );
	wxStaticText* staTxtName = new wxStaticText( this, wxID_ANY, wxT("Name"       ), wxDefaultPosition, wxDefaultSize, 0 ); staTxtName->Wrap( -1 );
	wxStaticText* staTxtSSyb = new wxStaticText( this, wxID_ANY, wxT("S. Syb."    ), wxDefaultPosition, wxDefaultSize, 0 ); staTxtSSyb->Wrap( -1 );
	wxStaticText* staTxtBw   = new wxStaticText( this, wxID_ANY, wxT("Bandwidth"  ), wxDefaultPosition, wxDefaultSize, 0 ); staTxtSSyb->Wrap( -1 );
	wxStaticText* staTxtOut  = new wxStaticText( this, wxID_ANY, wxT("Output File"), wxDefaultPosition, wxDefaultSize, 0 ); staTxtOut ->Wrap( -1 );

	m_fileIn  = new wxFilePickerCtrl( this, wxID_ANY, wxEmptyString , wxT("Input Master Pattern" ), wxT("*.h5"), wxDefaultPosition, wxDefaultSize, wxFLP_DEFAULT_STYLE|wxFLP_SMALL );
	m_txtForm = new wxTextCtrl      ( this, wxID_ANY, wxT("Unknown"), wxDefaultPosition, wxDefaultSize, 0        );
	m_txtName = new wxTextCtrl      ( this, wxID_ANY, wxEmptyString , wxDefaultPosition, wxDefaultSize, 0        );
	m_txtSSyb = new wxTextCtrl      ( this, wxID_ANY, wxEmptyString , wxDefaultPosition, wxDefaultSize, 0        );
	m_txtBw   = new wxTextCtrl      ( this, wxID_ANY, wxT("384"    ), wxDefaultPosition, wxDefaultSize, 0        );
	m_fileOut = new wxFilePickerCtrl( this, wxID_ANY, wxEmptyString , wxT("Output Master Pattern"), wxT("*.sht"), wxDefaultPosition, wxDefaultSize, wxFLP_OVERWRITE_PROMPT|wxFLP_SAVE|wxFLP_SMALL|wxFLP_USE_TEXTCTRL );

	m_button = new wxButton( this, wxID_ANY, wxT("Convert"), wxDefaultPosition, wxDefaultSize, 0 );

	m_button->SetDefault();
	m_button->Enable(false);

	fgSizer->Add( staTxtIn  , 0, wxALL|wxALIGN_CENTER_VERTICAL              , 5 );
	fgSizer->Add( m_fileIn  , 0, wxALL|wxALIGN_CENTER_VERTICAL|wxEXPAND     , 5 );
	fgSizer->Add( staTxtForm, 0, wxALL|wxALIGN_CENTER_VERTICAL              , 5 );
	fgSizer->Add( m_txtForm , 0, wxALL|wxALIGN_CENTER_VERTICAL|wxEXPAND     , 5 );
	fgSizer->Add( staTxtName, 0, wxALL|wxALIGN_CENTER_VERTICAL              , 5 );
	fgSizer->Add( m_txtName , 0, wxALL|wxALIGN_CENTER_VERTICAL|wxEXPAND     , 5 );
	fgSizer->Add( staTxtSSyb, 0, wxALL|wxALIGN_CENTER_VERTICAL              , 5 );
	fgSizer->Add( m_txtSSyb , 0, wxALL|wxALIGN_CENTER_VERTICAL|wxEXPAND     , 5 );
	fgSizer->Add( staTxtBw  , 0, wxALL|wxALIGN_CENTER_VERTICAL              , 5 );
	fgSizer->Add( m_txtBw   , 0, wxALL|wxALIGN_CENTER_VERTICAL|wxEXPAND     , 5 );
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
	m_txtBw  ->Connect( wxEVT_COMMAND_TEXT_UPDATED      , wxCommandEventHandler      ( MpConvertDialog::OnBwChanged      ), NULL, this );
	m_fileOut->Connect( wxEVT_COMMAND_FILEPICKER_CHANGED, wxFileDirPickerEventHandler( MpConvertDialog::OnOutputChanged  ), NULL, this );
	m_button ->Connect( wxEVT_COMMAND_BUTTON_CLICKED    , wxCommandEventHandler      ( MpConvertDialog::OnConvert        ), NULL, this );
}

MpConvertDialog::~MpConvertDialog() {
	// Disconnect Events
	m_fileIn ->Disconnect( wxEVT_COMMAND_FILEPICKER_CHANGED, wxFileDirPickerEventHandler( MpConvertDialog::OnInputChanged   ), NULL, this );
	m_txtForm->Disconnect( wxEVT_COMMAND_TEXT_UPDATED      , wxCommandEventHandler      ( MpConvertDialog::OnFormulaChanged ), NULL, this );
	m_txtBw  ->Disconnect( wxEVT_COMMAND_TEXT_UPDATED      , wxCommandEventHandler      ( MpConvertDialog::OnBwChanged      ), NULL, this );
	m_fileOut->Disconnect( wxEVT_COMMAND_FILEPICKER_CHANGED, wxFileDirPickerEventHandler( MpConvertDialog::OnOutputChanged  ), NULL, this );
	m_button ->Disconnect( wxEVT_COMMAND_BUTTON_CLICKED    , wxCommandEventHandler      ( MpConvertDialog::OnConvert        ), NULL, this );
}

#include <array>
#include <cmath>
#include <sstream>
#include <map>
#include "xtal/hm.hpp"

namespace gcd_detail {
	//@brief: compute the greatest common denominator of 2 numbers
	size_t gcd(size_t a, size_t b) {return b == 0 ? a : gcd(b, a % b);}
}

void MpConvertDialog::OnInputChanged  ( wxFileDirPickerEvent& event ) {
	//check if file exists
	m_button->Enable(false);
	if(!wxFileExists(event.GetPath())) return; 

	//read in crystal data
	int32_t sgNum, sgSet, nAtom;
	std::vector<int32_t> aTy;
	std::vector<float> aCd;
	try {
		H5::Group xtalData = H5::H5File(event.GetPath().ToStdString().c_str(), H5F_ACC_RDONLY).openGroup("CrystalData");//read only access
		xtalData.openDataSet("SpaceGroupNumber" ).read(&sgNum, H5::PredType::NATIVE_INT32);
		xtalData.openDataSet("SpaceGroupSetting").read(&sgSet, H5::PredType::NATIVE_INT32);
		xtalData.openDataSet("Natomtypes"       ).read(&nAtom, H5::PredType::NATIVE_INT32);
		aTy.resize(nAtom);
		aCd.resize(nAtom * 5);
		xtalData.openDataSet("AtomData" ).read(aCd.data(), H5::PredType::NATIVE_FLOAT);
		xtalData.openDataSet("Atomtypes").read(aTy.data(), H5::PredType::NATIVE_INT32);
	} catch (std::exception& e) {
		wxMessageDialog msgDlg(this, e.what(), "Error Reading Master Pattern");
		msgDlg.ShowModal();
		m_fileIn->SetPath("");
		return;
	} catch (H5::Exception& e) {
		std::ostringstream ss;
		ss << "H5 error attempting to read EMsoft master pattern:\n";
		ss << "\tfile - " << event.GetPath().ToStdString() << '\n';
		ss << "\tfunc - " << e.getCFuncName() << '\n';
		ss << "detailed message:\n" << e.getCDetailMsg();
		wxMessageDialog msgDlg(this, ss.str(), "Error Reading Master Pattern");
		msgDlg.ShowModal();
		m_fileIn->SetPath("");
		return;
	}

	//get the multiplicity of each site
	xtal::HermannMaguin hm;
	hm.fromNumber(sgNum, 2 == sgSet);
	std::vector<xtal::GenPos> symMat = xtal::GenPos::CloseSet(hm.generators());
	struct AtomPos {
		double xyz[3];//coordinates

		//@brief: check if 2 atomic coordinates are ~equal
		//@note : this logic assumes that {x,y,z} are always in [0,1]
		bool operator==(const AtomPos& rhs) const {
			static const double eps = std::sqrt(std::numeric_limits<double>::epsilon());//semi arbitrary cutoff
			for(size_t i = 0; i < 3; i++) {
				double d = std::fabs(xyz[i] - rhs.xyz[i]);
				if(d > 0.5) d = 1.0 - d;
				if(d > eps) return false;//~arbitrary cutoff of 0.1% difference in normalized position
			}
			return true;
		}

		//@brief: arbitrary sort function for set construction
		bool operator<(const AtomPos& rhs) const {
			if(operator==(rhs)) return false;
			return std::lexicographical_compare(xyz, xyz+3, rhs.xyz, rhs.xyz+3);
		}
	};

	std::vector< std::set<AtomPos> > aPos(nAtom);//a set of symmetrically equivalent atomic positions for each atom
	for(int32_t i = 0; i < nAtom; i++) {
		AtomPos p0;
		p0.xyz[0] = aCd[nAtom*0+i];
		p0.xyz[1] = aCd[nAtom*1+i];
		p0.xyz[2] = aCd[nAtom*2+i];
		aPos[i].insert(p0);
		for(const xtal::GenPos& g : symMat) {//loop over symmetry operators
			//get general position matrix elements
			double trs[3];
			int8_t const * mat3 = g.getMat3();
			g.getTransReal(trs);

			//compute symmetrically equivalent point
			AtomPos pSym;
			for(size_t j = 0; j < 3; j++) {
				pSym.xyz[j] = p0.xyz[0] * mat3[3*j+0]
				            + p0.xyz[1] * mat3[3*j+1]
				            + p0.xyz[2] * mat3[3*j+2]
				            + trs[j];
            	pSym.xyz[j] = std::fmod(pSym.xyz[j], 1);
            	if(pSym.xyz[j] == -0.0) pSym.xyz[j] = 0;
            	else if(pSym.xyz[j] < 0) pSym.xyz[j] += 1;
			}
			aPos[i].insert(pSym);
		}
	}

	//combine different atoms on the same site
	struct AtSite {
		size_t                  mult;//multiplicity of site
		std::map<size_t, float> occ ;//accumulated occupancy of site
	};
	std::vector<AtSite> sites;
	std::vector<bool> used(nAtom, false);
	for(int32_t i = 0; i < nAtom; i++) {
		if(used[i]) continue;//we already used this site
		AtSite s;
		s.occ[aTy[i]] = aCd[nAtom * 3 + i];//save first atom at this site
		s.mult = aPos[i].size();//save multiplicity
		for(int32_t j = i+1; j < nAtom; j++) {//check if there are other atoms at the same site
			if(used[j]) continue;//we already used this site
			if(aPos[i] == aPos[j]) {//same site coordinates
				s.occ[aTy[j]] = aCd[nAtom * 3 + j];//save additional atom at same site
				used[j] = true;//flag as used
			}
		}
		sites.push_back(s);
	}

	//finally combine site with like composition
	std::vector<AtSite> cmbSite;
	std::fill(used.begin(), used.end(), false);
	for(size_t i = 0; i < sites.size(); i++) {
		if(used[i]) continue;//we already used this site
		AtSite s = sites[i];//get site
		for(size_t j = i+1; j < sites.size(); j++) {
			if(used[j]) continue;//we already used this site
			if(s.occ == sites[j].occ) {//same site chemistry
				s.mult += sites[j].mult;//combine sites
				used[j] = true;//flag as used
			}
		}
		cmbSite.push_back(s);//save this site
	}

	//get greatest common denominator of multiplicities
	size_t gcd = cmbSite.front().mult;
	for(size_t i = 1; i < cmbSite.size(); i++) gcd = gcd_detail::gcd(cmbSite[i].mult, gcd);

	//now that we have the combined sites assemble the formula
	static const std::vector<std::string> AtSyb = {
		"H" , "He", "Li", "Be", "B" , "C" , "N" , "O" , "F" , "Ne", "Na", "Mg", "Al", "Si", "P" , "S" ,
		"Cl", "Ar", "K" , "Ca", "Sc", "Ti", "V" , "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge",
		"As", "Se", "Br", "Kr", "Rb", "Sr", "Y" , "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
		"In", "Sn", "Sb", "Te", "I" , "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd",
		"Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W" , "Re", "Os", "Ir", "Pt", "Au", "Hg",
		"Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U" , "Np", "Pu", "Am", "Cm",
		"Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn",
		"Nh", "Fl", "Mc", "Lv", "Ts", "Og", 
	};
	std::ostringstream ss;
	std::vector<bool> siteUsed(nAtom, false);
	for(const AtSite& s : cmbSite) {
		const bool mult = s.occ.size() > 1;//does this site contain more than 1 atom
		const size_t num = s.mult / gcd;//normalized site multiplicity
		if(mult) {
			ss << "(";
			for(const std::pair<size_t, float> p : s.occ) {
				ss << AtSyb[p.first-1];//print atomic symbol for element
				ss << p.second;
			}
			ss << ")";
			if(num > 1) ss << num;
		} else {
			ss << AtSyb[s.occ.cbegin()->first-1];//print atomic symbol for element
			const float f = s.occ.cbegin()->second;//single species in site, get occupancy
			if(1.0f == f) {//fully occupied
				if(num > 1) ss << num;
			} else {//partial occupancy
				ss << f * num;
			}
		}
	}

	//if we made it this far we were able to parse a formula, fill in, set default bandwidth, and set output file
	wxFileName fn = event.GetPath();
	fn.SetExt("sht");
	m_txtForm->ChangeValue(ss.str());//Overrides any manual formula input with formula constructed from EMsoft master pattern data file.
	//m_txtBw->ChangeValue("384");//Overrides any manual conversion bandwidth input with fixed value. Commented out to allow GUI selection of conversion BW.
	m_fileOut->SetFileName(fn);
	m_button->Enable(true);
}

#endif//_MP_CVT_DLG_H_
