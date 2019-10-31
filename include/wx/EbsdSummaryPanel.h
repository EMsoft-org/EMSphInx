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

#ifndef _EBSD_SUM_H_
#define _EBSD_SUM_H_

#include <wx/propgrid/propgrid.h>
#include <wx/sizer.h>
#include <wx/statbox.h>
#include <wx/panel.h>

#include "wx/ValidityPanel.h"

#include "modality/ebsd/nml.hpp"

///////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
/// Class EbsdSummaryPanel
///////////////////////////////////////////////////////////////////////////////
class EbsdSummaryPanel : public ValidityPanel {
	private:

	protected:
		wxPropertyGrid* m_propGrid    ;

		//input data
		wxPGProperty* m_propInData  ;//wxPropertyCategory
		wxPGProperty* m_propIpth    ;//wxDirProperty
		wxPGProperty* m_propPatFile ;//wxFileProperty
		wxPGProperty* m_propPatDset ;//wxStringProperty
		wxPGProperty* m_propMP      ;//wxArrayStringProperty
		wxPGProperty* m_propPsym    ;//wxFileProperty

		//pattern processing
		wxPGProperty* m_patProc     ;//wxPropertyCategory
		wxPGProperty* m_propPatSz   ;//wxStringProperty
		wxPGProperty* m_propPatW    ;//wxUIntProperty
		wxPGProperty* m_propPatH    ;//wxUIntProperty
		wxPGProperty* m_propCirc    ;//wxBoolProperty
		wxPGProperty* m_propBckg    ;//wxUIntProperty
		wxPGProperty* m_propNreg    ;//wxUIntProperty

		//camera calibration
		wxPGProperty* m_propCamCalib;//wxPropertyCategory
		wxPGProperty* m_propDlt     ;//wxFloatProperty
		wxPGProperty* m_propVen     ;//wxEnumProperty
		wxPGProperty* m_propPctr    ;//wxStringProperty
		wxPGProperty* m_propPcenX   ;//wxFloatProperty
		wxPGProperty* m_propPcenY   ;//wxFloatProperty
		wxPGProperty* m_propPcenZ   ;//wxFloatProperty
		wxPGProperty* m_propThtC    ;//wxFloatProperty

		//scan dimensions
		wxPGProperty* m_propScanInfo;//wxPropertyCategory
		wxPGProperty* m_propScnSz   ;//wxStringProperty
		wxPGProperty* m_propScnW    ;//wxUIntProperty
		wxPGProperty* m_propScnH    ;//wxUIntProperty
		wxPGProperty* m_propScnDx   ;//wxFloatProperty
		wxPGProperty* m_propScnDy   ;//wxFloatProperty
		wxPGProperty* m_propRoi     ;//wxStringProperty

		//indexing parameters
		wxPGProperty* m_propIdxPrm  ;//wxPropertyCategory
		wxPGProperty* m_propBw      ;//wxUIntProperty
		wxPGProperty* m_propNrm     ;//wxBoolProperty
		wxPGProperty* m_propRef     ;//wxBoolProperty
		wxPGProperty* m_propThd     ;//wxUIntProperty
		wxPGProperty* m_propBatSz   ;//wxUIntProperty

		//output data
		wxPGProperty* m_propOutData ;//wxPropertyCategory
		wxPGProperty* m_propOpth    ;//wxDirProperty
		wxPGProperty* m_propDatafile;//wxFileProperty
		wxPGProperty* m_propVenFile ;//wxFileProperty
		wxPGProperty* m_propIpf     ;//wxFileProperty
		wxPGProperty* m_propQual    ;//wxFileProperty

		void updatePctr() {}//TODO update values on change

		void PropChanged( wxPropertyGridEvent& event ) {if(m_propVen == event.GetProperty()) updatePctr();}

	public:

		void setIdx(long bw, bool nrm, bool ref) {
			m_propBw ->SetValue(WXVARIANT(bw ));
			m_propNrm->SetValue(WXVARIANT(nrm));
			m_propRef->SetValue(WXVARIANT(ref));
		}

		void setScnDims(long w, long h, double x, double y) {
			m_propScnW ->SetValue(WXVARIANT(w));
			m_propScnH ->SetValue(WXVARIANT(h));
			m_propScnDx->SetValue(WXVARIANT(x));
			m_propScnDy->SetValue(WXVARIANT(y));
		}

		void setRoi(emsphinx::RoiSelection const & roi) {
			m_propRoi->SetValue(WXVARIANT(wxString(roi.to_string())));
		}

		void setPatDims(long w, long h) {
			m_propPatW->SetValue(WXVARIANT(w));
			m_propPatH->SetValue(WXVARIANT(h));
		}

		bool setPctr(double x, double y, double z, wxString ven) {
			if     ("EMsoft" == ven) m_propVen->SetValue(WXVARIANT(0));
			else if("Bruker" == ven) m_propVen->SetValue(WXVARIANT(1));
			else if("Edax"   == ven) m_propVen->SetValue(WXVARIANT(2));
			else if("Oxford" == ven) m_propVen->SetValue(WXVARIANT(3));
			else return false;
			m_propPcenX->SetValue(WXVARIANT(x));
			m_propPcenY->SetValue(WXVARIANT(y));
			m_propPcenZ->SetValue(WXVARIANT(z));
			return true;
		}

		void setTltDlt(double tht, double dlt) {
			m_propThtC->SetValue(WXVARIANT(tht));
			m_propDlt ->SetValue(WXVARIANT(dlt));
		}

		void setImPrc(int cir, bool bckg, int nrg) {
			m_propCirc->SetValue(WXVARIANT(cir ));
			m_propBckg->SetValue(WXVARIANT(bckg));
			m_propNreg->SetValue(WXVARIANT(nrg ));
		}

		void setInputs(std::vector<std::string> mp, wxString pf, wxString aux) {
			wxArrayString arr;
			for(const std::string& str : mp) arr.Add(str);
			m_propMP     ->SetValue(WXVARIANT(arr));
			m_propPatFile->SetValue(WXVARIANT(pf ));
			m_propPatDset->SetValue(WXVARIANT(aux));
		}

		void setOutputs(wxString df, wxString vf, wxString ipf, wxString ci) {
			m_propDatafile->SetValue(WXVARIANT(df ));
			m_propVenFile ->SetValue(WXVARIANT(vf ));
			m_propIpf     ->SetValue(WXVARIANT(ipf));
			m_propQual    ->SetValue(WXVARIANT(ci ));
		}

		void setNamelist(emsphinx::ebsd::Namelist const & nml) {
			setIdx(nml.bw, nml.normed, nml.refine);
			setScnDims(nml.scanDims[0], nml.scanDims[1], nml.scanSteps[0], nml.scanSteps[1]);
			setRoi(nml.roi);
			setPatDims(nml.patDims[0], nml.patDims[1]);
			setPctr(nml.pctr[0], nml.pctr[1], nml.pctr[2], nml.ven);
			setImPrc(nml.circRad, nml.gausBckg, nml.nRegions);
			setTltDlt(nml.thetac, nml.delta);
			setInputs(nml.masterFiles, nml.patFile, nml.patName);
			setOutputs(nml.dataFile, nml.vendorFile, nml.ipfName, nml.qualName);
		}

		bool updateNamelist(emsphinx::ebsd::Namelist* nml) const {
			bool chg = false;
			long        vLong;
			bool        vBool;
			double      vDoub;
			std::string vStr ;

			vLong = m_propBw  ->DoGetValue().GetLong   (); if(vLong != nml->bw          ) {chg = true; nml->bw           = vLong;}
			vBool = m_propNrm ->DoGetValue().GetBool   (); if(vBool != nml->normed      ) {chg = true; nml->normed       = vBool;}
			vBool = m_propRef ->DoGetValue().GetBool   (); if(vBool != nml->refine      ) {chg = true; nml->refine       = vBool;}

			vLong = m_propScnW ->DoGetValue().GetLong  (); if(vLong != nml->scanDims [0]) {chg = true; nml->scanDims [0] = vLong;}
			vLong = m_propScnH ->DoGetValue().GetLong  (); if(vLong != nml->scanDims [1]) {chg = true; nml->scanDims [1] = vLong;}
			vDoub = m_propScnDx->DoGetValue().GetDouble(); if(vDoub != nml->scanSteps[0]) {chg = true; nml->scanSteps[0] = vDoub;}
			vDoub = m_propScnDy->DoGetValue().GetDouble(); if(vDoub != nml->scanSteps[1]) {chg = true; nml->scanSteps[1] = vDoub;}

			vStr = m_propRoi->DoGetValue().GetString().ToStdString(); if(vStr != nml->roi.to_string()) {chg = true; nml->roi.from_string(vStr);}

			vLong = m_propPatW->DoGetValue().GetLong   (); if(vLong != nml->patDims  [0]) {chg = true; nml->patDims  [0] = vLong;}
			vLong = m_propPatH->DoGetValue().GetLong   (); if(vLong != nml->patDims  [1]) {chg = true; nml->patDims  [1] = vLong;}

			vDoub = m_propPcenX->DoGetValue().GetDouble(); if(vDoub != nml->pctr     [0]) {chg = true; nml->pctr     [0] = vDoub;}
			vDoub = m_propPcenY->DoGetValue().GetDouble(); if(vDoub != nml->pctr     [1]) {chg = true; nml->pctr     [1] = vDoub;}
			vDoub = m_propPcenZ->DoGetValue().GetDouble(); if(vDoub != nml->pctr     [2]) {chg = true; nml->pctr     [2] = vDoub;}
			switch(m_propVen->DoGetValue().GetLong()) {
				case 0: vStr = "EMsoft"; break;
				case 1: vStr = "Bruker"; break;
				case 2: vStr = "Edax"  ; break;
				case 3: vStr = "Oxford"; break;
			}
			if(vStr != nml->ven) {chg = true; nml->ven = vStr;}

			vLong = m_propCirc ->DoGetValue().GetLong  (); if(vLong != nml->circRad     ) {chg = true; nml->circRad      = vLong;}
			vBool = m_propBckg ->DoGetValue().GetBool  (); if(vBool != nml->gausBckg    ) {chg = true; nml->gausBckg     = vBool;}
			vLong = m_propNreg ->DoGetValue().GetLong  (); if(vLong != nml->nRegions    ) {chg = true; nml->nRegions     = vLong;}

			vDoub = m_propThtC ->DoGetValue().GetDouble(); if(vDoub != nml->thetac      ) {chg = true; nml->thetac       = vDoub;}
			vDoub = m_propDlt  ->DoGetValue().GetDouble(); if(vDoub != nml->delta       ) {chg = true; nml->delta        = vDoub;}

			wxArrayString arr = m_propMP    ->DoGetValue().GetArrayString();
			bool sameMp = arr.GetCount() == nml->masterFiles.size();
			if(sameMp) {
				for(size_t i = 0; i < arr.GetCount(); i++) {
					if(nml->masterFiles[i] != arr[i].ToStdString()) {
						sameMp = false;
						break;
					}
				}
			}
			if(!sameMp) {
				chg = true;
				nml->masterFiles.clear();
				for(size_t i = 0; i < arr.GetCount(); i++) nml->masterFiles.push_back(arr[i].ToStdString());
			}
			
			vStr  = m_propPatFile->DoGetValue().GetString().ToStdString(); if(vStr != nml->patFile) {chg = true; nml->patFile = vStr;}
			vStr  = m_propPatDset->DoGetValue().GetString().ToStdString(); if(vStr != nml->patName) {chg = true; nml->patName = vStr;}

			vStr = m_propDatafile->DoGetValue().GetString().ToStdString(); if(vStr != nml->dataFile  ) {chg = true; nml->dataFile   = vStr;}
			vStr = m_propVenFile ->DoGetValue().GetString().ToStdString(); if(vStr != nml->vendorFile) {chg = true; nml->vendorFile = vStr;}
			vStr = m_propIpf     ->DoGetValue().GetString().ToStdString(); if(vStr != nml->ipfName   ) {chg = true; nml->ipfName    = vStr;}
			vStr = m_propQual    ->DoGetValue().GetString().ToStdString(); if(vStr != nml->qualName  ) {chg = true; nml->qualName   = vStr;}

			vLong = m_propThd  ->DoGetValue().GetLong(); if(vLong !=  nml->nThread  ) {chg = true;  nml->nThread   = vLong;}
			vLong = m_propBatSz->DoGetValue().GetLong(); if(vLong !=  nml->batchSize) {chg = true;  nml->batchSize = vLong;}

			return chg;
		}

		void EnableEditing(const bool enb) {
			m_propInData  ->Enable(enb);
			m_patProc     ->Enable(enb);
			m_propCamCalib->Enable(enb);
			m_propScanInfo->Enable(enb);
			m_propIdxPrm  ->Enable(enb);
			m_propOutData ->Enable(enb);
		}

		void EnableWizardEditing(const bool enb) {
			m_propInData  ->Enable(enb);
			m_propPatSz   ->Enable(enb);
			m_propVen     ->Enable(enb);
			m_propPctr    ->Enable(enb);
			m_propScanInfo->Enable(enb);
		}

		//@brief: get the current master pattern list
		wxArrayString getMP() const {return m_propMP->DoGetValue().GetArrayString();}

		EbsdSummaryPanel( wxWindow* parent, wxWindowID id = wxID_ANY, const wxPoint& pos = wxDefaultPosition, const wxSize& size = wxSize( 569,588 ), long style = wxTAB_TRAVERSAL, const wxString& name = wxEmptyString );
		~EbsdSummaryPanel();
};

///////////////////////////////////////////////////////////////////////////
// C++ code generated with wxFormBuilder (version Oct 26 2018)
// http://www.wxformbuilder.org/
//
// PLEASE DO *NOT* EDIT THIS FILE!
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////

EbsdSummaryPanel::EbsdSummaryPanel( wxWindow* parent, wxWindowID id, const wxPoint& pos, const wxSize& size, long style, const wxString& name ) : ValidityPanel( parent, id, pos, size, style, name ) {
	wxStaticBoxSizer* sbSum = new wxStaticBoxSizer( new wxStaticBox( this, wxID_ANY, wxT("Summary") ), wxVERTICAL );

	wxPGChoices chs;
	chs.Add("EMsoft", 0);
	chs.Add("Bruker", 1);
	chs.Add("Edax"  , 2);
	chs.Add("Oxford", 3);
	wxString cmp("<composed>");

	m_propGrid = new wxPropertyGrid(sbSum->GetStaticBox(), wxID_ANY, wxDefaultPosition, wxDefaultSize, wxPG_DEFAULT_STYLE|wxPG_SPLITTER_AUTO_CENTER|wxTAB_TRAVERSAL);

	//build input data
	m_propInData   = m_propGrid->Append  (                new wxPropertyCategory   ( wxT("Input Files"        ), wxPG_LABEL      ) );
	// m_propIpth     = m_propGrid->AppendIn(m_propInData  , new wxDirProperty        ( wxT("ipath"              ), wxPG_LABEL      ) );
	m_propPatFile  = m_propGrid->AppendIn(m_propInData  , new wxFileProperty       ( wxT("patfile"            ), wxPG_LABEL      ) );
	m_propPatDset  = m_propGrid->AppendIn(m_propInData  , new wxStringProperty     ( wxT("patdset"            ), wxPG_LABEL      ) );
	m_propMP       = m_propGrid->AppendIn(m_propInData  , new wxArrayStringProperty( wxT("masterfile"         ), wxPG_LABEL      ) );
	// m_propPsym     = m_propGrid->AppendIn(m_propInData  , new wxFileProperty       ( wxT("psymfile"           ), wxPG_LABEL      ) );

	//build pattern processing
	m_patProc      = m_propGrid->Append  (                new wxPropertyCategory   ( wxT("Pattern Processing" ), wxPG_LABEL      ) );
	m_propPatSz    = m_propGrid->AppendIn(m_patProc     , new wxStringProperty     ( wxT("patdims"            ), wxPG_LABEL, cmp ) );
	m_propPatW     = m_propGrid->AppendIn(m_propPatSz   , new wxUIntProperty       ( wxT("w"                  ), wxT("pat_w"   ) ) );
	m_propPatH     = m_propGrid->AppendIn(m_propPatSz   , new wxUIntProperty       ( wxT("h"                  ), wxT("pat_h"   ) ) );
	m_propCirc     = m_propGrid->AppendIn(m_patProc     , new wxIntProperty        ( wxT("circmask"           ), wxPG_LABEL      ) );
	m_propBckg     = m_propGrid->AppendIn(m_patProc     , new wxBoolProperty       ( wxT("gausbckg"           ), wxPG_LABEL      ) );
	m_propNreg     = m_propGrid->AppendIn(m_patProc     , new wxUIntProperty       ( wxT("nregions"           ), wxPG_LABEL      ) );

	//build camera calibration
	m_propCamCalib = m_propGrid->Append  (                new wxPropertyCategory   ( wxT("Camera Calibration" ), wxPG_LABEL      ) );
	m_propDlt      = m_propGrid->AppendIn(m_propCamCalib, new wxFloatProperty      ( wxT("delta"              ), wxPG_LABEL      ) );
	m_propVen      = m_propGrid->AppendIn(m_propCamCalib, new wxEnumProperty       ( wxT("vendor"             ), wxPG_LABEL, chs ) );
	m_propPctr     = m_propGrid->AppendIn(m_propCamCalib, new wxStringProperty     ( wxT("pctr"               ), wxPG_LABEL, cmp ) );
	m_propPcenX    = m_propGrid->AppendIn(m_propPctr    , new wxFloatProperty      ( wxT("x"                  ), wxT("pctr_x"  ) ) );
	m_propPcenY    = m_propGrid->AppendIn(m_propPctr    , new wxFloatProperty      ( wxT("y"                  ), wxT("pctr_y"  ) ) );
	m_propPcenZ    = m_propGrid->AppendIn(m_propPctr    , new wxFloatProperty      ( wxT("z"                  ), wxT("pctr_z"  ) ) );
	m_propThtC     = m_propGrid->AppendIn(m_propCamCalib, new wxFloatProperty      ( wxT("thetac"             ), wxPG_LABEL      ) );
	
	//build scan dimensions
	m_propScanInfo = m_propGrid->Append  (                new wxPropertyCategory   ( wxT("Scan Information"   ), wxPG_LABEL      ) );
	m_propScnSz    = m_propGrid->AppendIn(m_propScanInfo, new wxStringProperty     ( wxT("scandims"           ), wxPG_LABEL, cmp ) );
	m_propScnW     = m_propGrid->AppendIn(m_propScnSz   , new wxUIntProperty       ( wxT("w"                  ), wxT("scn_w"   ) ) );
	m_propScnH     = m_propGrid->AppendIn(m_propScnSz   , new wxUIntProperty       ( wxT("h"                  ), wxT("scn_h"   ) ) );
	m_propScnDx    = m_propGrid->AppendIn(m_propScnSz   , new wxFloatProperty      ( wxT("dx"                 ), wxT("scn_dx"  ) ) );
	m_propScnDy    = m_propGrid->AppendIn(m_propScnSz   , new wxFloatProperty      ( wxT("dy"                 ), wxT("scn_dy"  ) ) );
	m_propRoi      = m_propGrid->AppendIn(m_propScanInfo, new wxStringProperty     ( wxT("roimask"            ), wxPG_LABEL      ) );

	//build indexing parameters
	m_propIdxPrm   = m_propGrid->Append  (                new wxPropertyCategory   ( wxT("Indexing Parameters"), wxPG_LABEL      ) );
	m_propBw       = m_propGrid->AppendIn(m_propIdxPrm  , new wxUIntProperty       ( wxT("bw"                 ), wxPG_LABEL      ) );
	m_propNrm      = m_propGrid->AppendIn(m_propIdxPrm  , new wxBoolProperty       ( wxT("normed"             ), wxPG_LABEL      ) );
	m_propRef      = m_propGrid->AppendIn(m_propIdxPrm  , new wxBoolProperty       ( wxT("refine"             ), wxPG_LABEL      ) );
	m_propThd      = m_propGrid->AppendIn(m_propIdxPrm  , new wxUIntProperty       ( wxT("nthread"            ), wxPG_LABEL      ) );
	m_propBatSz    = m_propGrid->AppendIn(m_propIdxPrm  , new wxUIntProperty       ( wxT("batchsize"          ), wxPG_LABEL      ) );

	//build output data
	m_propOutData  = m_propGrid->Append  (                new wxPropertyCategory   ( wxT("Output Files"       ), wxPG_LABEL      ) );
	// m_propOpth     = m_propGrid->AppendIn(m_propOutData , new wxDirProperty        ( wxT("opath"              ), wxPG_LABEL      ) );
	m_propDatafile = m_propGrid->AppendIn(m_propOutData , new wxFileProperty       ( wxT("datafile"           ), wxPG_LABEL      ) );
	m_propVenFile  = m_propGrid->AppendIn(m_propOutData , new wxFileProperty       ( wxT("vendorfile"         ), wxPG_LABEL      ) );
	m_propIpf      = m_propGrid->AppendIn(m_propOutData , new wxFileProperty       ( wxT("ipfmap"             ), wxPG_LABEL      ) );
	m_propQual     = m_propGrid->AppendIn(m_propOutData , new wxFileProperty       ( wxT("qualmap"            ), wxPG_LABEL      ) );
	
	//formatting
	m_propGrid->SetPropertyAttributeAll(wxPG_BOOL_USE_CHECKBOX, true);//use checkboxes instead of true/false dropdowns

	//final assmebly
	sbSum->Add( m_propGrid, 1, wxALL|wxEXPAND, 5 );
	this->SetSizer( sbSum );
	this->Layout();
	
	//connect events
	m_propGrid->Connect( wxEVT_PG_CHANGED, wxPropertyGridEventHandler( EbsdSummaryPanel::PropChanged ), NULL, this );
}

EbsdSummaryPanel::~EbsdSummaryPanel() {
	m_propGrid->Disconnect( wxEVT_PG_CHANGED, wxPropertyGridEventHandler( EbsdSummaryPanel::PropChanged ), NULL, this );
}

#endif//_EBSD_SUM_H_
