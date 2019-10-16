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

#ifndef _EBSD_WIZARD_H
#define _EBSD_WIZARD_H

#include <wx/artprov.h>
#include <wx/xrc/xmlres.h>
#include <wx/gdicmn.h>
#include <wx/font.h>
#include <wx/colour.h>
#include <wx/settings.h>
#include <wx/string.h>
#include <wx/sizer.h>
#include <wx/bitmap.h>
#include <wx/image.h>
#include <wx/icon.h>
#include <wx/wizard.h>
#include <wx/dynarray.h>
#include <wx/fileconf.h>
#include <wx/msgdlg.h>
#include <wx/filedlg.h>

#include "PatternLoadPanel.h"
#include "MasterPatternSelectPanel.h"
#include "PatternCenterPanel.h"
#include "ScanDimsPanel.h"
#include "IdxParamPanel.h"
#include "EbsdSummaryPanel.h"

#include "ValidityWizard.h"

#include "BibtexDialog.h"
#include "IndexingFrame.h"

///////////////////////////////////////////////////////////////////////////

#include "util/sysnames.hpp" // could also use wxStandardPaths

#include "xtal/orientation_map.hpp"
#include "modality/ebsd/pattern.hpp"
#include "modality/ebsd/nml.hpp"

///////////////////////////////////////////////////////////////////////////////
/// Class EbsdNamelistWizard
///////////////////////////////////////////////////////////////////////////////
class EbsdNamelistWizard : public ValidityWizard
{
	private:
		emsphinx::ebsd::Namelist nml;

	protected:
		PatternLoadPanel        * m_patLoadPan   ;
		MasterPatternSelectPanel* m_mastPatSelPan;
		PatternCenterPanel      * m_patCenPan    ;
		ScanDimsPanel           * m_scnDimPan    ;
		IdxParamPanel           * m_idxPrmPan    ;
		EbsdSummaryPanel        * m_ebsdSumPan   ;

		virtual void OnNext( wxCommandEvent& event ) {
			if(5 == m_book->GetSelection()) {
				wxMessageDialog dlg(this, "Namelist Complete", "", wxYES_NO|wxCANCEL|wxCENTRE|wxICON_QUESTION);
				dlg.SetYesNoLabels ("Index Now", "Export Namelist");
				int res = dlg.ShowModal();

				if(wxID_NO == res) {//save namelist file
					wxFileDialog svDlg(this, _("Save Namelist file"), "", "", "Namelist files (*.nml)|*.nml", wxFD_SAVE|wxFD_OVERWRITE_PROMPT);
					if(svDlg.ShowModal() == wxID_CANCEL) return;//cancel
					std::ofstream os(svDlg.GetPath().ToStdString());
					os << nml.to_string();
				} else if(wxID_YES == res) {
					//create indexing window
					// Close();
					IndexingFrame* frm = new IndexingFrame(NULL);
					frm->setNml(nml);
					wxImage wim = m_scnDimPan->getMap();
					frm->setImage(wim);
					frm->Show();
					// return;
				} else {//don't do anything on wxID_CANCEL
					return;
				}

				//if we made it this far everything made it, ask users to cite papers
				if(!wxFileExists(getUserAppDataDir() + "silenced")) {
					BibtexDialog bDlg(this);
					bDlg.ShowModal();
					if(bDlg.Silence()) std::ofstream os(getUserAppDataDir() + "silenced");
				}
			}
			ValidityWizard::OnNext(event);
		}

		virtual void PageChanging( wxBookCtrlEvent& event ) {
			switch(event.GetOldSelection()) {//what page are we leaving
				case 0://pattern load
					if(!m_patLoadPan->LoadImages()) event.Veto();//make sure patterns are valid before allowing page turn
					nml.patFile = m_patLoadPan->getFile().ToStdString();
					nml.patName = m_patLoadPan->getAux();
					nml.patDims[0] = m_patLoadPan->getW();
					nml.patDims[1] = m_patLoadPan->getH();

					nml.circRad  = m_patLoadPan->getCirc();
					nml.gausBckg = m_patLoadPan->getBckg();
					nml.nRegions = m_patLoadPan->getNreg();

					m_patCenPan->setBinnedPix(m_patLoadPan->getW(), m_patLoadPan->getH());//save binned detector size
				break;

				case 1://master pattern
					nml.masterFiles = m_mastPatSelPan->getSelected();
				break;

				case 2://pattern center
					m_patCenPan->getPatternCenter(nml.pctr[0], nml.pctr[1], nml.pctr[2]);
					nml.ven = "EMsoft";
					nml.thetac = m_patCenPan->getDetTlt();
					nml.delta = m_patCenPan->getDelta();
				break;

				case 3://scan dimensions
					nml.scanDims[0] = m_scnDimPan->getW();
					nml.scanDims[1] = m_scnDimPan->getH();
					nml.scanSteps[0] = m_scnDimPan->getX();
					nml.scanSteps[1] = m_scnDimPan->getY();
					nml.roi = m_scnDimPan->getRoi();
				break;

				case 4://indexing paramters
					nml.bw     = m_idxPrmPan->getBw  ();
					nml.normed = m_idxPrmPan->getNorm();
					nml.refine = m_idxPrmPan->getRef ();
					nml.dataFile   = m_idxPrmPan->getDataFile  ();
					nml.vendorFile = m_idxPrmPan->getVendorFile();
					nml.ipfName    = m_idxPrmPan->getIpfFile   ();
					nml.qualName   = m_idxPrmPan->getCiFile    ();
				break;

				case 5://ebsd summary
				break;

			}

			if(5 == event.GetSelection()) {
				m_ebsdSumPan->setNamelist(&nml);
			}
		}

		void PatternFileChanged(wxFileDirPickerEvent& event) {
			//search for associated scan file and use to populate subsequent pages
			nml.patFile = m_patLoadPan->getFile().ToStdString();
			nml.patName = m_patLoadPan->getAux();
			std::shared_ptr< std::vector<float> > ciMap = std::make_shared< std::vector<float> >();
			std::shared_ptr< std::vector<float> > iqMap = std::make_shared< std::vector<float> >();
			if(nml.findScanFile(iqMap.get(), ciMap.get())) {
				m_patCenPan->setPatternCenter(nml.pctr[0], nml.pctr[1], nml.pctr[2], nml.ven);
				m_patCenPan->setDetTlt(nml.thetac);
				m_scnDimPan->setW(nml.scanDims[0]);
				m_scnDimPan->setH(nml.scanDims[1]);
				m_scnDimPan->setX(nml.scanSteps[0]);
				m_scnDimPan->setY(nml.scanSteps[1]);
			} else {
				m_patCenPan->clear();//clear any stored pattern center
			}

			m_scnDimPan->setMaps(iqMap, ciMap, m_patLoadPan->getNum());//set maps (empty or not) + pattern count
			m_scnDimPan->setPats( emsphinx::ebsd::PatternFile::Read(m_patLoadPan->getFile().ToStdString(), m_patLoadPan->getAux(), m_patLoadPan->getW(), m_patLoadPan->getH()) );
		}

	public:

		EbsdNamelistWizard( wxWindow* parent, wxWindowID id = wxID_ANY, const wxString& title = wxT("EBSD Indexing Wizard"), const wxBitmap& bitmap = wxNullBitmap, const wxPoint& pos = wxDefaultPosition, long style = wxDEFAULT_DIALOG_STYLE );
		~EbsdNamelistWizard();
};

///////////////////////////////////////////////////////////////////////////
// C++ code generated with wxFormBuilder (version Oct 26 2018)
// http://www.wxformbuilder.org/
//
// PLEASE DO *NOT* EDIT THIS FILE!
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////

EbsdNamelistWizard::EbsdNamelistWizard( wxWindow* parent, wxWindowID id, const wxString& title, const wxBitmap& bitmap, const wxPoint& pos, long style ) : 
ValidityWizard( parent, id, title, pos )
{
	m_patLoadPan    = new PatternLoadPanel        (m_book);
	m_mastPatSelPan = new MasterPatternSelectPanel(m_book);
	m_patCenPan     = new PatternCenterPanel      (m_book);
	m_scnDimPan     = new ScanDimsPanel           (m_book);
	m_idxPrmPan     = new IdxParamPanel           (m_book);
	m_ebsdSumPan    = new EbsdSummaryPanel        (m_book);
	m_ebsdSumPan->EnableEditing(false);//editin shoulw be done through wizard

	AddPage(m_patLoadPan   );
	AddPage(m_mastPatSelPan);
	AddPage(m_patCenPan    );
	AddPage(m_scnDimPan    );
	AddPage(m_idxPrmPan    );
	AddPage(m_ebsdSumPan   );
	m_book->SetSelection(0);
	this->SetClientSize(m_book->GetBestSize());

	//read known master patterns from file
	std::vector<wxString> masterLib;
	wxFileConfig config( wxEmptyString, wxEmptyString, getUserAppDataDir() + "MasterLib.ini");
	size_t numEtr = config.GetNumberOfEntries();
	for(size_t i = 0; i < numEtr; i++) {
		std::ostringstream ss;
		ss << i;
		wxString str;
		config.Read(ss.str(), &str);
		masterLib.push_back(str);
	}
	m_mastPatSelPan->setLibrary(masterLib);

	m_patLoadPan->Connect( wxEVT_COMMAND_FILEPICKER_CHANGED, wxFileDirPickerEventHandler( EbsdNamelistWizard::PatternFileChanged ), NULL, this );
}

EbsdNamelistWizard::~EbsdNamelistWizard() {
	m_patLoadPan->Disconnect( wxEVT_COMMAND_FILEPICKER_CHANGED, wxFileDirPickerEventHandler( EbsdNamelistWizard::PatternFileChanged ), NULL, this );

	//save known master patterns to config file
	std::vector<wxString> masterLib = m_mastPatSelPan->getLibrary();
	wxFileConfig config( wxEmptyString, wxEmptyString, getUserAppDataDir() + "MasterLib.ini");
	config.DeleteAll();
	for(size_t i = 0; i < masterLib.size(); i++) {
		std::ostringstream ss;
		ss << i;
		config.Write(ss.str(), masterLib[i]);
	}
	config.Flush();
}

#endif//_EBSD_WIZARD_H