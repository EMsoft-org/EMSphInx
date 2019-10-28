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

#ifndef _IDX_FRAME_H
#define _IDX_FRAME_H

#include <wx/sizer.h>
#include <wx/gauge.h>
#include <wx/button.h>
#include <wx/statusbr.h>
#include <wx/frame.h>
#include <wx/splitter.h>
#include <wx/progdlg.h>
#include <wx/msgdlg.h>

#include "ImagePanel.h"
#include "EbsdSummaryPanel.h"

#include "modality/ebsd/nml.hpp"
#include "modality/ebsd/idx.hpp"
#include "util/timer.hpp"

///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
/// Class IndexingFrame
///////////////////////////////////////////////////////////////////////////////

BEGIN_DECLARE_EVENT_TYPES()
	DECLARE_EVENT_TYPE(wxEVT_IndexingThread, -1)//type for thread to communicate with frame
END_DECLARE_EVENT_TYPES()

class IndexingFrame : public wxFrame {
	private:
		std::thread                  idxThd ;//thread for indexing
		std::atomic_flag             flg    ;//flag to tell worker thread to stop
		emsphinx::ebsd::Namelist     nml    ;//namelist
		wxImage                      imOrig ;//original image (e.g. iq)
		wxImage                      imPaint;//ipf painted image

		enum IdxStatus {
			IdxError = -  1,
			IdxExit  = -  2,
			IdxCmplt =  101,
		};


	protected:
		wxSplitterWindow* m_split  ;
		wxImagePanel    * m_imPan  ;
		EbsdSummaryPanel* m_sumPan ;
		wxGauge         * m_prog   ;
		wxButton        * m_btn    ;
		wxStatusBar     * m_statBar;

		virtual void OnBtn( wxCommandEvent& event ) {
			if("Start" == m_btn->GetLabel()) startIdx();
			else stopIdx();
		}

		void OnThread(wxCommandEvent& event) {
			//update the status bar and get event type
			m_statBar->SetStatusText(event.GetString());
			int i = event.GetInt();

			//handle event type
			if(IdxError == i) {//an error
				//nothing else to do for now (wait for exit message)
				wxMessageDialog msgDlg(this, event.GetString(), "Error Indexing");
				msgDlg.ShowModal();//make sure user sees message (since it will be overwritten with "Indexing Stopped")
			} else if(IdxCmplt == i || IdxExit == i) {//(un)successful indexing completion
				idxThd.join();//wait for the thread to wrap up
				m_sumPan->EnableEditing(true);//enable the parameter panel
				m_btn->SetLabel("Start");//change stop button to a start button again
				m_btn->Enable(true);
			} else {//status update
				m_btn->Enable(true);
				m_prog->SetValue(event.GetInt());//update progress bar
				m_imPan->setImage(imPaint);//update image
				m_imPan->Refresh();//repaint
				m_imPan->paintNow();
			}
		}

	public:

		void startIdx();

		void stopIdx();

		void setNml(emsphinx::ebsd::Namelist& n) {nml = n; m_sumPan->setNamelist(&nml);}

		void setImage(wxImage& im) {m_imPan->setImage(im); imOrig = im; imPaint = im.Copy();}

		IndexingFrame( wxWindow* parent, wxWindowID id = wxID_ANY, const wxString& title = "EMSphInx", const wxPoint& pos = wxDefaultPosition, const wxSize& size = wxSize( 720,540 ), long style = wxDEFAULT_FRAME_STYLE|wxTAB_TRAVERSAL );

		~IndexingFrame();

		wxDECLARE_EVENT_TABLE();
};

BEGIN_EVENT_TABLE(IndexingFrame, wxFrame)
    EVT_COMMAND(wxID_ANY, wxEVT_IndexingThread, IndexingFrame::OnThread)
END_EVENT_TABLE()

///////////////////////////////////////////////////////////////////////////
// C++ code generated with wxFormBuilder (version Oct 26 2018)
// http://www.wxformbuilder.org/
//
// PLEASE DO *NOT* EDIT THIS FILE!
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////

DEFINE_EVENT_TYPE(wxEVT_IndexingThread)

void IndexingFrame::startIdx() {
	imPaint = imOrig.Copy();
	m_sumPan->EnableEditing(false);//disable the parameter panel
	m_sumPan->updateNamelist();
	m_btn->SetLabel("Stop");//change start button to stop button
	m_btn->Enable(false);

	flg.test_and_set();
	idxThd = std::thread([&](){
		//create event for dispatching
		wxCommandEvent evt(wxEVT_IndexingThread);
		evt.SetInt(0);
		evt.SetString("Initializing Indexers");
		wxPostEvent(this, evt);//update image + progress bar
		try {
			emsphinx::ebsd::IndexingData<double> idxData(nml, (char*)imPaint.GetData());//read inputs and build indexers
			ThreadPool pool(idxData.threadCount);//build thread pool
			try {
				//queue parallel indexing
				time_t tmStart = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
				Timer t;
				size_t batches = idxData.pat->numPat() / nml.batchSize;//how many batches are needed
				if(batches * nml.batchSize < idxData.pat->numPat()) ++batches;//extra batch for leftovers
				for(size_t i = 0; i < batches; i++) {//loop over batches
					const size_t start = i * nml.batchSize;//first pattern
					const size_t end = std::min(start + nml.batchSize, idxData.pat->numPat());//last pattern
					pool.schedule(std::bind(idxData.workItem, std::placeholders::_1));//queue indexing
				}

				//wait for completion posting updates periodically
				const std::chrono::milliseconds uptFreq(250);//milliseconds between updates (too long -> unresponsive, too short -> resource intensive), 0.25 ~human reaction time
				while(!pool.waitAll(uptFreq)) {
					if(!flg.test_and_set()) {
						//post error status update
						evt.SetString("Stop Requested");
						wxPostEvent(this, evt);
						pool.clear();//clear all unstarted items
						evt.SetString("Indexing Stopped");
						evt.SetInt(0);
						wxPostEvent(this, evt);//update image + progress bar
						evt.SetInt(IdxExit);
						wxPostEvent(this, evt);//mark exit
						return;
					}

					//get the time elapsed since we started (without resetting the reference point)
					const double elapsed = t.poll(false);
					const double percent = double(idxData.idxCtr) / idxData.numIdx;
					const double rate = elapsed / percent;//estimate seconds per %
					const double remaining = rate * (1.0 - percent);//estimate seconds left

					//print update
					std::ostringstream ss;
					Timer::PrintSeconds(elapsed  , ss);
					ss << " elapsed, ";
					Timer::PrintSeconds(remaining, ss);
					ss << " remaining   ";
					evt.SetString(ss.str());
					evt.SetInt(std::min((int)(percent * 100 + 0.5), 100));
					wxPostEvent(this, evt);
				}
				
				//compute indexing time and save
				const double total = t.poll();
				time_t tmEnd = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
				idxData.save(tmStart, tmEnd, total);
			
				//post completion status
				evt.SetString("Indexing Complete");
				evt.SetInt(100);
				wxPostEvent(this, evt);//post status update with 100% completion
				evt.SetInt(IdxCmplt);
				wxPostEvent(this, evt);//post actual completion message
				return;
			 } catch (std::exception& e) {//clear items and post
				pool.clear();//clear all unstarted items
				//post error message
				evt.SetString(e.what());
				evt.SetInt(IdxError);
				wxPostEvent(this, evt);
			}
			//pool goes out of scope here (waiting for completion)
		} catch (std::exception& e) {
			//post error message
			evt.SetString(e.what());
			evt.SetInt(IdxError);
			wxPostEvent(this, evt);
		} catch (...) {
			evt.SetString("Unhandled exception type during indexing, please contact developers");
			evt.SetInt(IdxError);
			wxPostEvent(this, evt);
		}

		//post exit message
		evt.SetString("Indexing Stopped");
		evt.SetInt(0);
		wxPostEvent(this, evt);//update image + progress bar
		evt.SetInt(IdxExit);
		wxPostEvent(this, evt);//mark exit
	});
}

void IndexingFrame::stopIdx() {
	flg.clear();//tell any running indexers to stop
	m_sumPan->EnableEditing(true);//enable the parameter panel
	m_btn->Enable(false);//disable the start/stop button
}

IndexingFrame::IndexingFrame( wxWindow* parent, wxWindowID id, const wxString& title, const wxPoint& pos, const wxSize& size, long style ) : wxFrame( parent, id, title, pos, size, style ) {
	this->SetSizeHints( wxDefaultSize, wxDefaultSize );

	//create splitter window
	m_split = new wxSplitterWindow( this, wxID_ANY, wxDefaultPosition, wxDefaultSize, wxSP_3DSASH|wxSP_LIVE_UPDATE );
	m_split->SetSashGravity( 1 );//put all growth into left side
	m_split->SetMinimumPaneSize( 50 );

	//create image + indexing parameters and add to splitter
	m_imPan  = new wxImagePanel    ( m_split, wxID_ANY,              wxDefaultPosition, wxDefaultSize, wxTAB_TRAVERSAL );
	m_sumPan = new EbsdSummaryPanel( m_split, wxID_ANY,              wxDefaultPosition, wxDefaultSize, wxTAB_TRAVERSAL );
	m_split->SplitVertically(m_imPan, m_sumPan);

	//create sizers and assemble
	wxBoxSizer* vSizer = new wxBoxSizer( wxVERTICAL   );
	wxBoxSizer* hSizer = new wxBoxSizer( wxHORIZONTAL );
	vSizer->Add( m_split, 1, wxEXPAND, 5 );
	vSizer->Add( hSizer , 0, wxEXPAND, 5 );

	//build and assemble controls
	m_prog   = new wxGauge         ( this, wxID_ANY, 100        , wxDefaultPosition, wxDefaultSize, wxGA_HORIZONTAL );
	m_btn    = new wxButton        ( this, wxID_ANY, wxT("Start"), wxDefaultPosition, wxDefaultSize, 0 );
	m_prog->SetValue( 0 );
	m_prog->SetRange( 100 );
	hSizer->Add( m_prog  , 1, wxEXPAND, 0 );
	hSizer->Add( m_btn   , 0, wxEXPAND, 0 );

	//layout
	this->SetSizer( vSizer );
	this->Layout();
	m_statBar = this->CreateStatusBar( 1, wxSTB_ELLIPSIZE_END|wxSTB_SIZEGRIP, wxID_ANY );
	this->Centre( wxBOTH );

	m_btn->Connect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( IndexingFrame::OnBtn ), NULL, this );
}

IndexingFrame::~IndexingFrame() {
	m_btn->Disconnect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( IndexingFrame::OnBtn ), NULL, this );
}

#endif//_IDX_FRAME_H

