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

#ifndef _ROI_SEL_PANEL_H_
#define _ROI_SEL_PANEL_H_

#include <wx/panel.h>
#include <wx/frame.h>
#include <wx/dcbuffer.h>

#include "idx/roi.h"
#include "wx/ImagePanel.h"

#include <wx/event.h>

//@brief: an event to be emmited when the shape changes
class DrawShapeEvent : public wxCommandEvent {
	std::vector< std::pair<int, int> > const* m_shape;
	public:
		DrawShapeEvent(wxEventType commandEventType = wxEVT_NULL, int id = 0) : wxCommandEvent(commandEventType, id) {}
		void SetShape(const std::vector< std::pair<int, int> >& shp) {m_shape = &shp;}
		std::vector< std::pair<int, int> > const* GetShape() {return m_shape;}
		virtual wxEvent *Clone() const { return new DrawShapeEvent(*this); }
};
wxDEFINE_EVENT(DRAW_SHAPE, DrawShapeEvent);
typedef void (wxEvtHandler::*DrawShapeEventFunction)(DrawShapeEvent&);
#define DrawShapeEventHandler(func) wxEVENT_HANDLER_CAST(DrawShapeEventFunction, func)
#define EVT_DRAW_SHAPE(id, func) wx__DECLARE_EVT1(DRAW_SHAPE, id, DrawShapeEventHandler(func))

//@brief: a panel to display + manipulate a region of interest
class ImageRoiEditor : public wxImagePanel {
	emsphinx::RoiSelection sel       ;
	std::vector<char>      alphaMask ;
	int                    stroke    ;
	int                    cpRd      ;//capture radius for nearby point determination
	wxPoint                lstMouse  ;//location of last mouse movement for click+drag
	wxColor                penClr    ;

	//@brief: update mask + refresh
	void updateMask();

	//@brief: emit a shape changed event
	void emitShape();

	// event handling
	void leftClick    (wxMouseEvent& event);//start new shape, move existing handle, or click + drag shape for outside shape, inside shape, or on handle; add new handle if double click on handle
	void rightClick   (wxMouseEvent& event);//finish polygon or remove polygon node
	void mouseMoved   (wxMouseEvent& event);//move shape or handle, update cursor for user hint near handles
	void mouseReleased(wxMouseEvent& event);//stop modifying
	void keyPressed   (wxKeyEvent  & event);//fix aspect ratio/line direction while shift is held, delete nodes with backspace/escape
	void keyReleased  (wxKeyEvent  & event);//clear shift held
	void onSize       (wxSizeEvent & event);//resize image + selection then refresh
	void paintEvent   (wxPaintEvent& event);//render

	DECLARE_EVENT_TABLE()

	public:

		//@brief       : constructor
		//@param parent: parent frame
		ImageRoiEditor(wxWindow* parent, wxWindowID id=wxID_ANY) : wxImagePanel(parent, id), stroke(2), cpRd(4), lstMouse(-1, -1) {this->SetBackgroundStyle(wxBG_STYLE_PAINT);}
		
		//@brief   : update the image to draw
		//@param im: new image
		void setImage(wxImage& im) {wxImagePanel::setImage(im); updateMask();}

		//@brief   : change the drawing mode
		//@param dm: new drawing mode
		void setDrawMode(const emsphinx::DrawMode dm) {if(sel.changeMode(dm)) {updateMask();}}
		
		//@brief    : set if the selected area is included or excluded
		//@param inv: true/false to select inside/outside shape
		void setInvert(const bool inv) {sel.setInv(inv); if(!alphaMask.empty()) {updateMask();}}
		
		//@brief    : set the ROI outline line color
		//@param clr: color to draw lines with
		void setColor(wxColor clr) {penClr = clr; Refresh();}

		//@brief    : set the outline thickness
		//@param stk: thickness
		void setStroke(int stk) {stroke = stk; cpRd = 3 * stk / 2; Refresh();}

		//@brief   : render the drawing
		//@param dc: context to render onto
		void render(wxDC& dc);

		//@brief    : update coordinates of a point
		//@param idx: point to update
		//@param x  : new x coordiante
		//@param y  : new y coordiante
		void movePoint(const size_t idx, const int x, const int y) {sel.movePoint(idx, x, y); updateMask();}

		//@brief : get ROI
		//@return: roi
		emsphinx::RoiSelection getRoi() const {return sel;}

		//@brief    : set ROI
		//@param roi: roi
		void setRoi(emsphinx::RoiSelection& roi) {sel = roi; updateMask(); Refresh();}
};

BEGIN_EVENT_TABLE(ImageRoiEditor, wxPanel)
	EVT_LEFT_DOWN  (ImageRoiEditor::leftClick    )
	EVT_LEFT_DCLICK(ImageRoiEditor::leftClick    )
	EVT_RIGHT_DOWN (ImageRoiEditor::rightClick   )
	EVT_MOTION     (ImageRoiEditor::mouseMoved   )
	EVT_LEFT_UP    (ImageRoiEditor::mouseReleased)
	EVT_KEY_DOWN   (ImageRoiEditor::keyPressed   )
	EVT_KEY_UP     (ImageRoiEditor::keyReleased  )
	EVT_SIZE       (ImageRoiEditor::onSize       )
	EVT_PAINT      (ImageRoiEditor::paintEvent   )
END_EVENT_TABLE()

#include <wx/choice.h>
#include <wx/sizer.h>
#include <wx/checkbox.h>
#include <wx/clrpicker.h>
#include <wx/stattext.h>
#include <wx/statline.h>
#include <wx/spinctrl.h>
#include <wx/grid.h>

#include <sstream>

enum {
	ID_CHOICE = 1000,
	ID_COLOR        ,
	ID_SPIN         ,
	ID_CHECK        ,
	ID_DRAW         ,
};

//@brief: panel to hold ImageRoiEditorl with additional controls
class RoiSelectionPanel : public wxPanel {

	wxGrid            * grid    ;
	ImageRoiEditor    * drawPane;
	wxChoice          * choice  ;
	wxColourPickerCtrl* clrPick ;
	wxSpinCtrl        * spinCtrl;
	wxCheckBox        * checkBox;

	std::vector< std::pair<int, int> > lastShape;

	// event handling
	void choiceChanged( wxCommandEvent     & event );
	void colorChanged ( wxColourPickerEvent& event ) {drawPane->setColor (event    .GetColour  ());}
	void spinChanged  ( wxSpinEvent        & event ) {drawPane->setStroke(event    .GetPosition());}
	void invToggled   ( wxCommandEvent     & event ) {drawPane->setInvert(checkBox->GetValue   ());}
	void shapeChanged ( DrawShapeEvent     & event );
	void gridChanged  ( wxGridEvent        & event );

	public:
		RoiSelectionPanel( wxWindow* parent );

		void setImage(wxImage& im);

		//@brief : get ROI
		//@return: roi
		emsphinx::RoiSelection getRoi() const {return drawPane->getRoi();}

		//@brief    : set ROI
		//@param roi: roi
		void setRoi(emsphinx::RoiSelection& roi) {checkBox->SetValue(roi.getInv()); drawPane->setRoi(roi);}

		DECLARE_EVENT_TABLE()
};

BEGIN_EVENT_TABLE(RoiSelectionPanel, wxPanel)
	EVT_CHOICE              (ID_CHOICE, RoiSelectionPanel::choiceChanged)
	EVT_COLOURPICKER_CHANGED(ID_COLOR , RoiSelectionPanel::colorChanged )
	EVT_SPINCTRL            (ID_SPIN  , RoiSelectionPanel::spinChanged  )
	EVT_CHECKBOX            (ID_CHECK , RoiSelectionPanel::invToggled   )
    EVT_DRAW_SHAPE          (ID_DRAW  , RoiSelectionPanel::shapeChanged )
    EVT_GRID_CELL_CHANGED   (           RoiSelectionPanel::gridChanged  )
END_EVENT_TABLE()

///////////////////////////////////////////////////////////////////////
//                          ImageRoiEditor                           //
///////////////////////////////////////////////////////////////////////

//@brief: update mask + refresh
void ImageRoiEditor::updateMask() {
	//build mask
	if(wxImagePanel::hasImage()) {
		const size_t w = wxImagePanel::GetSize().GetWidth ();
		const size_t h = wxImagePanel::GetSize().GetHeight();
		alphaMask = sel.buildMask(w, h);
		if(!sel.hasShape()) std::fill(alphaMask.begin(), alphaMask.end(), sel.getInv() ? 0 : 1);
		for(size_t j = 0; j < h; j++) {
			for(size_t i = 0; i < w; i++) {
				bool in = alphaMask[j*w+i] == 1;
				wxImagePanel::SetAlpha(i, j, in ? 0xFF : 0x40);
			}
		}
	} else {
		alphaMask.clear();
	}
	emitShape();
	Refresh();
}

//@brief: emit a shape changed event
void ImageRoiEditor::emitShape() {
	DrawShapeEvent event(DRAW_SHAPE, GetId());
	event.SetEventObject(this);
	event.SetShape(sel.getPts());
	ProcessWindowEvent(event);//send
}

void ImageRoiEditor::leftClick(wxMouseEvent& event) {
	wxPoint xy(event.GetX(), event.GetY());
	wxImagePanel::dc2im(xy.x, xy.y);
	if(sel.hasSelection()) {//not currently manipulating an existing point
		if(sel.trySelect(xy.x, xy.y, (int)std::round(1.0 / wxImagePanel::getScl2() * cpRd * cpRd))) {
			//did click on an existing point 
			switch(sel.getMode()) {
				case emsphinx::DrawMode::Rectangle:
				case emsphinx::DrawMode::Ellipse  : break;//we found it, we're done
				case emsphinx::DrawMode::Polygon  :
					if(event.LeftDClick()) {
						sel.dupPolyPt((int)std::round(2.0 / wxImagePanel::getScl() * cpRd));//duplicate selected point
						updateMask();
					} else {
						//just selecting point was enough
					}
				break;
				sel.startPoly(xy.x, xy.y); updateMask(); break;
			}
		} else if(sel.inside(xy.x, xy.y) != sel.getInv()) {
			lstMouse = wxPoint(xy.x, xy.y);
		} else {//user didn't try click on an existing point or inside existing shape, start a new shape
			switch(sel.getMode()) {
				case emsphinx::DrawMode::Rectangle:
				case emsphinx::DrawMode::Ellipse  : sel.startRect(xy.x, xy.y); updateMask(); break;
				case emsphinx::DrawMode::Polygon  : sel.startPoly(xy.x, xy.y); updateMask(); break;
			}
			Refresh();
		}
		this->SetCursor(*wxCROSS_CURSOR);//use cross during drawing
	} else {//currently manipulating an existing point
		switch(sel.getMode()) {
			case emsphinx::DrawMode::Rectangle:
			case emsphinx::DrawMode::Ellipse  : break;
			case emsphinx::DrawMode::Polygon  : if(sel.addPolyPt(xy.x, xy.y, (int)std::round(1.0 / wxImagePanel::getScl2() * cpRd * cpRd)) ) updateMask(); break;
		}
	}
	event.Skip();
}

void ImageRoiEditor::rightClick(wxMouseEvent& event) {
	wxPoint xy(event.GetX(), event.GetY());
	wxImagePanel::dc2im(xy.x, xy.y);
	switch(sel.getMode()) {
		case emsphinx::DrawMode::Rectangle:
		case emsphinx::DrawMode::Ellipse  : break;
		case emsphinx::DrawMode::Polygon  : 
			if(sel.buildingPoly()) {
				sel.addPolyPt(xy.x, xy.y, cpRd * cpRd);//add the current point
				sel.finishPoly();//close polygon
				updateMask();//update mask
				// this->SetCursor(*wxSTANDARD_CURSOR);//switch cursor to indicate close
			} else if(sel.trySelect(xy.x, xy.y, (int)std::round(1.0 / wxImagePanel::getScl2() * cpRd * cpRd)) ) {
				sel.remPolyPt();
				updateMask();
			}
		break;
	}
	event.Skip();
}

void ImageRoiEditor::mouseMoved(wxMouseEvent& event) {
	wxPoint xy(event.GetX(), event.GetY());
	wxImagePanel::dc2im(xy.x, xy.y);
	if(-1 != lstMouse.x && -1 != lstMouse.y) {//clicking + dragging
		sel.translatePoints(xy.x - lstMouse.x, xy.y - lstMouse.y);//drag selection
		lstMouse.x = xy.x;//update most recent x
		lstMouse.y = xy.y;//update most recent y
		updateMask();
	} else {
		switch(sel.getMode()) {
			case emsphinx::DrawMode::Rectangle:
			case emsphinx::DrawMode::Ellipse  : if(sel.resizeRect(xy.x, xy.y)) updateMask(); break;
			case emsphinx::DrawMode::Polygon  : if(sel.resizePoly(xy.x, xy.y)) updateMask(); break;
		}
		if(sel.hasSelection()) {//we're not currently building/manipulating
			this->SetCursor(-1 == sel.nearPt(xy.x, xy.y, (int)std::round(1.0 / wxImagePanel::getScl2() * cpRd * cpRd)) ? *wxSTANDARD_CURSOR : *wxCROSS_CURSOR);//use cross to indicate near enough to grab handle
		} else if (emsphinx::DrawMode::Polygon == sel.getMode() && sel.buildingPoly()) {//handle in construction polygon specially
			std::pair<int, int> pt = sel.getPts().front();
			wxImagePanel::im2dc(pt.first, pt.second);
			wxPoint dlt(event.GetX() - pt.first, event.GetY() - pt.second);
			const bool close = dlt.x*dlt.x + dlt.y*dlt.y < cpRd * cpRd;
			this->SetCursor(close ? *wxSTANDARD_CURSOR : *wxCROSS_CURSOR);//use regular cursor to indicate near enough to close polygon
		}
	}
	event.Skip();
}

void ImageRoiEditor::mouseReleased(wxMouseEvent& event) {
	switch(sel.getMode()) {
		case emsphinx::DrawMode::Rectangle:
		case emsphinx::DrawMode::Ellipse  : sel.ungrabRect(); updateMask(); break;
		case emsphinx::DrawMode::Polygon  : sel.ungrabPoly(); updateMask(); break;
	}
	this->SetCursor(sel.hasSelection() ? *wxSTANDARD_CURSOR : *wxCROSS_CURSOR);//switch to regular cursor when done editing
	lstMouse = wxPoint(-1, -1);
	event.Skip();
}

void ImageRoiEditor::keyPressed(wxKeyEvent& event) {
	switch(event.GetKeyCode()) {
		case WXK_DELETE:
		case WXK_BACK  :
			switch(sel.getMode()) {
				case emsphinx::DrawMode::Rectangle:
				case emsphinx::DrawMode::Ellipse  : break;
				case emsphinx::DrawMode::Polygon  : if(sel.remPolyPt()) {emitShape(); Refresh();} break;
			}
		break;

		case WXK_ESCAPE: if(sel.clear()) updateMask(); break;

		case WXK_SHIFT : sel.setFixedAspect(true); emitShape(); event.Skip(); break;

		default: event.Skip();
	}
}

void ImageRoiEditor::keyReleased(wxKeyEvent& event) {
	switch(event.GetKeyCode()) {
		case WXK_SHIFT: sel.setFixedAspect(false); emitShape(); event.Skip(); break;

		default: event.Skip();
	}
}

//@brief      : event for image resizing (refresh panel)
//@param event: resize event
void ImageRoiEditor::onSize(wxSizeEvent& event) {
	lstMouse = wxPoint(-1, -1);//invalidate click+drag
	Refresh();
	// event.Skip();//skip the event.
}

void ImageRoiEditor::paintEvent(wxPaintEvent & event) {
	wxAutoBufferedPaintDC dc(this);
	render(dc);
}

void ImageRoiEditor::render(wxDC& dc) {
	//scale points
	std::vector<wxPoint> pts;
	for(const std::pair<int, int>& p : sel.getPts()) {
		int x = p.first ;
		int y = p.second;
		wxImagePanel::im2dc(x, y);
		pts.push_back(wxPoint(x, y));
	}

	//draw based on mode
	dc.Clear();
	wxImagePanel::render(dc);//draw the image first
	dc.SetBrush( *wxTRANSPARENT_BRUSH );
	dc.SetPen( wxPen(penClr , stroke) );
	switch(sel.getMode()) {
		case emsphinx::DrawMode::Rectangle:
		case emsphinx::DrawMode::Ellipse  : {
			if(!pts.empty()) {
				//draw rectangle
				const int x0 = std::min(pts.front().x, pts.back().x);
				const int y0 = std::min(pts.front().y, pts.back().y);
				const int w  = std::max(pts.front().x, pts.back().x) + 1 - x0;
				const int h = std::max(pts.front().y, pts.back().y) + 1 - y0;
				if(emsphinx::DrawMode::Rectangle == sel.getMode() || 0 == w || 0 == h)//gtk doesn't like drawing empty ellipses
					dc.DrawRectangle( x0, y0, w, h );
				else
					dc.DrawEllipse  ( x0, y0, w, h );
			}
		} break;

		case emsphinx::DrawMode::Polygon  :
			if(!pts.empty()) dc.DrawLines(pts.size(), pts.data());
		break;
	}

	//draw handles
	dc.SetBrush( *wxWHITE_BRUSH );
	dc.SetPen( *wxBLACK_PEN );
	for(const auto& pt : pts) dc.DrawCircle(pt, cpRd);
}

///////////////////////////////////////////////////////////////////////
//                         RoiSelectionPanel                         //
///////////////////////////////////////////////////////////////////////

void RoiSelectionPanel::choiceChanged( wxCommandEvent& event ) {
	switch(choice->GetSelection()) {
		case 0: drawPane->setDrawMode(emsphinx::DrawMode::Rectangle); break;
		case 1: drawPane->setDrawMode(emsphinx::DrawMode::Ellipse  ); break;
		case 2: drawPane->setDrawMode(emsphinx::DrawMode::Polygon  ); break;
	}
	event.Skip();
}

void RoiSelectionPanel::shapeChanged ( DrawShapeEvent     & event ) {
	//first resize grid if needed
	size_t numRow = grid->GetNumberRows();//current grid rows
	size_t tarRow = event.GetShape()->size();//target grid rows
	if(tarRow > 3) {//only possible for polygons
		if(event.GetShape()->front() == event.GetShape()->back()) --tarRow;//dont display duplicate first/last point
	}
	if(tarRow < numRow) {
		grid->DeleteRows(tarRow, numRow - tarRow);
	} else if(tarRow > numRow) {
		grid->AppendRows(tarRow - numRow);
	}

	//next update grid entries / saved shape
	lastShape.resize(tarRow, std::pair<int,int>(-1,-1));
	for(size_t i = 0; i < lastShape.size(); i++) {
		if(lastShape[i].first  != event.GetShape()->at(i).first ) {
			lastShape[i].first  = event.GetShape()->at(i).first;
			std::ostringstream ss;
			ss << lastShape[i].first ;
			grid->SetCellValue(i, 0, ss.str().c_str());
		}
		if(lastShape[i].second != event.GetShape()->at(i).second) {
			lastShape[i].second  = event.GetShape()->at(i).second;
			std::ostringstream ss;
			ss << lastShape[i].second;
			grid->SetCellValue(i, 1, ss.str().c_str());
			
		}

	}
}

void RoiSelectionPanel::gridChanged  ( wxGridEvent        & event ) {
	//get row
	int idx = event.GetRow();
	if(0 == event.GetCol())
		lastShape[idx].first  = atoi(grid->GetCellValue(idx, 0));
	else
		lastShape[idx].second = atoi(grid->GetCellValue(idx, 1));
	drawPane->movePoint(idx, lastShape[idx].first, lastShape[idx].second);
}

void RoiSelectionPanel::setImage(wxImage& im) {
	//get the minimum size of the panel (from contrl bar) and update minimum height to reflect images
	const wxSize& sz = this->GetSizer()->GetMinSize();
	double numIm = double(sz.GetWidth() - grid->GetSize().GetWidth()) / im.GetSize().GetWidth();//this is how many times we could fit the full image
	int imH = std::round(numIm * im.GetSize().GetHeight());//this is how tall the images should be for the widhts to nicely fill the frame
	this->GetSizer()->SetMinSize(wxSize(sz.GetWidth(), sz.GetHeight() + imH));//make images fill the frame
	this->SetMinSize(wxSize(sz.GetWidth(), sz.GetHeight() + imH));//make images fill the frame
	this->GetSizer()->Fit(GetParent());//force the parent window to grow to accomodate our new size

	//actually update the image
	drawPane->setImage(im);
}

RoiSelectionPanel::RoiSelectionPanel( wxWindow* parent ) : wxPanel(parent) {
	
	//build sizers
	wxBoxSizer* vSizer    = new wxBoxSizer(wxVERTICAL  );
	wxBoxSizer* drwSizer  = new wxBoxSizer(wxHORIZONTAL);
	wxBoxSizer* cntrSizer = new wxBoxSizer(wxHORIZONTAL);

	//build draw elements
	drawPane          = new ImageRoiEditor    ( this, ID_DRAW );

	grid = new wxGrid( this, wxID_ANY, wxDefaultPosition, wxDefaultSize, 0 );

	// Grid
	grid->CreateGrid        ( 0, 2  );
	// grid->EnableEditing     ( false );
	grid->EnableGridLines   ( true  );
	grid->EnableDragGridSize( false );
	grid->SetMargins        ( 0, 0 );
	grid->SetTabBehaviour   ( wxGrid::TabBehaviour::Tab_Wrap );

	// rows/columns
	grid->EnableDragColMove( false );
	grid->SetColLabelSize( 30 );
	grid->SetRowLabelSize( 30 );
	grid->SetColLabelAlignment( wxALIGN_CENTER, wxALIGN_CENTER );
	grid->SetRowLabelAlignment( wxALIGN_CENTER, wxALIGN_CENTER );
	grid->SetColLabelValue( 0, wxT("X") );
	grid->SetColLabelValue( 1, wxT("Y") );

	// Label Appearance

	// Cell Defaults
	grid->SetDefaultCellAlignment( wxALIGN_LEFT, wxALIGN_TOP );
	grid->SetDefaultEditor(new wxGridCellNumberEditor());


	//build control elements
	wxString choiceStrs[] = {"Rectangle", "Ellipse", "Polygon" };
	choice            = new wxChoice          ( this, ID_CHOICE,                  wxDefaultPosition, wxDefaultSize, 3, choiceStrs, 0 );
	clrPick           = new wxColourPickerCtrl( this, ID_COLOR , *wxRED);
	checkBox          = new wxCheckBox        ( this, ID_CHECK , wxT("Inverted"));
	wxStaticText* txt = new wxStaticText      ( this, wxID_ANY , wxT("Stroke")  , wxDefaultPosition, wxDefaultSize, 0 );//spinner label
	spinCtrl          = new wxSpinCtrl        ( this, ID_SPIN  , wxEmptyString  , wxDefaultPosition, wxDefaultSize, wxSP_ARROW_KEYS, 0, 16, 3);

	//set some defaults
	drawPane->setDrawMode(emsphinx::DrawMode::Rectangle);
	drawPane->setColor(clrPick->GetColour());
	drawPane->setStroke(spinCtrl->GetValue());
	drawPane->setInvert (checkBox->GetValue());
	choice->SetSelection( 0 );

	//build the draw pane
	drwSizer->Add(grid    , 0, wxALL|wxEXPAND, 5 );
	drwSizer->Add(drawPane, 1,       wxEXPAND, 0 );

	//build the controls
	cntrSizer->Add(choice  , 0, wxALL|wxALIGN_CENTER_VERTICAL, 5);
	cntrSizer->Add(checkBox, 0, wxALL|wxALIGN_CENTER_VERTICAL, 5);
	cntrSizer->Add(new wxStaticLine(this, wxID_ANY, wxDefaultPosition, wxDefaultSize, wxLI_VERTICAL), 0, wxALL|wxEXPAND, 5);
	cntrSizer->Add(clrPick , 0, wxALL|wxALIGN_CENTER_VERTICAL, 5);
	cntrSizer->Add(txt     , 0, wxALL|wxALIGN_CENTER_VERTICAL, 5);
	cntrSizer->Add(spinCtrl, 0, wxALL|wxALIGN_CENTER_VERTICAL, 5);

	//and draw + control sizers to main sizer
	vSizer->Add(drwSizer , 1, wxEXPAND                       , 0 );
	vSizer->Add(cntrSizer, 0, wxALL|wxALIGN_CENTER_HORIZONTAL, 5 );
	
	//layout frame
	this->SetSizer(vSizer);
	this->SetAutoLayout(true);
}

#endif//_ROI_SEL_PANEL_H_
