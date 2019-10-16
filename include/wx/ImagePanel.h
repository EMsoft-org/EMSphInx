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

#ifndef _IMAGE_PANEL_H
#define _IMAGE_PANEL_H

#include <wx/frame.h>
#include <wx/panel.h>
#include <wx/sizer.h>
#include <wx/dcclient.h>

//@brief: helper class to dynamically draw a rescaled image into a frame (preserving aspect ratio)
//@note : based on //https://wiki.wxwidgets.org/An_image_panel
class wxImagePanel : public wxPanel {
	wxImage  image     ;//the original image
	wxBitmap resized   ;//the resized image (to render)
	int      curW, curH;//the current width/height of the resized image
	bool     stale     ;//is a redraw required even if the size hasn't changed (e.g. the image was changed)
	wxPoint  ori       ;//current origin of image
	
	public:
		//@brief       : constructor, mirrors signature of wxPanel constructor
		//@param parent: The parent window.
		//@param id    : An identifier for the panel. wxID_ANY is taken to mean a default.
		//@param pos   : The panel position. The value wxDefaultPosition indicates a default position, chosen by either the windowing system or wxWidgets, depending on platform.
		//@param size  : The panel size. The value wxDefaultSize indicates a default size, chosen by either the windowing system or wxWidgets, depending on platform.
		//@param style : The window style. See wxPanel.
		//@param name  : Window name.
		wxImagePanel(wxWindow *parent, wxWindowID id=wxID_ANY, const wxPoint &pos=wxDefaultPosition, const wxSize &size=wxDefaultSize, long style=wxTAB_TRAVERSAL, const wxString &name=wxPanelNameStr);
		
		//@brief   : update the image to draw
		//@param im: new image
		void setImage(wxImage& im) {image = im.Copy(); if(!im.HasAlpha()) image.SetAlpha(); stale = true;}

		//@brief : check if there is an image
		//@return: true if there is an image, false otherwise
		bool hasImage() const {return image.IsOk();}

		//@brief  : update a single pixel of the image
		//@param x: x index of pixel to update
		//@param y: y index of pixel to update
		//@param r: red value
		//@param g: green value
		//@param b: blue value
		void SetRGB(const size_t x, const size_t y, const char r, const char g, const char b) {image.SetRGB(x, y, r, g, b); stale = true;}

		//@brief  : update a single pixel of the image
		//@param x: x index of pixel to update
		//@param y: y index of pixel to update
		//@param v: gray value
		void SetGray(const size_t x, const size_t y, const char v) {SetRGB(x, y, v, v, v);}

		//@brief  : update a single pixel of the image
		//@param x: x index of pixel to update
		//@param y: y index of pixel to update
		//@param v: alpha value
		void SetAlpha(const size_t x, const size_t y, const char v) {image.SetAlpha(x, y, v); stale = true;}

		//@brief: manually mark the image as outdated
		void markStale() {stale = true;}

		//@brief : get the underlying rgb buffer
		//@return: rgb buffer (you are responsible for calling markStale() if you change values)
		unsigned char* GetRGB() {return image.GetData();}

		//@brief : get the underlying alpha buffer
		//@return: alpha buffer (you are responsible for calling markStale() if you change values)
		unsigned char* GetAlpha() {return image.GetAlpha();}

		//@brief : get original (unscaled) image size
		//@return: image size
		wxSize GetSize() const {return image.IsOk() ? image.GetSize() : wxSize(0, 0);}

		//@brief : get the current origin of the image
		//@return: image origin
		wxPoint getOri() const {return ori;}

		//@brief : get the current scale of the image
		//@return: image scaling
		double getScl() const {return image.IsOk() ? double(curW) / image.GetSize().GetWidth() : 0.0;}

		//@brief : get the current scale of the image^2
		//@return: image scaling^2
		double getScl2() const {return image.IsOk() ? double(curW * curH) / (image.GetSize().GetWidth() * image.GetSize().GetHeight()): 0.0;}

		//@brief  : convert from image coordinates to DC coordinates
		//@param x: x coordinate to convert from image -> panel coordinates
		//@param y: x coordinate to convert from image -> panel coordinates
		void im2dc(int& x, int& y) const;

		//@brief: convert from image coordinates to DC coordinates
		//@param x: x coordinate to convert from panel -> image coordinates
		//@param y: x coordinate to convert from panel -> image coordinates
		void dc2im(int& x, int& y) const;

		//@brief   : actual code to draw the image
		//@param dc: device context to draw with
		//@note    : seperate function to handle different context types
		void render(wxDC& dc);
		
		//@brief     : called by wxWidgets event loop to redraw
		//@param evt: redraw event
		//@note     : triggered by Refresh()/Update()
		void paintEvent(wxPaintEvent & evt);

		//@brief: manually redraw
		void paintNow();

		//@brief      : event for image resizing (refresh panel)
		//@param event: resize event
		void OnSize(wxSizeEvent& event);

		DECLARE_EVENT_TABLE()
};

BEGIN_EVENT_TABLE(wxImagePanel, wxPanel)
	EVT_PAINT(wxImagePanel::paintEvent)// paint events
	EVT_SIZE (wxImagePanel::OnSize    )// size event
END_EVENT_TABLE()

//@brief       : image panel constructor, mirrors signature of wxPanel constructor
//@param parent: The parent window.
//@param id    : An identifier for the panel. wxID_ANY is taken to mean a default.
//@param pos   : The panel position. The value wxDefaultPosition indicates a default position, chosen by either the windowing system or wxWidgets, depending on platform.
//@param size  : The panel size. The value wxDefaultSize indicates a default size, chosen by either the windowing system or wxWidgets, depending on platform.
//@param style : The window style. See wxPanel.
//@param name  : Window name.
wxImagePanel::wxImagePanel(wxWindow *parent, wxWindowID id, const wxPoint &pos, const wxSize &size, long style, const wxString &name) : wxPanel(parent, id, pos, size, style, name), image(), curW(-1), curH(-1), stale(true) {
	SetBackgroundColour( parent->GetBackgroundColour() );//transparent background
}

//@brief  : convert from image coordinates to DC coordinates
//@param x: x coordinate to convert from image -> panel coordinates
//@param y: x coordinate to convert from image -> panel coordinates
void wxImagePanel::im2dc(int& x, int& y) const {
	x = (int)std::round( double(x * curW) / image.GetSize().GetWidth () ) + ori.x;
	y = (int)std::round( double(y * curH) / image.GetSize().GetHeight() ) + ori.y;
}

//@brief: convert from image coordinates to DC coordinates
//@param x: x coordinate to convert from panel -> image coordinates
//@param y: x coordinate to convert from panel -> image coordinates
void wxImagePanel::dc2im(int& x, int& y)const  {
	x = (int) std::round( double( (x - ori.x) * image.GetSize().GetWidth () ) / curW );
	y = (int) std::round( double( (y - ori.y) * image.GetSize().GetHeight() ) / curH );
}

//@brief   : actual code to draw the image
//@param dc: device context to draw with
//@note    : seperate function to handle different context types
void wxImagePanel::render(wxDC& dc) {
	//we we don't have an image we're done
	if(!image.IsOk()) return;

	//get frame size
	int dcW, dcH;
	wxWindowDC(this).GetSize( &dcW, &dcH );
	// dc.GetSize( &dcW, &dcH );//this doesn't awlays update correctly for all sizes on windows

	//get image size
	int imW = image.GetSize().GetWidth ();
	int imH = image.GetSize().GetHeight();

	//get max image size that fits in frame with same aspect ratio
	int newW = dcW;
	int newH = dcH;
	float rW = float(newW) / imW;
	float rH = float(newH) / imH;
	if(rW <= rH) {
		newH = (int) std::round(rW * imH);
	} else {
		newW = (int) std::round(rH * imW);
	}

	//update bitmap if needed
	if(newW > 0 && newH > 0) {

		if( newW != curW || newH != curH || stale) {
			if(0 == newW || 0 == newH) {
				resized = wxBitmap();//rescaled to nothing
			} else {
				resized = wxBitmap( image.Scale( newW, newH /*, wxIMAGE_QUALITY_HIGH*/ ) );
			}
			curW = newW;
			curH = newH;
			stale = false;
		}
	}

	//draw
	int dx = (dcW - newW) / 2;
	int dy = (dcH - newH) / 2;
	ori.x = dx;
	ori.y = dy;
	if(curW == 0 || curH == 0) return;//nothing to draw
	dc.DrawBitmap( resized, dx, dy, false );
}

//@brief     : called by wxWidgets event loop to redraw
//@param evt: redraw event
//@note     : triggered by Refresh()/Update()
void wxImagePanel::paintEvent(wxPaintEvent & evt) {
	// depending on your system you may need to look at double-buffered dcs
	wxPaintDC dc(this);
	render(dc);
}

//@brief: manually redraw
void wxImagePanel::paintNow() {
	// depending on your system you may need to look at double-buffered dcs
	wxClientDC dc(this);
	render(dc);
}

//@brief      : event for image resizing (refresh panel)
//@param event: resize event
void wxImagePanel::OnSize(wxSizeEvent& event){
	Refresh();
	event.Skip();//skip the event.
}

#endif//_IMAGE_PANEL_H
