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

#ifndef _VALIDITY_EVENT_H_
#define _VALIDITY_EVENT_H_

#include <string>

#include <wx/event.h>

//@brief: an event to be emmited when the shape changes
class ValidityChangedEvent : public wxCommandEvent {
	bool        m_vld;
	std::string m_msg;
	public:
		ValidityChangedEvent(wxEventType commandEventType = wxEVT_NULL, int id = 0) : wxCommandEvent(commandEventType, id), m_vld(false) {}
		void SetValid(const bool vld) {m_vld = vld;}
		bool GetValid() const {return m_vld;}
		void SetMsg(const std::string msg) {m_msg = msg;}
		std::string GetMsg() const {return m_msg;}
		virtual wxEvent *Clone() const { return new ValidityChangedEvent(*this); }
};
wxDEFINE_EVENT(VALID_CHANGED, ValidityChangedEvent);
typedef void (wxEvtHandler::*ValidityChangedEventFunction)(ValidityChangedEvent&);
#define ValidityChangedEventHandler(func) wxEVENT_HANDLER_CAST(ValidityChangedEventFunction, func)
#define EVT_VALID_CHANGED(id, func) wx__DECLARE_EVT1(VALID_CHANGED, id, ValidityChangedEventHandler(func))


class ValidityPanel : public wxPanel {
	protected:
		bool        m_first ;//has at least one test been processed
		bool        m_valid ;//is the panel currently valid
		std::string m_lstMsg;//last message for validity

		//@brief: test for and essentially emit a validity changed event
		virtual void testValid();

		//@brief: force update the status message
		void updateStatus(std::string str);

	public:
		//@brief : sanity check the current state
		//@return: true if the values parsed from the panel are reasonable, false otherwise
		virtual bool isValid() const {return validMsg().empty();}
		virtual std::string validMsg() const {return "";}

		ValidityPanel( wxWindow* parent, wxWindowID id = wxID_ANY, const wxPoint& pos = wxDefaultPosition, const wxSize& size = wxSize( 400,400 ), long style = wxTAB_TRAVERSAL, const wxString& name = wxEmptyString ) : wxPanel( parent, id, pos, size, style, name ), m_first(true), m_valid(true) {}
		
};

//@brief: emit a validity changed event
void ValidityPanel::testValid() {
	const bool        vld = isValid ();
	const std::string msg = validMsg();
	if(vld != m_valid || msg != m_lstMsg || m_first) {
		m_first = false;
		m_valid  = vld;
		m_lstMsg = msg;
		ValidityChangedEvent event(VALID_CHANGED, GetId());
		event.SetEventObject(this);
		event.SetValid(vld);
		event.SetMsg(msg);
		ProcessWindowEvent(event);//send
	}
}

//@brief: force update the status message
void ValidityPanel::updateStatus(std::string str) {ValidityChangedEvent event(VALID_CHANGED, GetId());
	event.SetEventObject(this);
	event.SetValid(isValid());
	event.SetMsg(str);
	ProcessWindowEvent(event);//send
}

#endif//_VALIDITY_EVENT_H_
