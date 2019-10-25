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

#ifndef _MASTERFILELIST_H_
#define _MASTERFILELIST_H_

#include <wx/listctrl.h>

#include "PeriodicTablePanel.h"

//@brief  : get the address of the function pointer stored by a std::function
//@param f: function to get address of
//@return : function address
template<typename T, typename... U>
size_t getAddress(std::function<T(U...)> f) {
	typedef T(fnType)(U...);
	fnType ** fnPointer = f.template target<fnType*>();
	return (size_t) *fnPointer;
}

//@brief: class to hold information about a master pattern file
struct MasterFile {
	static const wxString SgNames[230];//could do this with hm.hpp but this is easier and probably faster

	//actual data
	wxString    pth;//full path
	wxString    frm;//formula
	wxString    nam;//name
	wxString    syb;//structure symbol
	float       kv ;//voltage
	float       tlt;//tilt (degrees)
	ElementMask els;//elements
	int         sg ;//space group number

	//string versions of data
	wxString    sEl;//element

	std::bitset<4> flg;//flags for external use (I'll use them to mark the filtered and checked status)

	//@brief: comparison op
	bool operator< (const MasterFile& rhs) const {return pth <  rhs.pth;}
	bool operator==(const MasterFile& rhs) const {return pth == rhs.pth;}

	//@brief: attempt to read the file currently stored in pth
	//@return: true if successful
	bool read();

	//@brief: convert various elements to strings nicely
	wxString system() const;
	wxString sgName() const {return sg < 1 || sg > 230 ? _("?") : SgNames[sg-1];}

	//@brief: constructor takes file path
	MasterFile(wxString path = "") : pth(path) {}

	//@brief: helper functions to sort a MasterFile by various fields
	//@note : most sort functions use pth as tiebreak
	typedef std::function<bool(const MasterFile&, const MasterFile&)> CompFunc;
	static bool SameComp(const CompFunc lhs, const CompFunc rhs) {return getAddress(lhs) == getAddress(rhs);}
	static bool PthAsc (const MasterFile& lhs, const MasterFile& rhs) {return lhs.pth < rhs.pth;}//operator<
	static bool PthDes (const MasterFile& lhs, const MasterFile& rhs) {return lhs.pth > rhs.pth;}
	static bool FrmAsc (const MasterFile& lhs, const MasterFile& rhs) {return lhs.frm == rhs.frm ? lhs.pth < rhs.pth : lhs.frm < rhs.frm;}
	static bool FrmDes (const MasterFile& lhs, const MasterFile& rhs) {return lhs.frm == rhs.frm ? lhs.pth > rhs.pth : lhs.frm > rhs.frm;}
	static bool NamAsc (const MasterFile& lhs, const MasterFile& rhs) {return lhs.nam == rhs.nam ? lhs.pth < rhs.pth : lhs.nam < rhs.nam;}
	static bool NamDes (const MasterFile& lhs, const MasterFile& rhs) {return lhs.nam == rhs.nam ? lhs.pth > rhs.pth : lhs.nam > rhs.nam;}
	static bool SybAsc (const MasterFile& lhs, const MasterFile& rhs) {return lhs.syb == rhs.syb ? lhs.pth < rhs.pth : lhs.syb < rhs.syb;}
	static bool SybDes (const MasterFile& lhs, const MasterFile& rhs) {return lhs.syb == rhs.syb ? lhs.pth > rhs.pth : lhs.syb > rhs.syb;}
	static bool KvAsc  (const MasterFile& lhs, const MasterFile& rhs) {return lhs.kv  == rhs.kv  ? lhs.pth < rhs.pth : lhs.kv  < rhs.kv ;}
	static bool KvDes  (const MasterFile& lhs, const MasterFile& rhs) {return lhs.kv  == rhs.kv  ? lhs.pth > rhs.pth : lhs.kv  > rhs.kv ;}
	static bool TltAsc (const MasterFile& lhs, const MasterFile& rhs) {return lhs.tlt == rhs.tlt ? lhs.pth < rhs.pth : lhs.tlt < rhs.tlt;}
	static bool TltDes (const MasterFile& lhs, const MasterFile& rhs) {return lhs.tlt == rhs.tlt ? lhs.pth > rhs.pth : lhs.tlt > rhs.tlt;}
	static bool ElsAsc (const MasterFile& lhs, const MasterFile& rhs) {return lhs.sEl == rhs.sEl ? lhs.pth < rhs.pth : lhs.sEl < rhs.sEl;}
	static bool ElsDes (const MasterFile& lhs, const MasterFile& rhs) {return lhs.sEl == rhs.sEl ? lhs.pth > rhs.pth : lhs.sEl > rhs.sEl;}
	static bool SgAsc  (const MasterFile& lhs, const MasterFile& rhs) {return lhs.sg  == rhs.sg  ? lhs.pth < rhs.pth : lhs.sg  < rhs.sg ;}
	static bool SgDes  (const MasterFile& lhs, const MasterFile& rhs) {return lhs.sg  == rhs.sg  ? lhs.pth > rhs.pth : lhs.sg  > rhs.sg ;}
	static bool SgNmAsc(const MasterFile& lhs, const MasterFile& rhs) {return (lhs.sg > 0 && lhs.sg < 231 && rhs.sg > 0 && rhs.sg < 231 && lhs.sg == rhs.sg) ? lhs.pth < rhs.pth : SgNames[lhs.sg-1] < SgNames[rhs.sg-1];}
	static bool SgNmDes(const MasterFile& lhs, const MasterFile& rhs) {return (lhs.sg > 0 && lhs.sg < 231 && rhs.sg > 0 && rhs.sg < 231 && lhs.sg == rhs.sg) ? lhs.pth > rhs.pth : SgNames[lhs.sg-1] > SgNames[rhs.sg-1];}
	static bool SysAsc (const MasterFile& lhs, const MasterFile& rhs) {return lhs.system() == rhs.system() ? ( lhs.sg == rhs.sg ? lhs.pth < rhs.pth : lhs.sg  < rhs.sg ) : lhs.system() < rhs.system();}
	static bool SysDes (const MasterFile& lhs, const MasterFile& rhs) {return lhs.system() == rhs.system() ? ( lhs.sg == rhs.sg ? lhs.pth > rhs.pth : lhs.sg  > rhs.sg ) : lhs.system() > rhs.system();}

	//@brief     : set the filter flag using bounds
	//@param fKv : upper and lower bounds for kv filtering
	//@param fTlt: upper and lower bounds for tilt filtering
	//@param fEl : bitmask of ~TODO should all or any elements be required?
	//@param fSg : upper and lower bounds for space group filtering
	void filter(std::pair<float, float> fKv, std::pair<float, float> fTlt, ElementMask fEl, std::pair<int  , int  > fSg) {
		flg.set(2, kv  >= fKv .first && kv  <= fKv .second &&
		           tlt >= fTlt.first && tlt <= fTlt.second && 
		           (fEl.none() ? true : (fEl & els) == fEl)&&
		           sg  >= fSg .first && sg  <= fSg .second);
	}

	//@brief    : set the search flag using search string
	//@param str: search string
	void search(wxString str) {flg.set(3, wxNOT_FOUND != pth.Lower().Find(str) || 
	                                      wxNOT_FOUND != frm.Lower().Find(str) ||
	                                      wxNOT_FOUND != nam.Lower().Find(str) ||
	                                      wxNOT_FOUND != syb.Lower().Find(str) );}

	//@brief : check if this item has passed the filter and search tests
	//@return: true if both passed, false otherwise
	bool display() const {return flg.test(2) && flg.test(3);}

	//@brief    : (un)check this item
	//@param chk: check state
	void check(bool chk) {flg.set(1, chk);}

	//@brief : get check state
	//@return: check state
	bool isChecked() const {return flg.test(1);}

	//@brief : check if the file has been parsed
	//@return: true if read() has previously returned true
	bool isRead() const {return flg.test(0);}
};

//@brief: gui element to display a a list of master files
class MasterFileList: public wxListCtrl{
	std::vector<MasterFile>                             m_allFiles;//list of all files (m_srtFnc exists)
	MasterFile::CompFunc                                m_srtFnc  ;//function to sort file list with (or nullptr for unsorted)
	std::vector< std::vector<MasterFile>::iterator >    m_disFiles;//files to display (filtered)

	std::pair<float, float> m_kvFilt   ;//kv filtering range
	std::pair<float, float> m_tltFilt  ;//tilt filtering range
	ElementMask             m_elFlt    ;//element filtering
	std::pair<int  , int  > m_sgFilt   ;//space group filtering

	wxString                m_strSrch  ;//current search string

	void populateDis();//@brief: build list of visible items

	void ColClicked   (wxListEvent& Event);//sort on column click
	void ItemChecked  (wxListEvent& Event) {m_disFiles[Event.GetIndex()]->check(true ); RefreshItem(Event.GetIndex());}// Event.Skip();}//track check/uncheck state
	void ItemUnChecked(wxListEvent& Event) {m_disFiles[Event.GetIndex()]->check(false); RefreshItem(Event.GetIndex());}// Event.Skip();}//track check/uncheck state

	public:
		//@brief : get the indices of currently selected items
		//@return: indices
		std::vector<long> GetSelection();
		
		//@brief    : set the sort function + update
		//@param func: new sort function (or nullptr to stop sorting)
		void setSort(MasterFile::CompFunc func);//set the sort function
		
		static std::pair<float, float> AllKv () {return std::pair<float, float>(0.0f   , 1.0e20 );}
		static std::pair<float, float> AllTlt() {return std::pair<float, float>(-400.0f, 400.0f );}
		static ElementMask             AllEl () {return ElementMask()                            ;}
		static std::pair<int  , int  > AllSg () {return std::pair<int  , int  >(1      , 230    );}

		//@brief    : set filter bounds
		//@param kv : upper and lower bounds for kv filtering
		//@param tlt: upper and lower bounds for tilt filtering
		//@param el : bitmask of ~TODO should all or any elements be required?
		//@param sg : upper and lower bounds for space group filtering
		void setFilters(std::pair<float, float> kv = AllKv(), std::pair<float, float> tlt = AllTlt(), ElementMask el = AllEl(), std::pair<int  , int  > sg = AllSg());

		//@brief    : get filter bounds
		//@param kv : location to write upper and lower bounds for kv filtering
		//@param tlt: location to write upper and lower bounds for tilt filtering
		//@param el : location to write bitmask of ~TODO should all or any elements be required?
		//@param sg : location to write upper and lower bounds for space group filtering
		void getFilters(std::pair<float, float>& kv, std::pair<float, float>& tlt, ElementMask& el, std::pair<int  , int  >& sg) const;

		//@brief    : set the string search term
		//@param str: new search term
		void setSearch(wxString str);

		//@brief     : add a new file to the list
		//@param item: file to add
		//@param chk : should the item be ticked
		//@note      : only added if item.isRead() ? true : item.read()
		//@return    : true if item was added (false otherwise)
		bool AddItem(MasterFile item, bool chk);

		//@brief    : check if the list has a master file (even if it isn't displayed)
		//@param pth: path of file to check for
		//@return   : true if the item is already included
		bool HasItem(wxString pth) {for(const MasterFile& f : m_allFiles) {if(f.pth == pth) return true;} return false;}

		//@brief     : remove an item
		//@param item: index of item to remove
		void RemoveItem(long item);

		MasterFileList(wxWindow* parent, const wxWindowID id = wxID_ANY, const wxPoint& pos = wxDefaultPosition, const wxSize& size = wxDefaultSize);
		~MasterFileList();

		wxString OnGetItemText(long item, long column) const;
		bool OnGetItemIsChecked(long item) const {return m_disFiles[item]->isChecked();}

		MasterFile GetItem(const size_t item) {return *m_disFiles[item];}

		//@brief: remove all currently selected items
		void RemoveSelected();

		//@brief   : move all currently selected items one position up or down
		//@param up: up or down
		void MoveSelected(const bool up);

		//@brief    : enable multiple item selection
		//@param enb: enable/disable
		void EnableMultipleSelection(const bool enb = true) {SetSingleStyle (wxLC_SINGLE_SEL, !enb);}

		//@brief : get all files in list
		//@return: all files (not just displayed ones)
		std::vector<wxString> getAllFiles() const;

};

////////////////////////////////////////////////////////////////////////
//                           MasterFileList                           //
////////////////////////////////////////////////////////////////////////

//@brief: build list of visible items
void MasterFileList::populateDis() {
	//rebuild display list
	m_disFiles.clear();
	SetItemCount(0);
	for(std::vector<MasterFile>::iterator iter = m_allFiles.begin(); iter != m_allFiles.end(); ++iter) {
		if(iter->display()) {
			m_disFiles.push_back(iter);
			RefreshItem(m_disFiles.size()-1);
		}
	}
	SetItemCount(m_disFiles.size());

	for(size_t i = 0; i < m_disFiles.size(); i++) {
		CheckItem(i, m_disFiles[i]->isChecked());
		// SetItemState(i, m_disFiles[i]->isChecked(), wxLIST_STATE_SELECTED);
		// RefreshItem(i);
	}
}

void MasterFileList::ColClicked(wxListEvent& Event) {
	if(m_srtFnc) {
		switch(Event.GetColumn()) {
			case 0: setSort( MasterFile::SameComp(m_srtFnc, MasterFile::PthAsc ) ? MasterFile::PthDes  : MasterFile::PthAsc ); break;
			case 1: setSort( MasterFile::SameComp(m_srtFnc, MasterFile::FrmAsc ) ? MasterFile::FrmDes  : MasterFile::FrmAsc ); break;
			case 2: setSort( MasterFile::SameComp(m_srtFnc, MasterFile::NamAsc ) ? MasterFile::NamDes  : MasterFile::NamAsc ); break;
			case 3: setSort( MasterFile::SameComp(m_srtFnc, MasterFile::SybAsc ) ? MasterFile::SybDes  : MasterFile::SybAsc ); break;
			case 4: setSort( MasterFile::SameComp(m_srtFnc, MasterFile::KvAsc  ) ? MasterFile::KvDes   : MasterFile::KvAsc  ); break;
			case 5: setSort( MasterFile::SameComp(m_srtFnc, MasterFile::TltAsc ) ? MasterFile::TltDes  : MasterFile::TltAsc ); break;
			case 6: setSort( MasterFile::SameComp(m_srtFnc, MasterFile::ElsAsc ) ? MasterFile::ElsDes  : MasterFile::ElsAsc ); break;
			case 7: setSort( MasterFile::SameComp(m_srtFnc, MasterFile::SysAsc ) ? MasterFile::SysDes  : MasterFile::SysAsc ); break;
			case 8: setSort( MasterFile::SameComp(m_srtFnc, MasterFile::SgAsc  ) ? MasterFile::SgDes   : MasterFile::SgAsc  ); break;
		}
	}
}

//@brief : get the indices of currently selected items
//@return: indices
std::vector<long> MasterFileList::GetSelection() {
	//build list of selected indices
	long idx = -1;
	std::vector<long> sel;
	while ((idx = GetNextItem(idx, wxLIST_NEXT_ALL, wxLIST_STATE_SELECTED)) != wxNOT_FOUND) sel.push_back(idx);
	return sel;
}

//@brief    : set the sort function + update
//@param func: new sort function (or nullptr to stop sorting)
void MasterFileList::setSort(MasterFile::CompFunc func) {
	bool same = func && m_srtFnc ? MasterFile::SameComp(func, m_srtFnc) : bool(func) == bool(m_srtFnc);
	if(same) return;
	m_srtFnc = func;
	if(m_srtFnc) {
		std::sort(m_allFiles.begin(), m_allFiles.end(), m_srtFnc);
		populateDis();
	}
}

void MasterFileList::setFilters(std::pair<float, float> kv, std::pair<float, float> tlt, ElementMask el, std::pair<int  , int  > sg) {
	if(kv != m_kvFilt || tlt != m_tltFilt || el != m_elFlt || sg != m_sgFilt ) {//filters actually changed
		//save filters and clear filtered list
		m_kvFilt  = kv ;
		m_tltFilt = tlt;
		m_elFlt   = el ;
		m_sgFilt  = sg ;

		//loop over all files applying new filters
		for(std::vector<MasterFile>::iterator iter = m_allFiles.begin(); iter != m_allFiles.end(); ++iter) iter->filter(m_kvFilt, m_tltFilt, m_elFlt, m_sgFilt);
		populateDis();
	}
}

//@brief    : get filter bounds
//@param kv : location to write upper and lower bounds for kv filtering
//@param tlt: location to write upper and lower bounds for tilt filtering
//@param el : location to write bitmask of ~TODO should all or any elements be required?
//@param sg : location to write upper and lower bounds for space group filtering
void MasterFileList::getFilters(std::pair<float, float>& kv, std::pair<float, float>& tlt, ElementMask& el, std::pair<int  , int  >& sg) const {
	kv  = m_kvFilt ;
	tlt = m_tltFilt;
	el  = m_elFlt  ;
	sg  = m_sgFilt ;
}

//@brief    : set the string search term
//@param str: new search term
void MasterFileList::setSearch(wxString str) {
	if(str.Lower() != m_strSrch) {
		m_strSrch = str.Lower();//case insensitive

		//loop over all files applying new filters
		for(std::vector<MasterFile>::iterator iter = m_allFiles.begin(); iter != m_allFiles.end(); ++iter) iter->search(m_strSrch);
		populateDis();
	}
}

//@brief     : add a new file to the list
//@param item: file to add
//@note      : only added if item.read() == true
//@param chk : should the item be ticked
//@return    : true if item was added (false otherwise)
bool MasterFileList::AddItem(MasterFile item, bool chk) {
	if(!(item.isRead() ? true : item.read())) return false;
	if(std::find(m_allFiles.cbegin(), m_allFiles.cend(), item) == m_allFiles.cend()) {
		//compute item filtering and add to list
		item.filter(m_kvFilt, m_tltFilt, m_elFlt, m_sgFilt);
		item.search(m_strSrch);
		item.check(chk);

		if(m_srtFnc) {
			std::vector<MasterFile>::iterator iter = std::upper_bound(m_allFiles.begin(), m_allFiles.end(), item, m_srtFnc);//find sorted insert location
			m_allFiles.insert(iter, item);//insert item into list
		} else {
			m_allFiles.push_back(item);
		}
		populateDis();//rebuild display list
		return true;
	}
	return false;
}

//@brief     : remove an item
//@param item: index of item to remove
void MasterFileList::RemoveItem(long item) {
	m_allFiles.erase(m_disFiles[item]);
	populateDis();
}

//Constructor, sets up virtual report list with 3 columns
MasterFileList::MasterFileList(wxWindow* parent, const wxWindowID id, const wxPoint& pos, const wxSize& size):
	wxListCtrl(parent, id, pos, size, wxLC_REPORT|wxLC_VIRTUAL|wxLC_HRULES|wxLC_SINGLE_SEL),
	m_srtFnc(MasterFile::PthAsc)
{
	// add columns
	AppendColumn("File"    , wxLIST_FORMAT_LEFT, 70);
	AppendColumn("Formula" , wxLIST_FORMAT_LEFT, 80);
	AppendColumn("Name"    , wxLIST_FORMAT_LEFT, 95);
	AppendColumn("S.Syb"   , wxLIST_FORMAT_LEFT, 50);
	AppendColumn("kV"      , wxLIST_FORMAT_LEFT, 35);
	AppendColumn("Tilt"    , wxLIST_FORMAT_LEFT, 35);
	AppendColumn("Els"     , wxLIST_FORMAT_LEFT, 45);
	AppendColumn("Laue"    , wxLIST_FORMAT_LEFT, 55);
	AppendColumn("SG #"    , wxLIST_FORMAT_LEFT, 40);
	// AppendColumn("SG Name" , wxLIST_FORMAT_LEFT, 70);
	EnableCheckBoxes (true);//used to be EnableCheckboxes
	SetItemCount(0);
	setFilters();//default filters
	EnableAlternateRowColours();

	// Connect Events
	this->Connect( wxEVT_COMMAND_LIST_COL_CLICK, wxListEventHandler( MasterFileList::ColClicked ), NULL, this );
	this->Connect( wxEVT_LIST_ITEM_UNCHECKED, wxListEventHandler( MasterFileList::ItemUnChecked ), NULL, this );
	this->Connect( wxEVT_LIST_ITEM_CHECKED, wxListEventHandler( MasterFileList::ItemChecked ), NULL, this );
}

MasterFileList::~MasterFileList()
{
	// Disconnect Events
	this->Disconnect( wxEVT_COMMAND_LIST_COL_CLICK, wxListEventHandler( MasterFileList::ColClicked ), NULL, this );
	this->Disconnect( wxEVT_LIST_ITEM_UNCHECKED, wxListEventHandler( MasterFileList::ItemUnChecked ), NULL, this );
	this->Disconnect( wxEVT_LIST_ITEM_CHECKED, wxListEventHandler( MasterFileList::ItemChecked ), NULL, this );
}

//Overload virtual method of wxListCtrl to provide text data for virtual list
wxString MasterFileList::OnGetItemText(long item, long column) const {
	const std::vector<MasterFile>::iterator& mf = m_disFiles[item];
	switch(column) {
		case 0: return mf->pth;
		case 1: return mf->frm;
		case 2: return mf->nam;
		case 3: return mf->syb;
		case 4: return wxString::Format(wxT("%.1f"), mf->kv );
		case 5: return wxString::Format(wxT("%.1f"), mf->tlt);
		case 6: return mf->sEl;
		case 7: return mf->system();
		case 8: return wxString::Format(wxT("%i"), mf->sg);
		// case 9: return mf->sgName();
		default: return _("?");
	}
	return _("?");
}

//@brief: remove all currently selected items
void MasterFileList::RemoveSelected() {
	//build list of selected indices and sort do we don't invalidate
	std::vector<long> sel = GetSelection();
	std::sort(sel.begin(), sel.end(), [this](const long& lhs, const long& rhs) {return std::distance(this->m_allFiles.begin(), this->m_disFiles[lhs]) < std::distance(this->m_allFiles.begin(), this->m_disFiles[rhs]);});

	//remove back to front to prevent invalidating unprocessed items and update list
	for(std::vector<long>::reverse_iterator iter = sel.rbegin(); iter != sel.rend(); ++iter) m_allFiles.erase(m_disFiles[*iter]);
	populateDis();
}

//@brief   : move all currently selected items one position up or down
//@param up: up or down
void MasterFileList::MoveSelected(const bool up) {
	//build list of selected indices
	std::vector<long> sel = GetSelection();
	if(!up) std::reverse(sel.begin(), sel.end());//order depends on movement direction
	bool haveMoved = false;//have we moved at least one item up
	for(size_t i = 0; i < sel.size(); i++) {
		bool move;
		if(up)
			move = sel[i] > 0;//don't move up if we're already at the top
		else
			move = m_disFiles.size() != sel[i] + 1;//don't move up if we're already at the top

		if(i > 0) move = haveMoved;//make sure there is space above this to move
		if(move) {
			haveMoved = true;
			SetItemState(sel[i], 0, wxLIST_STATE_SELECTED);
			if(up) {
				SetItemState(sel[i]-1, wxLIST_STATE_SELECTED, wxLIST_STATE_SELECTED);
				std::swap(m_disFiles[sel[i]-1], m_disFiles[sel[i]]);
			} else {
				SetItemState(sel[i]+1, wxLIST_STATE_SELECTED, wxLIST_STATE_SELECTED);
				std::swap(m_disFiles[sel[i]+1], m_disFiles[sel[i]]);
			}
		}
	}
	for(size_t i = 0; i < m_disFiles.size(); i++) RefreshItem(i);
}

//@brief : get all files in list
//@return: all files (not just displayed ones)
std::vector<wxString> MasterFileList::getAllFiles() const {
	std::vector<wxString> files;
	for(const MasterFile& mf : m_allFiles) files.push_back(mf.pth);
	return files;
}

////////////////////////////////////////////////////////////////////////
//                           MasterFileList                           //
////////////////////////////////////////////////////////////////////////

const wxString MasterFile::SgNames[230] = {//could do this with hm.hpp but this is easier and probably faster
	_("P 1"       ), _("P -1"      ), _("P 1 2 1"   ), _("P 1 21 1"  ), _("C 1 2 1"   ), _("P 1 m 1"   ), _("P 1 c 1"   ), _("C 1 m 1"   ),
	_("C 1 c 1"   ), _("P 1 2/m 1" ), _("P 1 21/m 1"), _("C 1 2/m 1" ), _("P 1 2/c 1" ), _("P 1 21/c 1"), _("C 1 2/c 1" ), _("P 2 2 2"   ),
	_("P 2 2 21"  ), _("P 21 21 2" ), _("P 21 21 21"), _("C 2 2 21"  ), _("C 2 2 2"   ), _("F 2 2 2"   ), _("I 2 2 2"   ), _("I 21 21 21"),
	_("P m m 2"   ), _("P m c 21"  ), _("P c c 2"   ), _("P m a 2"   ), _("P c a 21"  ), _("P n c 2"   ), _("P m n 21"  ), _("P b a 2"   ),
	_("P n a 21"  ), _("P n n 2"   ), _("C m m 2"   ), _("C m c 21"  ), _("C c c 2"   ), _("A m m 2"   ), _("A e m 2"   ), _("A m a 2"   ),
	_("A e a 2"   ), _("F m m 2"   ), _("F d d 2"   ), _("I m m 2"   ), _("I b a 2"   ), _("I m a 2"   ), _("P m m m"   ), _("P n n n"   ),
	_("P c c m"   ), _("P b a n"   ), _("P m m a"   ), _("P n n a"   ), _("P m n a"   ), _("P c c a"   ), _("P b a m"   ), _("P c c n"   ),
	_("P b c m"   ), _("P n n m"   ), _("P m m n"   ), _("P b c n"   ), _("P b c a"   ), _("P n m a"   ), _("C m c m"   ), _("C m c e"   ),
	_("C m m m"   ), _("C c c m"   ), _("C m m e"   ), _("C c c e"   ), _("F m m m"   ), _("F d d d"   ), _("I m m m"   ), _("I b a m"   ),
	_("I b c a"   ), _("I m m a"   ), _("P 4"       ), _("P 41"      ), _("P 42"      ), _("P 43"      ), _("I 4"       ), _("I 41"      ),
	_("P -4"      ), _("I -4"      ), _("P 4/m"     ), _("P 42/m"    ), _("P 4/n"     ), _("P 42/n"    ), _("I 4/m"     ), _("I 41/a"    ),
	_("P 4 2 2"   ), _("P 4 21 2"  ), _("P 41 2 2"  ), _("P 41 21 2" ), _("P 42 2 2"  ), _("P 42 21 2" ), _("P 43 2 2"  ), _("P 43 21 2" ),
	_("I 4 2 2"   ), _("I 41 2 2"  ), _("P 4 m m"   ), _("P 4 b m"   ), _("P 42 c m"  ), _("P 42 n m"  ), _("P 4 c c"   ), _("P 4 n c"   ),
	_("P 42 m c"  ), _("P 42 b c"  ), _("I 4 m m"   ), _("I 4 c m"   ), _("I 41 m d"  ), _("I 41 c d"  ), _("P -4 2 m"  ), _("P -4 2 c"  ),
	_("P -4 21 m" ), _("P -4 21 c" ), _("P -4 m 2"  ), _("P -4 c 2"  ), _("P -4 b 2"  ), _("P -4 n 2"  ), _("I -4 m 2"  ), _("I -4 c 2"  ),
	_("I -4 2 m"  ), _("I -4 2 d"  ), _("P 4/m m m" ), _("P 4/m c c" ), _("P 4/n b m" ), _("P 4/n n c" ), _("P 4/m b m" ), _("P 4/m n c" ),
	_("P 4/n m m" ), _("P 4/n c c" ), _("P 42/m m c"), _("P 42/m c m"), _("P 42/n b c"), _("P 42/n n m"), _("P 42/m b c"), _("P 42/m n m"),
	_("P 42/n m c"), _("P 42/n c m"), _("I 4/m m m" ), _("I 4/m c m" ), _("I 41/a m d"), _("I 41/a c d"), _("P 3"       ), _("P 31"      ),
	_("P 32"      ), _("R 3"       ), _("P -3"      ), _("R -3"      ), _("P 3 1 2"   ), _("P 3 2 1"   ), _("P 31 1 2"  ), _("P 31 2 1"  ),
	_("P 32 1 2"  ), _("P 32 2 1"  ), _("R 3 2"     ), _("P 3 m 1"   ), _("P 3 1 m"   ), _("P 3 c 1"   ), _("P 3 1 c"   ), _("R 3 m"     ),
	_("R 3 c"     ), _("P -3 1 m"  ), _("P -3 1 c"  ), _("P -3 m 1"  ), _("P -3 c 1"  ), _("R -3 m"    ), _("R -3 c"    ), _("P 6"       ),
	_("P 61"      ), _("P 65"      ), _("P 62"      ), _("P 64"      ), _("P 63"      ), _("P -6"      ), _("P 6/m"     ), _("P 63/m"    ),
	_("P 6 2 2"   ), _("P 61 2 2"  ), _("P 65 2 2"  ), _("P 62 2 2"  ), _("P 64 2 2"  ), _("P 63 2 2"  ), _("P 6 m m"   ), _("P 6 c c"   ),
	_("P 63 c m"  ), _("P 63 m c"  ), _("P -6 m 2"  ), _("P -6 c 2"  ), _("P -6 2 m"  ), _("P -6 2 c"  ), _("P 6/m m m" ), _("P 6/m c c" ),
	_("P 63/m c m"), _("P 63/m m c"), _("P 2 3"     ), _("F 2 3"     ), _("I 2 3"     ), _("P 21 3"    ), _("I 21 3"    ), _("P m 3"     ),
	_("P n 3"     ), _("F m 3"     ), _("F d 3"     ), _("I m 3"     ), _("P a 3"     ), _("I a 3"     ), _("P 4 3 2"   ), _("P 42 3 2"  ),
	_("F 4 3 2"   ), _("F 41 3 2"  ), _("I 4 3 2"   ), _("P 43 3 2"  ), _("P 41 3 2"  ), _("I 41 3 2"  ), _("P -4 3 m"  ), _("F -4 3 m"  ),
	_("I -4 3 m"  ), _("P -4 3 n"  ), _("F -4 3 c"  ), _("I -4 3 d"  ), _("P m 3 m"   ), _("P n 3 n"   ), _("P m 3 n"   ), _("P n 3 m"   ),
	_("F m 3 m"   ), _("F m 3 c"   ), _("F d 3 m"   ), _("F d 3 c"   ), _("I m 3 m"   ), _("I a 3 d"   ), 
};

wxString MasterFile::system() const {
	if     (sg < 1  ) return _("?"    );
	else if(sg < 3  ) return _("-1"   );
	else if(sg < 16 ) return _("2/m"  );
	else if(sg < 75 ) return _("mmm"  );
	else if(sg < 89 ) return _("4/m"  );
	else if(sg < 143) return _("4/mmm");
	else if(sg < 149) return _("-3"   );
	else if(sg < 168) return _("3m"   );
	else if(sg < 177) return _("6/m"  );
	else if(sg < 195) return _("6/mmm");
	else if(sg < 207) return _("m-3"  );
	else if(sg < 231) return _("m-3m" );
	// else if(sg < 3  ) return _("Tric");//"Triclinic"   
	// else if(sg < 16 ) return _("Mono");//"Monoclinic"  
	// else if(sg < 75 ) return _("Orth");//"Orthorhombic"
	// else if(sg < 143) return _("Tet" );//"Tetragonal"  
	// else if(sg < 168) return _("Trig");//"Trigonal"    
	// else if(sg < 195) return _("Hex" );//"Hexagonal"   
	// else if(sg < 231) return _("Cub" );//"Cubic"       
	else return _("?");
}

#include "sht_file.hpp"
#include <fstream>

//@brief: attempt to read the file currently stored in pth
//@return: true if successful
bool MasterFile::read() {
	frm = "";
	nam = "";
	syb = "";
	kv  = 0;
	tlt = 0;
	els = ElementMask();
	sg  = 0;
	flg.set(0, false);

	if(wxFileExists(pth)) {//first check if file exists
		if("sht" == wxFileName(pth).GetExt().Lower()) {//make sure the file extension is correct
			try {
				//read file and make sure it is an EBSD master pattern spectra
				sht::File file;
				std::ifstream is(pth.ToStdString(), std::ios::in | std::ios::binary);
				file.read(is);
				if(sht::Modality::EBSD != file.header.modality()) return false;//some other type of master pattern

				//save meta data
				kv  = file.header.beamEnergy();
				tlt = file.header.primaryAngle();
				sg  = file.mpData.sgEff();

				//build up bitmask of elements
				for(const sht::CrystalData& xtal : file.mpData.xtals) {
					std::string xFrm = xtal.form.substr(0, xtal.formulaLen  ());
					std::string xNam = xtal.name.substr(0, xtal.matNameLen  ());
					std::string xSyb = xtal.symb.substr(0, xtal.structSymLen());
					if(!xFrm.empty()) {
						if(!frm.IsEmpty()) frm += '+';
						frm += xFrm;
					}
					if(!xNam.empty()) {
						if(!nam.IsEmpty()) nam += '+';
						nam += xNam;
					}
					if(!xSyb.empty()) {
						if(!syb.IsEmpty()) syb += '+';
						syb += xSyb;
					}

					for(const sht::AtomData& at : xtal.atoms) {
						size_t z = at.atZ();
						if(z > 0 && z < 119) {
							els.set(z, true);
						} else {
							return false;
						}
					}
				}

				//convert element bitmask to string and return
				sEl = els.str();
				flg.set(0, true);
				return true;//read everything
			} catch (...) {
				return false;
			}
		}
	}
	return false;
}

#endif//_MASTERFILELIST_H_

