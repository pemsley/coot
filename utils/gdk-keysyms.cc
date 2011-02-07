/* GDK - The GIMP Drawing Kit
 * Copyright (C) 1995-1997 Peter Mattis, Spencer Kimball and Josh MacDonald
 * Copyright (C) 2005 GNOME Foundation
 * Copyright 2008 by The University of Oxford
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 3 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */

#include "coot-utils.hh"
namespace coot {
   namespace util {
      std::vector<std::pair<std::string, int> > key_sym_vec_1();
      std::vector<std::pair<std::string, int> > key_sym_vec_2();
      std::vector<std::pair<std::string, int> > key_sym_vec_3();
      std::vector<std::pair<std::string, int> > key_sym_vec_4();
      std::vector<std::pair<std::string, int> > key_sym_vec_5();
   }
} 

// Return -1 on failure to find symbol
int 
coot::util::decode_keysym(const std::string &s) {

   int r = -1;
   std::vector<std::pair<std::string, int> > v = coot::util::key_sym_vec();
   for (unsigned int i=0; i<v.size(); i++) {
      if (s == v[i].first) { 
	 r = v[i].second;
	 break;
      }
   }
   return r;
}


std::vector<std::pair<std::string, int> >
coot::util::key_sym_vec_1() {

   std::vector<std::pair<std::string, int> > a; 
   a.push_back(std::pair<std::string,int>("BackSpace", 65288));
   a.push_back(std::pair<std::string,int>("Tab", 65289));
   a.push_back(std::pair<std::string,int>("Linefeed", 65290));
   a.push_back(std::pair<std::string,int>("Clear", 65291));
   a.push_back(std::pair<std::string,int>("Return", 65293));
   a.push_back(std::pair<std::string,int>("Pause", 65299));
   a.push_back(std::pair<std::string,int>("Scroll_Lock", 65300));
   a.push_back(std::pair<std::string,int>("Sys_Req", 65301));
   a.push_back(std::pair<std::string,int>("Escape", 65307));
   a.push_back(std::pair<std::string,int>("Delete", 65535));
   a.push_back(std::pair<std::string,int>("Multi_key", 65312));
   a.push_back(std::pair<std::string,int>("Home", 65360));
   a.push_back(std::pair<std::string,int>("Left", 65361));
   a.push_back(std::pair<std::string,int>("Up", 65362));
   a.push_back(std::pair<std::string,int>("Right", 65363));
   a.push_back(std::pair<std::string,int>("Down", 65364));
   a.push_back(std::pair<std::string,int>("Prior", 65365));
   a.push_back(std::pair<std::string,int>("Page_Up", 65365));
   a.push_back(std::pair<std::string,int>("Next", 65366));
   a.push_back(std::pair<std::string,int>("Page_Down", 65366));
   a.push_back(std::pair<std::string,int>("End", 65367));
   a.push_back(std::pair<std::string,int>("Begin", 65368));
   a.push_back(std::pair<std::string,int>("Select", 65376));
   a.push_back(std::pair<std::string,int>("Print", 65377));
   a.push_back(std::pair<std::string,int>("Execute", 65378));
   a.push_back(std::pair<std::string,int>("Insert", 65379));
   a.push_back(std::pair<std::string,int>("Undo", 65381));
   a.push_back(std::pair<std::string,int>("Redo", 65382));
   a.push_back(std::pair<std::string,int>("Menu", 65383));
   a.push_back(std::pair<std::string,int>("Find", 65384));
   a.push_back(std::pair<std::string,int>("Cancel", 65385));
   a.push_back(std::pair<std::string,int>("Help", 65386));
   a.push_back(std::pair<std::string,int>("Break", 65387));
   a.push_back(std::pair<std::string,int>("Mode_switch", 65406));
   a.push_back(std::pair<std::string,int>("script_switch", 65406));
   a.push_back(std::pair<std::string,int>("Num_Lock", 65407));
   a.push_back(std::pair<std::string,int>("KP_Space", 65408));
   a.push_back(std::pair<std::string,int>("KP_Tab", 65417));
   a.push_back(std::pair<std::string,int>("KP_Enter", 65421));
   a.push_back(std::pair<std::string,int>("KP_F1", 65425));
   a.push_back(std::pair<std::string,int>("KP_F2", 65426));
   a.push_back(std::pair<std::string,int>("KP_F3", 65427));
   a.push_back(std::pair<std::string,int>("KP_F4", 65428));
   a.push_back(std::pair<std::string,int>("KP_Home", 65429));
   a.push_back(std::pair<std::string,int>("KP_Left", 65430));
   a.push_back(std::pair<std::string,int>("KP_Up", 65431));
   a.push_back(std::pair<std::string,int>("KP_Right", 65432));
   a.push_back(std::pair<std::string,int>("KP_Down", 65433));
   a.push_back(std::pair<std::string,int>("KP_Prior", 65434));
   a.push_back(std::pair<std::string,int>("KP_Page_Up", 65434));
   a.push_back(std::pair<std::string,int>("KP_Next", 65435));
   a.push_back(std::pair<std::string,int>("KP_Page_Down", 65435));
   a.push_back(std::pair<std::string,int>("KP_End", 65436));
   a.push_back(std::pair<std::string,int>("KP_Begin", 65437));
   a.push_back(std::pair<std::string,int>("KP_Insert", 65438));
   a.push_back(std::pair<std::string,int>("KP_Delete", 65439));
   a.push_back(std::pair<std::string,int>("KP_Equal", 65469));
   a.push_back(std::pair<std::string,int>("KP_Multiply", 65450));
   a.push_back(std::pair<std::string,int>("KP_Add", 65451));
   a.push_back(std::pair<std::string,int>("KP_Separator", 65452));
   a.push_back(std::pair<std::string,int>("KP_Subtract", 65453));
   a.push_back(std::pair<std::string,int>("KP_Decimal", 65454));
   a.push_back(std::pair<std::string,int>("KP_Divide", 65455));
   a.push_back(std::pair<std::string,int>("KP_0", 65456));
   a.push_back(std::pair<std::string,int>("KP_1", 65457));
   a.push_back(std::pair<std::string,int>("KP_2", 65458));
   a.push_back(std::pair<std::string,int>("KP_3", 65459));
   a.push_back(std::pair<std::string,int>("KP_4", 65460));
   a.push_back(std::pair<std::string,int>("KP_5", 65461));
   a.push_back(std::pair<std::string,int>("KP_6", 65462));
   a.push_back(std::pair<std::string,int>("KP_7", 65463));
   a.push_back(std::pair<std::string,int>("KP_8", 65464));
   a.push_back(std::pair<std::string,int>("KP_9", 65465));
   return a;
} 


std::vector<std::pair<std::string, int> >
coot::util::key_sym_vec_2() {

   std::vector<std::pair<std::string, int> > a;
   {
   a.push_back(std::pair<std::string,int>("F1", 65470));
   a.push_back(std::pair<std::string,int>("F2", 65471));
   a.push_back(std::pair<std::string,int>("F3", 65472));
   a.push_back(std::pair<std::string,int>("F4", 65473));
   a.push_back(std::pair<std::string,int>("F5", 65474));
   a.push_back(std::pair<std::string,int>("F6", 65475));
   a.push_back(std::pair<std::string,int>("F7", 65476));
   a.push_back(std::pair<std::string,int>("F8", 65477));
   a.push_back(std::pair<std::string,int>("F9", 65478));
   a.push_back(std::pair<std::string,int>("F10", 65479));
   a.push_back(std::pair<std::string,int>("F11", 65480));
   a.push_back(std::pair<std::string,int>("L1", 65480));
   a.push_back(std::pair<std::string,int>("F12", 65481));
   a.push_back(std::pair<std::string,int>("L2", 65481));
   a.push_back(std::pair<std::string,int>("F13", 65482));
   a.push_back(std::pair<std::string,int>("L3", 65482));
   a.push_back(std::pair<std::string,int>("F14", 65483));
   a.push_back(std::pair<std::string,int>("L4", 65483));
   a.push_back(std::pair<std::string,int>("F15", 65484));
   a.push_back(std::pair<std::string,int>("L5", 65484));
   a.push_back(std::pair<std::string,int>("F16", 65485));
   a.push_back(std::pair<std::string,int>("L6", 65485));
   a.push_back(std::pair<std::string,int>("F17", 65486));
   a.push_back(std::pair<std::string,int>("L7", 65486));
   a.push_back(std::pair<std::string,int>("F18", 65487));
   a.push_back(std::pair<std::string,int>("L8", 65487));
   a.push_back(std::pair<std::string,int>("F19", 65488));
   a.push_back(std::pair<std::string,int>("L9", 65488));
   a.push_back(std::pair<std::string,int>("F20", 65489));
   a.push_back(std::pair<std::string,int>("L10", 65489));
   a.push_back(std::pair<std::string,int>("F21", 65490));
   a.push_back(std::pair<std::string,int>("R1", 65490));
   a.push_back(std::pair<std::string,int>("F22", 65491));
   a.push_back(std::pair<std::string,int>("R2", 65491));
   a.push_back(std::pair<std::string,int>("F23", 65492));
   a.push_back(std::pair<std::string,int>("R3", 65492));
   a.push_back(std::pair<std::string,int>("F24", 65493));
   a.push_back(std::pair<std::string,int>("R4", 65493));
   a.push_back(std::pair<std::string,int>("F25", 65494));
   a.push_back(std::pair<std::string,int>("R5", 65494));
   a.push_back(std::pair<std::string,int>("F26", 65495));
   a.push_back(std::pair<std::string,int>("R6", 65495));
   a.push_back(std::pair<std::string,int>("F27", 65496));
   a.push_back(std::pair<std::string,int>("R7", 65496));
   a.push_back(std::pair<std::string,int>("F28", 65497));
   a.push_back(std::pair<std::string,int>("R8", 65497));
   a.push_back(std::pair<std::string,int>("F29", 65498));
   a.push_back(std::pair<std::string,int>("R9", 65498));
   a.push_back(std::pair<std::string,int>("F30", 65499));
   a.push_back(std::pair<std::string,int>("R10", 65499));
   a.push_back(std::pair<std::string,int>("F31", 65500));
   a.push_back(std::pair<std::string,int>("R11", 65500));
   a.push_back(std::pair<std::string,int>("F32", 65501));
   a.push_back(std::pair<std::string,int>("R12", 65501));
   a.push_back(std::pair<std::string,int>("F33", 65502));
   a.push_back(std::pair<std::string,int>("R13", 65502));
   a.push_back(std::pair<std::string,int>("F34", 65503));
   a.push_back(std::pair<std::string,int>("R14", 65503));
   a.push_back(std::pair<std::string,int>("F35", 65504));
   a.push_back(std::pair<std::string,int>("R15", 65504));
   a.push_back(std::pair<std::string,int>("Shift_L", 65505));
   a.push_back(std::pair<std::string,int>("Shift_R", 65506));
   a.push_back(std::pair<std::string,int>("Control_L", 65507));
   a.push_back(std::pair<std::string,int>("Control_R", 65508));
   a.push_back(std::pair<std::string,int>("Caps_Lock", 65509));
   a.push_back(std::pair<std::string,int>("Shift_Lock", 65510));
   a.push_back(std::pair<std::string,int>("Meta_L", 65511));
   a.push_back(std::pair<std::string,int>("Meta_R", 65512));
   a.push_back(std::pair<std::string,int>("Alt_L", 65513));
   a.push_back(std::pair<std::string,int>("Alt_R", 65514));
   a.push_back(std::pair<std::string,int>("Super_L", 65515));
   a.push_back(std::pair<std::string,int>("Super_R", 65516));
   a.push_back(std::pair<std::string,int>("Hyper_L", 65517));
   a.push_back(std::pair<std::string,int>("Hyper_R", 65518));
   }
   return a;
}

std::vector<std::pair<std::string, int> >
coot::util::key_sym_vec_3() {

   std::vector<std::pair<std::string, int> > a;
   { 
   a.push_back(std::pair<std::string,int>("ISO_Lock", 65025));
   a.push_back(std::pair<std::string,int>("ISO_Level2_Latch", 65026));
   a.push_back(std::pair<std::string,int>("ISO_Level3_Shift", 65027));
   a.push_back(std::pair<std::string,int>("ISO_Level3_Latch", 65028));
   a.push_back(std::pair<std::string,int>("ISO_Level3_Lock", 65029));
   a.push_back(std::pair<std::string,int>("ISO_Group_Shift", 65406));
   a.push_back(std::pair<std::string,int>("ISO_Group_Latch", 65030));
   a.push_back(std::pair<std::string,int>("ISO_Group_Lock", 65031));
   a.push_back(std::pair<std::string,int>("ISO_Next_Group", 65032));
   a.push_back(std::pair<std::string,int>("ISO_Next_Group_Lock", 65033));
   a.push_back(std::pair<std::string,int>("ISO_Prev_Group", 65034));
   a.push_back(std::pair<std::string,int>("ISO_Prev_Group_Lock", 65035));
   a.push_back(std::pair<std::string,int>("ISO_First_Group", 65036));
   a.push_back(std::pair<std::string,int>("ISO_First_Group_Lock", 65037));
   a.push_back(std::pair<std::string,int>("ISO_Last_Group", 65038));
   a.push_back(std::pair<std::string,int>("ISO_Last_Group_Lock", 65039));
   a.push_back(std::pair<std::string,int>("ISO_Left_Tab", 65056));
   a.push_back(std::pair<std::string,int>("ISO_Move_Line_Up", 65057));
   a.push_back(std::pair<std::string,int>("ISO_Move_Line_Down", 65058));
   a.push_back(std::pair<std::string,int>("ISO_Partial_Line_Up", 65059));
   a.push_back(std::pair<std::string,int>("ISO_Partial_Line_Down", 65060));
   a.push_back(std::pair<std::string,int>("ISO_Partial_Space_Left", 65061));
   a.push_back(std::pair<std::string,int>("ISO_Partial_Space_Right", 65062));
   a.push_back(std::pair<std::string,int>("ISO_Set_Margin_Left", 65063));
   a.push_back(std::pair<std::string,int>("ISO_Set_Margin_Right", 65064));
   a.push_back(std::pair<std::string,int>("ISO_Release_Margin_Left", 65065));
   a.push_back(std::pair<std::string,int>("ISO_Release_Margin_Right", 65066));
   a.push_back(std::pair<std::string,int>("ISO_Release_Both_Margins", 65067));
   a.push_back(std::pair<std::string,int>("ISO_Fast_Cursor_Left", 65068));
   a.push_back(std::pair<std::string,int>("ISO_Fast_Cursor_Right", 65069));
   a.push_back(std::pair<std::string,int>("ISO_Fast_Cursor_Up", 65070));
   a.push_back(std::pair<std::string,int>("ISO_Fast_Cursor_Down", 65071));
   a.push_back(std::pair<std::string,int>("ISO_Continuous_Underline", 65072));
   a.push_back(std::pair<std::string,int>("ISO_Discontinuous_Underline", 65073));
   a.push_back(std::pair<std::string,int>("ISO_Emphasize", 65074));
   a.push_back(std::pair<std::string,int>("ISO_Center_Object", 65075));
   a.push_back(std::pair<std::string,int>("ISO_Enter", 65076));
   a.push_back(std::pair<std::string,int>("dead_grave", 65104));
   a.push_back(std::pair<std::string,int>("dead_acute", 65105));
   a.push_back(std::pair<std::string,int>("dead_circumflex", 65106));
   a.push_back(std::pair<std::string,int>("dead_tilde", 65107));
   a.push_back(std::pair<std::string,int>("dead_macron", 65108));
   a.push_back(std::pair<std::string,int>("dead_breve", 65109));
   a.push_back(std::pair<std::string,int>("dead_abovedot", 65110));
   a.push_back(std::pair<std::string,int>("dead_diaeresis", 65111));
   a.push_back(std::pair<std::string,int>("dead_abovering", 65112));
   a.push_back(std::pair<std::string,int>("dead_doubleacute", 65113));
   a.push_back(std::pair<std::string,int>("dead_caron", 65114));
   a.push_back(std::pair<std::string,int>("dead_cedilla", 65115));
   a.push_back(std::pair<std::string,int>("dead_ogonek", 65116));
   a.push_back(std::pair<std::string,int>("dead_iota", 65117));
   a.push_back(std::pair<std::string,int>("dead_voiced_sound", 65118));
   a.push_back(std::pair<std::string,int>("dead_semivoiced_sound", 65119));
   a.push_back(std::pair<std::string,int>("dead_belowdot", 65120));
   a.push_back(std::pair<std::string,int>("dead_hook", 65121));
   a.push_back(std::pair<std::string,int>("dead_horn", 65122));
   a.push_back(std::pair<std::string,int>("First_Virtual_Screen", 65232));
   a.push_back(std::pair<std::string,int>("Prev_Virtual_Screen", 65233));
   a.push_back(std::pair<std::string,int>("Next_Virtual_Screen", 65234));
   a.push_back(std::pair<std::string,int>("Last_Virtual_Screen", 65236));
   a.push_back(std::pair<std::string,int>("Terminate_Server", 65237));
   a.push_back(std::pair<std::string,int>("AccessX_Enable", 65136));
   a.push_back(std::pair<std::string,int>("AccessX_Feedback_Enable", 65137));
   a.push_back(std::pair<std::string,int>("RepeatKeys_Enable", 65138));
   a.push_back(std::pair<std::string,int>("SlowKeys_Enable", 65139));
   a.push_back(std::pair<std::string,int>("BounceKeys_Enable", 65140));
   a.push_back(std::pair<std::string,int>("StickyKeys_Enable", 65141));
   a.push_back(std::pair<std::string,int>("MouseKeys_Enable", 65142));
   a.push_back(std::pair<std::string,int>("MouseKeys_Accel_Enable", 65143));
   a.push_back(std::pair<std::string,int>("Overlay1_Enable", 65144));
   a.push_back(std::pair<std::string,int>("Overlay2_Enable", 65145));
   a.push_back(std::pair<std::string,int>("AudibleBell_Enable", 65146));
   a.push_back(std::pair<std::string,int>("Pointer_Left", 65248));
   a.push_back(std::pair<std::string,int>("Pointer_Right", 65249));
   a.push_back(std::pair<std::string,int>("Pointer_Up", 65250));
   a.push_back(std::pair<std::string,int>("Pointer_Down", 65251));
   a.push_back(std::pair<std::string,int>("Pointer_UpLeft", 65252));
   a.push_back(std::pair<std::string,int>("Pointer_UpRight", 65253));
   a.push_back(std::pair<std::string,int>("Pointer_DownLeft", 65254));
   a.push_back(std::pair<std::string,int>("Pointer_DownRight", 65255));
   a.push_back(std::pair<std::string,int>("Pointer_Button_Dflt", 65256));
   a.push_back(std::pair<std::string,int>("Pointer_Button1", 65257));
   a.push_back(std::pair<std::string,int>("Pointer_Button2", 65258));
   a.push_back(std::pair<std::string,int>("Pointer_Button3", 65259));
   a.push_back(std::pair<std::string,int>("Pointer_Button4", 65260));
   a.push_back(std::pair<std::string,int>("Pointer_Button5", 65261));
   a.push_back(std::pair<std::string,int>("Pointer_DblClick_Dflt", 65262));
   a.push_back(std::pair<std::string,int>("Pointer_DblClick1", 65263));
   a.push_back(std::pair<std::string,int>("Pointer_DblClick2", 65264));
   a.push_back(std::pair<std::string,int>("Pointer_DblClick3", 65265));
   a.push_back(std::pair<std::string,int>("Pointer_DblClick4", 65266));
   a.push_back(std::pair<std::string,int>("Pointer_DblClick5", 65267));
   a.push_back(std::pair<std::string,int>("Pointer_Drag_Dflt", 65268));
   a.push_back(std::pair<std::string,int>("Pointer_Drag1", 65269));
   a.push_back(std::pair<std::string,int>("Pointer_Drag2", 65270));
   a.push_back(std::pair<std::string,int>("Pointer_Drag3", 65271));
   a.push_back(std::pair<std::string,int>("Pointer_Drag4", 65272));
   a.push_back(std::pair<std::string,int>("Pointer_Drag5", 65277));
   a.push_back(std::pair<std::string,int>("Pointer_EnableKeys", 65273));
   a.push_back(std::pair<std::string,int>("Pointer_Accelerate", 65274));
   a.push_back(std::pair<std::string,int>("Pointer_DfltBtnNext", 65275));
   a.push_back(std::pair<std::string,int>("Pointer_DfltBtnPrev", 65276));
   }
   return a;
}


std::vector<std::pair<std::string, int> >
coot::util::key_sym_vec_4() {


   std::vector<std::pair<std::string, int> > a;

   a.push_back(std::pair<std::string,int>("space", 32));
   a.push_back(std::pair<std::string,int>("exclam", 33));
   a.push_back(std::pair<std::string,int>("quotedbl", 34));
   a.push_back(std::pair<std::string,int>("numbersign", 35));
   a.push_back(std::pair<std::string,int>("dollar", 36));
   a.push_back(std::pair<std::string,int>("percent", 37));
   a.push_back(std::pair<std::string,int>("ampersand", 38));
   a.push_back(std::pair<std::string,int>("apostrophe", 39));
   a.push_back(std::pair<std::string,int>("quoteright", 39));
   a.push_back(std::pair<std::string,int>("parenleft", 40));
   a.push_back(std::pair<std::string,int>("parenright", 41));
   a.push_back(std::pair<std::string,int>("asterisk", 42));
   a.push_back(std::pair<std::string,int>("plus", 43));
   a.push_back(std::pair<std::string,int>("comma", 44));
   a.push_back(std::pair<std::string,int>("minus", 45));
   a.push_back(std::pair<std::string,int>("period", 46));
   a.push_back(std::pair<std::string,int>("slash", 47));
   a.push_back(std::pair<std::string,int>("0", 48));
   a.push_back(std::pair<std::string,int>("1", 49));
   a.push_back(std::pair<std::string,int>("2", 50));
   a.push_back(std::pair<std::string,int>("3", 51));
   a.push_back(std::pair<std::string,int>("4", 52));
   a.push_back(std::pair<std::string,int>("5", 53));
   a.push_back(std::pair<std::string,int>("6", 54));
   a.push_back(std::pair<std::string,int>("7", 55));
   a.push_back(std::pair<std::string,int>("8", 56));
   a.push_back(std::pair<std::string,int>("9", 57));
   a.push_back(std::pair<std::string,int>("colon", 58));
   a.push_back(std::pair<std::string,int>("semicolon", 59));
   a.push_back(std::pair<std::string,int>("less", 60));
   a.push_back(std::pair<std::string,int>("equal", 61));
   a.push_back(std::pair<std::string,int>("greater", 62));
   a.push_back(std::pair<std::string,int>("question", 63));
   a.push_back(std::pair<std::string,int>("at", 64));
   a.push_back(std::pair<std::string,int>("A", 65));
   a.push_back(std::pair<std::string,int>("B", 66));
   a.push_back(std::pair<std::string,int>("C", 67));
   a.push_back(std::pair<std::string,int>("D", 68));
   a.push_back(std::pair<std::string,int>("E", 69));
   a.push_back(std::pair<std::string,int>("F", 70));
   a.push_back(std::pair<std::string,int>("G", 71));
   a.push_back(std::pair<std::string,int>("H", 72));
   a.push_back(std::pair<std::string,int>("I", 73));
   a.push_back(std::pair<std::string,int>("J", 74));
   a.push_back(std::pair<std::string,int>("K", 75));
   a.push_back(std::pair<std::string,int>("L", 76));
   a.push_back(std::pair<std::string,int>("M", 77));
   a.push_back(std::pair<std::string,int>("N", 78));
   a.push_back(std::pair<std::string,int>("O", 79));
   a.push_back(std::pair<std::string,int>("P", 80));
   a.push_back(std::pair<std::string,int>("Q", 81));
   a.push_back(std::pair<std::string,int>("R", 82));
   a.push_back(std::pair<std::string,int>("S", 83));
   a.push_back(std::pair<std::string,int>("T", 84));
   a.push_back(std::pair<std::string,int>("U", 85));
   a.push_back(std::pair<std::string,int>("V", 86));
   a.push_back(std::pair<std::string,int>("W", 87));
   a.push_back(std::pair<std::string,int>("X", 88));
   a.push_back(std::pair<std::string,int>("Y", 89));
   a.push_back(std::pair<std::string,int>("Z", 90));
   a.push_back(std::pair<std::string,int>("bracketleft", 91)); // not paren
   a.push_back(std::pair<std::string,int>("backslash", 92));
   a.push_back(std::pair<std::string,int>("bracketright", 93));
   a.push_back(std::pair<std::string,int>("asciicircum", 94));
   a.push_back(std::pair<std::string,int>("underscore", 95));
   a.push_back(std::pair<std::string,int>("grave", 96));
   a.push_back(std::pair<std::string,int>("quoteleft", 96)); // same as above
   a.push_back(std::pair<std::string,int>("a", 97));
   a.push_back(std::pair<std::string,int>("b", 98));
   a.push_back(std::pair<std::string,int>("c", 99));
   a.push_back(std::pair<std::string,int>("d", 100));
   a.push_back(std::pair<std::string,int>("e", 101));
   a.push_back(std::pair<std::string,int>("f", 102));
   a.push_back(std::pair<std::string,int>("g", 103));
   a.push_back(std::pair<std::string,int>("h", 104));
   a.push_back(std::pair<std::string,int>("i", 105));
   a.push_back(std::pair<std::string,int>("j", 106));
   a.push_back(std::pair<std::string,int>("k", 107));
   a.push_back(std::pair<std::string,int>("l", 108));
   a.push_back(std::pair<std::string,int>("m", 109));
   a.push_back(std::pair<std::string,int>("n", 110));
   a.push_back(std::pair<std::string,int>("o", 111));
   a.push_back(std::pair<std::string,int>("p", 112));
   a.push_back(std::pair<std::string,int>("q", 113));
   a.push_back(std::pair<std::string,int>("r", 114));
   a.push_back(std::pair<std::string,int>("s", 115));
   a.push_back(std::pair<std::string,int>("t", 116));
   a.push_back(std::pair<std::string,int>("u", 117));
   a.push_back(std::pair<std::string,int>("v", 118));
   a.push_back(std::pair<std::string,int>("w", 119));
   a.push_back(std::pair<std::string,int>("x", 120));
   a.push_back(std::pair<std::string,int>("y", 121));
   a.push_back(std::pair<std::string,int>("z", 122));
   a.push_back(std::pair<std::string,int>("braceleft", 123));
   a.push_back(std::pair<std::string,int>("bar", 124));
   a.push_back(std::pair<std::string,int>("braceright", 125));
   a.push_back(std::pair<std::string,int>("asciitilde", 126));
   a.push_back(std::pair<std::string,int>("nobreakspace", 160));
   a.push_back(std::pair<std::string,int>("exclamdown", 161));
   a.push_back(std::pair<std::string,int>("cent", 162));
   a.push_back(std::pair<std::string,int>("sterling", 163));
   a.push_back(std::pair<std::string,int>("currency", 164));

   return a; 
}
   
std::vector<std::pair<std::string, int> >
coot::util::key_sym_vec() {
   
   std::vector<std::pair<std::string, int> > a;

   std::vector<std::pair<std::string, int> > a_1 = coot::util::key_sym_vec_1();
   std::vector<std::pair<std::string, int> > a_2 = coot::util::key_sym_vec_2();
   std::vector<std::pair<std::string, int> > a_3 = coot::util::key_sym_vec_3();
   std::vector<std::pair<std::string, int> > a_4 = coot::util::key_sym_vec_4();

   for (unsigned int i=0; i<a_1.size(); i++)
      a.push_back(a_1[i]);
   for (unsigned int i=0; i<a_2.size(); i++)
      a.push_back(a_2[i]);
   for (unsigned int i=0; i<a_3.size(); i++)
      a.push_back(a_3[i]);
   for (unsigned int i=0; i<a_4.size(); i++)
      a.push_back(a_4[i]);

   // added by hand (from gdk-keysym.awk)
   
   a.push_back(std::pair<std::string,int>(":", 58));
   a.push_back(std::pair<std::string,int>(";", 59));
   a.push_back(std::pair<std::string,int>("<", 60));
   a.push_back(std::pair<std::string,int>("=", 61));
   a.push_back(std::pair<std::string,int>(">", 62));
   a.push_back(std::pair<std::string,int>("?", 63));
   a.push_back(std::pair<std::string,int>("@", 64));
   a.push_back(std::pair<std::string,int>("!", 33));
   a.push_back(std::pair<std::string,int>("$", 36));
   a.push_back(std::pair<std::string,int>("%", 37));
   a.push_back(std::pair<std::string,int>("&", 38));
   a.push_back(std::pair<std::string,int>("*", 42));
   a.push_back(std::pair<std::string,int>("+", 43));
   a.push_back(std::pair<std::string,int>(",", 44));
   a.push_back(std::pair<std::string,int>("-", 45));
   a.push_back(std::pair<std::string,int>(".", 46));
   a.push_back(std::pair<std::string,int>("/", 47));
   a.push_back(std::pair<std::string,int>("[", 91));
   a.push_back(std::pair<std::string,int>("]", 93));
   a.push_back(std::pair<std::string,int>("_", 95)); 
   a.push_back(std::pair<std::string,int>("|", 124)); 
   a.push_back(std::pair<std::string,int>("~", 126)); 

   return a;
}
