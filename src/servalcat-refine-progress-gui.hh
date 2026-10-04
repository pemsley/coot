/*
 * src/servalcat-refine-progress-gui.hh
 *
 * Copyright 2026 by Medical Research Council
 * Author: Paul Emsley
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copies of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 * See http://www.gnu.org/licenses/
 */

#ifndef SERVALCAT_REFINE_PROGRESS_GUI_HH
#define SERVALCAT_REFINE_PROGRESS_GUI_HH

#include <string>

//! Show the servalcat refinement-progress overlay and clear any previous data.
//!
//! Called at the start of an asynchronous servalcat refinement. The stats file
//! (servalcat's "<prefix>_stats.json") is read on each subsequent call to
//! servalcat_refine_progress_update().
void servalcat_refine_progress_show(const std::string &stats_file_name);

//! Re-read the stats file into the shared progress data and redraw the graphs.
//!
//! Safe to call repeatedly (e.g. from the refinement's idle/timeout function) -
//! a stats file that is missing or caught mid-write leaves the current graphs in
//! place.
void servalcat_refine_progress_update();

//! Hide the overlay (the data is kept).
void servalcat_refine_progress_hide();

#endif // SERVALCAT_REFINE_PROGRESS_GUI_HH
