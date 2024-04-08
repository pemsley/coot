/* src/gtk-utils.hh
 * -*-c++-*-
 *
 * Copyright 2023 by Global Phasing Ltd
 *
 * Author: Jakub Smulski
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
 * You should have received a copy of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 */

#ifndef GTK_UTILS_HH
#define GTK_UTILS_HH
#include <gtk/gtk.h>
#include <memory>

/// Utility class for creating progress bar popups, transient for the main window
class ProgressBarPopUp {
	GtkWindow* window;
	GtkProgressBar* progress_bar;

	public:
    ProgressBarPopUp(const std::string &title, const std::string &description) noexcept;
    ProgressBarPopUp(const ProgressBarPopUp&) = delete;
    ProgressBarPopUp(ProgressBarPopUp&&) noexcept;

	void pulse() noexcept;
	void set_fraction(float frac) noexcept;
        void set_text(const std::string &text) noexcept;
    /// The popup gets automatically closed upon deallocation
	~ProgressBarPopUp();
};

/// A thread-safe wrapper around `ProgressBarPopUp`.
/// It can be freely cloned and used across threads
class ProgressNotifier {
	std::shared_ptr<ProgressBarPopUp> progress_bar_popup;

	public:

	ProgressNotifier(std::shared_ptr<ProgressBarPopUp> popup) noexcept;
	/// This allows for reporting progress in a thread-safe way
	void update_progress(float frac) noexcept;
	/// This allows for reporting progress in a thread-safe way
	void pulse() noexcept;
	/// This allows for changing the label in a thread-safe way
        void set_text(const std::string &text) noexcept;

	~ProgressNotifier();

};

#endif
