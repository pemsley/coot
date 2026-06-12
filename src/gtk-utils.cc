/* src/gtk-utils.cc
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

#include "src/gtk-utils.hh"
#include <string>
#include "gtk/gtk.h"
#include "src/graphics-info.h"

ProgressBarPopUp::ProgressBarPopUp(const std::string &title, const std::string &description) noexcept
   : frame(nullptr), owns_window(true) {
   this->window = (GtkWindow*) gtk_window_new();
   this->progress_bar = (GtkProgressBar*) gtk_progress_bar_new();

   gtk_window_set_title(this->window, title.c_str());
   gtk_window_set_deletable(this->window, FALSE);
   graphics_info_t::set_transient_for_main_window(GTK_WIDGET(this->window));
   auto* box = gtk_box_new(GTK_ORIENTATION_VERTICAL, 10);
   gtk_widget_set_margin_bottom(box, 30);
   gtk_widget_set_margin_top(box,    30);
   gtk_widget_set_margin_start(box,  30);
   gtk_widget_set_margin_end(box,    30);

   gtk_window_set_child(this->window, box);
   gtk_box_append(GTK_BOX(box), gtk_label_new(description.c_str()));
   gtk_box_append(GTK_BOX(box), GTK_WIDGET(this->progress_bar));
   gtk_window_present(this->window);
}

ProgressBarPopUp::ProgressBarPopUp(GtkWidget *frame, GtkProgressBar *progress_bar) noexcept
   : window(nullptr), progress_bar(progress_bar), frame(frame), owns_window(false) {
}

ProgressBarPopUp::ProgressBarPopUp(ProgressBarPopUp&& other) noexcept {
   this->window = other.window;
   this->progress_bar = other.progress_bar;
   this->frame = other.frame;
   this->owns_window = other.owns_window;
   other.window = nullptr;
   other.progress_bar = nullptr;
   other.frame = nullptr;
}

void ProgressBarPopUp::pulse() noexcept {
   if(this->progress_bar) {
      gtk_progress_bar_pulse(this->progress_bar);
   }
}

void ProgressBarPopUp::set_fraction(float frac) noexcept {

   // closing thie window using the window manager button
   // deletes the widgets, I think, so this pointer goes bad.
   if(this->progress_bar) {
      if (GTK_IS_PROGRESS_BAR(this->progress_bar))
         gtk_progress_bar_set_fraction(this->progress_bar, frac);
   }
}

void ProgressBarPopUp::set_text(const std::string &text) noexcept {
   if(this->progress_bar) {
      gtk_progress_bar_set_show_text(this->progress_bar, TRUE);
      gtk_progress_bar_set_text(this->progress_bar, text.c_str());
   }
}

ProgressBarPopUp::~ProgressBarPopUp() {
   if (owns_window) {
      if (this->window) {
         gtk_window_close(this->window);
      }
   } else {
      if (this->frame) {
         g_idle_add(+[](gpointer user_data) {
            GtkWidget *frame = GTK_WIDGET(user_data);
            gtk_widget_set_visible(frame, FALSE);
            return FALSE;
         }, this->frame);
      }
   }
}

ProgressNotifier::ProgressNotifier(std::shared_ptr<ProgressBarPopUp> popup) noexcept
:progress_bar_popup(std::move(popup)),
 cancel_flag(std::make_shared<std::atomic<bool>>(false)) { }

void ProgressNotifier::update_progress(float frac) noexcept {

   struct callback_data {
      std::shared_ptr<ProgressBarPopUp> popup;
      float frac;
   };
   callback_data* data = new callback_data{this->progress_bar_popup, frac};
   // This guarantes execution from the main thread
   g_idle_add(+[](gpointer user_data){
      callback_data* data = (callback_data*) user_data;
      data->popup->set_fraction(data->frac);
      delete data;
      return FALSE;
   }, data);
}

void ProgressNotifier::pulse() noexcept {

   struct callback_data {
      std::shared_ptr<ProgressBarPopUp> popup;
   };
   callback_data* data = new callback_data{this->progress_bar_popup};
   // This guarantes execution from the main thread
   g_idle_add(+[](gpointer user_data){
      callback_data* data = (callback_data*) user_data;
      data->popup->pulse();
      delete data;
      return FALSE;
   }, data);
}

void ProgressNotifier::set_text(const std::string &text) noexcept {

   struct callback_data {
      std::shared_ptr<ProgressBarPopUp> popup;
      std::string text;
   };
   callback_data* data = new callback_data{this->progress_bar_popup, std::string(text)};
   // This guarantes execution from the main thread
   g_idle_add(+[](gpointer user_data){
      callback_data* data = (callback_data*) user_data;
      data->popup->set_text(data->text.c_str());
      delete data;
      return FALSE;
   }, data);
}

void ProgressNotifier::cancel() noexcept {
   cancel_flag->store(true);
}

bool ProgressNotifier::is_cancelled() const noexcept {
   return cancel_flag->load();
}

ProgressNotifier::~ProgressNotifier() {

   struct callback_data {
      std::shared_ptr<ProgressBarPopUp> popup;
   };
   callback_data* data = new callback_data{std::move(this->progress_bar_popup)};
   // This guarantes that de-allocation happens on the main thread
   g_idle_add(+[](gpointer user_data) {
      callback_data* data = (callback_data*) user_data;
      delete data;
      return FALSE;
   }, data);
}
