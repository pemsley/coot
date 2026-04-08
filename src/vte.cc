/*
 * src/vte.cc
 *
 * Copyright 2025 by Medical Research Council
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
 *
 */

#ifdef HAVE_VTE

#include <gtk/gtk.h>
#include <vte/vte.h>
#include <glib-unix.h>
#include <sys/socket.h>
#include <unistd.h>
#include <fcntl.h>

#include <string>

#include <coot-utils/json.hpp>
#include "python-results-container.hh"
#include "utils/coot-utils.hh"
#include "widget-from-builder.hh"
#include "vte.hh"

using json = nlohmann::json;

// Static state for the VTE terminal
static int vte_parent_socket_fd = -1;
static GtkWidget *vte_terminal_widget = nullptr;
static GtkWidget *vte_paned_widget = nullptr;
static std::string vte_read_buffer;
static bool vte_helper_spawned = false;

static gboolean
on_vte_code_received(gint fd, GIOCondition condition, gpointer user_data) {

   if (condition & (G_IO_HUP | G_IO_ERR)) {
      g_warning("VTE REPL helper process disconnected");
      return G_SOURCE_REMOVE;
   }

   char buf[4096];
   ssize_t n = read(fd, buf, sizeof(buf) - 1);
   if (n <= 0) {
      g_warning("VTE socket read error or EOF");
      return G_SOURCE_REMOVE;
   }
   buf[n] = '\0';
   vte_read_buffer.append(buf, n);

   // Process complete newline-delimited messages
   size_t pos;
   while ((pos = vte_read_buffer.find('\n')) != std::string::npos) {
      std::string message = vte_read_buffer.substr(0, pos);
      vte_read_buffer.erase(0, pos + 1);

      if (message.empty()) continue;

      try {
         json request = json::parse(message);

         // Handle tab-completion requests
         if (request.contains("complete")) {
            std::string text = request.value("complete", "");

            // Use Python's rlcompleter in the embedded interpreter
            std::string completions_code =
               "import rlcompleter, __main__\n"
               "_c = rlcompleter.Completer(__main__.__dict__)\n"
               "_results = []\n"
               "_i = 0\n"
               "while True:\n"
               "    _m = _c.complete('" + text + "', _i)\n"
               "    if _m is None: break\n"
               "    _results.append(_m)\n"
               "    _i += 1\n";
            PyRun_SimpleString(completions_code.c_str());

            // Retrieve the results list
            PyObject *main_module = PyImport_AddModule("__main__");
            PyObject *global_dict = PyModule_GetDict(main_module);
            PyObject *results_obj = PyDict_GetItemString(global_dict, "_results");

            json response;
            response["completions"] = json::array();
            if (results_obj && PyList_Check(results_obj)) {
               Py_ssize_t size = PyList_Size(results_obj);
               for (Py_ssize_t i=0; i<size; i++) {
                  PyObject *item = PyList_GetItem(results_obj, i);
                  const char *s = PyUnicode_AsUTF8(item);
                  if (s) response["completions"].push_back(s);
               }

               // If there is exactly one match, fetch its docstring
               if (size == 1) {
                  PyObject *item = PyList_GetItem(results_obj, 0);
                  const char *match_str = PyUnicode_AsUTF8(item);
                  if (match_str) {
                     // Strip trailing '(' if rlcompleter added it
                     std::string match_name(match_str);
                     if (!match_name.empty() && match_name.back() == '(')
                        match_name.pop_back();
                     // Evaluate the name to get the object
                     PyObject *obj = PyRun_String(match_name.c_str(),
                                                  Py_eval_input, global_dict, global_dict);
                     if (obj) {
                        std::string doc_text;
                        if (PyCallable_Check(obj)) {
                           // Try to get the signature via inspect
                           PyObject *inspect = PyImport_ImportModule("inspect");
                           if (inspect) {
                              PyObject *sig_func = PyObject_GetAttrString(inspect, "signature");
                              if (sig_func) {
                                 PyObject *sig = PyObject_CallOneArg(sig_func, obj);
                                 if (sig) {
                                    PyObject *sig_str = PyObject_Str(sig);
                                    if (sig_str) {
                                       const char *s = PyUnicode_AsUTF8(sig_str);
                                       if (s) doc_text = match_name + std::string(s);
                                       Py_DECREF(sig_str);
                                    }
                                    Py_DECREF(sig);
                                 } else {
                                    PyErr_Clear();
                                 }
                                 Py_DECREF(sig_func);
                              }
                              Py_DECREF(inspect);
                           }
                        }
                        // Get the first line of the docstring
                        PyObject *doc_attr = PyObject_GetAttrString(obj, "__doc__");
                        if (doc_attr && doc_attr != Py_None) {
                           const char *ds = PyUnicode_AsUTF8(doc_attr);
                           if (ds) {
                              std::string first_line(ds);
                              auto nl = first_line.find('\n');
                              if (nl != std::string::npos)
                                 first_line = first_line.substr(0, nl);
                              if (!first_line.empty()) {
                                 if (!doc_text.empty())
                                    doc_text += "\n" + first_line;
                                 else
                                    doc_text = first_line;
                              }
                           }
                        }
                        if (doc_attr) Py_DECREF(doc_attr);
                        if (!doc_text.empty())
                           response["doc"] = doc_text;
                        Py_DECREF(obj);
                     } else {
                        PyErr_Clear();
                     }
                  }
               }
            }
            // Clean up temporary variables
            PyRun_SimpleString("del _c, _results, _i, _m");

            std::string response_str = response.dump() + "\n";
            write(fd, response_str.c_str(), response_str.size());

         } else {
            // Handle code execution requests
            std::string code = request.value("code", "");
            bool multiline = request.value("multiline", false);

            execute_python_results_container_t rc;
            if (multiline)
               rc = execute_python_multiline_code_with_result_internal(code);
            else
               rc = execute_python_code_with_result_internal(code);

            json response;
            response["stdout"] = rc.stdout;
            response["error"] = rc.error_message;

            if (rc.result && rc.result != Py_None) {
               PyObject *str_obj = PyObject_Repr(rc.result);
               if (str_obj) {
                  const char *str = PyUnicode_AsUTF8(str_obj);
                  if (str) response["result"] = str;
                  Py_DECREF(str_obj);
               }
            } else {
               response["result"] = "";
            }

            std::string response_str = response.dump() + "\n";
            ssize_t written = write(fd, response_str.c_str(), response_str.size());
            if (written < 0) {
               g_warning("VTE socket write error");
            }
         }

      } catch (const json::parse_error &e) {
         g_warning("VTE REPL: JSON parse error: %s", e.what());
         std::string err_response = "{\"stdout\":\"\",\"result\":\"\",\"error\":\"Internal JSON parse error\"}\n";
         write(fd, err_response.c_str(), err_response.size());
      }
   }

   return G_SOURCE_CONTINUE;
}

static void
on_vte_spawn_callback(VteTerminal *terminal, GPid pid, GError *error, gpointer user_data) {

   if (error) {
      g_warning("VTE REPL: Failed to spawn helper: %s", error->message);
      return;
   }
   if (pid != -1) {
      g_info("VTE REPL helper spawned with PID %d", pid);
   }
}

static void
on_vte_child_exited(VteTerminal *terminal, gint status, gpointer user_data) {

   g_info("VTE REPL helper exited with status %d", status);
}

// Spawn the Python helper process. Called on first reveal of the terminal,
// so that the banner and prompt are visible to the user.
static void spawn_vte_helper() {

   if (vte_helper_spawned) return;
   if (!vte_terminal_widget) return;
   vte_helper_spawned = true;

   VteTerminal *terminal = VTE_TERMINAL(vte_terminal_widget);

   // Create a socket pair for IPC
   int socket_fds[2];
   if (socketpair(AF_UNIX, SOCK_STREAM, 0, socket_fds) < 0) {
      g_warning("VTE REPL: Failed to create socket pair");
      return;
   }
   vte_parent_socket_fd = socket_fds[0]; // parent keeps this end
   int child_socket_fd = socket_fds[1];  // child gets this end

   // VTE requires fds passed to spawn_with_fds_async to have close-on-exec set.
   // VTE then explicitly maps them into the child process at the requested fd numbers.
   fcntl(child_socket_fd, F_SETFD, FD_CLOEXEC);

   // Watch the parent socket for incoming code execution requests
   g_unix_fd_add(vte_parent_socket_fd, G_IO_IN, on_vte_code_received, nullptr);

   // Find the helper script
   std::string helper_script = coot::package_data_dir() + "/python/coot_vte_repl.py";

   // The child will receive the socket fd.
   // We pass it as fd 3 in the child process.
   int child_fd_map = 3;
   std::string fd_str = std::to_string(child_fd_map);

   const char *argv[] = {"python3", helper_script.c_str(), fd_str.c_str(), nullptr};

   vte_terminal_spawn_with_fds_async(
      terminal,
      VTE_PTY_DEFAULT,
      nullptr,          // working directory (inherit)
      const_cast<char **>(argv),
      nullptr,          // envv (inherit)
      &child_socket_fd, 1,    // fds to pass to child
      &child_fd_map, 1,       // map child_socket_fd to fd 3
      G_SPAWN_SEARCH_PATH,
      nullptr, nullptr, nullptr, // child_setup, data, destroy
      -1,               // timeout (default)
      nullptr,          // cancellable
      on_vte_spawn_callback,
      nullptr           // user_data
   );

   // Close our copy of the child's socket end
   close(child_socket_fd);
}

void setup_python_vte_terminal() {

   // Create the VTE terminal widget
   vte_terminal_widget = vte_terminal_new();
   if (!vte_terminal_widget) {
      g_error("Failed to create VTE terminal widget");
      return;
   }

   // Configure the terminal
   VteTerminal *terminal = VTE_TERMINAL(vte_terminal_widget);
   PangoFontDescription *font_desc = pango_font_description_from_string("Monospace 10");
   vte_terminal_set_font(terminal, font_desc);
   pango_font_description_free(font_desc);
   vte_terminal_set_scrollback_lines(terminal, 1000);
   vte_terminal_set_scroll_on_output(terminal, TRUE);

   // Set dark-ish colors suitable for embedding
   GdkRGBA fg = {0.9, 0.9, 0.9, 1.0};
   GdkRGBA bg = {0.15, 0.15, 0.15, 1.0};
   vte_terminal_set_color_foreground(terminal, &fg);
   vte_terminal_set_color_background(terminal, &bg);

   g_signal_connect(terminal, "child-exited", G_CALLBACK(on_vte_child_exited), nullptr);

   // Add Cmd+V (macOS) / Ctrl+Shift+V paste support
   GtkEventController *key_controller = gtk_event_controller_key_new();
   g_signal_connect(key_controller, "key-pressed",
                    G_CALLBACK(+[](GtkEventControllerKey *controller,
                                   guint keyval, guint keycode,
                                   GdkModifierType state,
                                   gpointer user_data) -> gboolean {
                       if (keyval == GDK_KEY_v &&
                           (state & GDK_META_MASK || state & GDK_CONTROL_MASK)) {
                          VteTerminal *t = VTE_TERMINAL(user_data);
                          vte_terminal_paste_clipboard(t);
                          return TRUE;
                       }
                       return FALSE;
                    }),
                    terminal);
   gtk_widget_add_controller(vte_terminal_widget, key_controller);

   // Now reparent the graphics overlay into a GtkPaned with the VTE terminal
   GtkWidget *overlay = widget_from_builder("main_window_graphics_overlay");
   GtkWidget *parent = gtk_widget_get_parent(overlay);

   if (!parent) {
      g_error("VTE REPL: main_window_graphics_overlay has no parent");
      return;
   }

   // Create a vertical paned
   vte_paned_widget = gtk_paned_new(GTK_ORIENTATION_VERTICAL);
   gtk_paned_set_wide_handle(GTK_PANED(vte_paned_widget), TRUE);

   // Create a scrolled window for the terminal
   GtkWidget *scrolled = gtk_scrolled_window_new();
   gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scrolled),
                                  GTK_POLICY_AUTOMATIC, GTK_POLICY_AUTOMATIC);
   gtk_scrolled_window_set_child(GTK_SCROLLED_WINDOW(scrolled), vte_terminal_widget);
   gtk_widget_set_size_request(scrolled, -1, 150); // minimum height

   // Reparent: remove overlay from its parent box, insert paned at the same position,
   // then add overlay and terminal as children of the paned.
   // The parent is main_window_hbox (a GtkBox). The overlay is its first child —
   // we must preserve this ordering so the vertical toolbar stays on the right.
   g_object_ref(overlay); // prevent destruction during reparent
   gtk_box_remove(GTK_BOX(parent), overlay);

   gtk_paned_set_start_child(GTK_PANED(vte_paned_widget), overlay);
   gtk_paned_set_resize_start_child(GTK_PANED(vte_paned_widget), TRUE);
   gtk_paned_set_shrink_start_child(GTK_PANED(vte_paned_widget), FALSE);

   gtk_paned_set_end_child(GTK_PANED(vte_paned_widget), scrolled);
   gtk_paned_set_resize_end_child(GTK_PANED(vte_paned_widget), FALSE);
   gtk_paned_set_shrink_end_child(GTK_PANED(vte_paned_widget), FALSE);

   g_object_unref(overlay);

   // Prepend to preserve child order (overlay was the first child of the hbox)
   gtk_box_prepend(GTK_BOX(parent), vte_paned_widget);

   // Start with VTE hidden
   gtk_widget_set_visible(scrolled, FALSE);

   // Store the scrolled window so we can toggle visibility later
   g_object_set_data(G_OBJECT(vte_paned_widget), "vte_scrolled_window", scrolled);

   // Hide the old GtkEntry revealer
   GtkWidget *old_revealer = widget_from_builder("python_scripting_revealer");
   if (old_revealer)
      gtk_widget_set_visible(old_revealer, FALSE);
}

void toggle_vte_terminal_visibility() {

   if (!vte_paned_widget) return;
   GtkWidget *scrolled = GTK_WIDGET(g_object_get_data(G_OBJECT(vte_paned_widget),
                                                       "vte_scrolled_window"));
   if (!scrolled) return;

   gboolean visible = gtk_widget_get_visible(scrolled);
   gtk_widget_set_visible(scrolled, !visible);

   if (!visible) {
      spawn_vte_helper(); // spawn on first show
      if (vte_terminal_widget)
         gtk_widget_grab_focus(vte_terminal_widget);
   }
}

void show_vte_terminal() {

   if (!vte_paned_widget) return;
   GtkWidget *scrolled = GTK_WIDGET(g_object_get_data(G_OBJECT(vte_paned_widget),
                                                       "vte_scrolled_window"));
   if (!scrolled) return;
   gtk_widget_set_visible(scrolled, TRUE);
   spawn_vte_helper(); // spawn on first show
   if (vte_terminal_widget)
      gtk_widget_grab_focus(vte_terminal_widget);
}

#endif // HAVE_VTE
