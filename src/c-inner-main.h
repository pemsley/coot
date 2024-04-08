/*
 * src/c-inner-main.h
 *
 * Copyright 2007 by University of York
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

#ifndef BEGIN_C_DECLS

#ifdef __cplusplus
#define BEGIN_C_DECLS extern "C" {
#define END_C_DECLS }

#else
#define BEGIN_C_DECLS extern
#define END_C_DECLS     
#endif

#endif /* BEGIN_C_DECLS */

BEGIN_C_DECLS

// Bernhard turns these off for now.
#ifdef UNDERSTAND_GUILE_ON_WINDOWS
/* debugging with CYGWIN voodoo */
#if (defined(_WIN32) || defined(__CYGWIN__))
__declspec(dllimport) extern scm_t_option scm_debug_opts[];
__declspec(dllimport) extern scm_t_option scm_read_opts[];
__declspec(dllimport) extern scm_t_option scm_evaluator_trap_table[];
#endif
#endif

void c_inner_main(void *closure, int argc, char** argv); 
void c_wrapper_scm_boot_guile(int argc, char** argv);
char* does_file_exist (const char     *directory,
		       const char     *filename); 
void handle_command_line_data_argc_argv(int argc, char **argv);

END_C_DECLS
