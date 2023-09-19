/* layla/notifier.hpp
 * 
 * Copyright 2023 by Global Phasing Ltd.
 * Author: Jakub Smulski
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */

#ifndef LAYLA_NOTIFIER_HPP
#define LAYLA_NOTIFIER_HPP
#include <glib-object.h>

// GObject declaration 
G_BEGIN_DECLS   

#define COOT_LAYLA_NOTIFIER_TYPE (coot_layla_notifier_get_type ())
G_DECLARE_FINAL_TYPE  (CootLaylaNotifier, coot_layla_notifier, COOT, COOT_LAYLA_NOTIFIER, GObject)

G_END_DECLS

extern "C" {
    CootLaylaNotifier* coot_layla_notifier_new();
}

void coot_layla_notifier_report_cif_file_generated(CootLaylaNotifier* self, const gchar* filepath);

#endif // LAYLA_NOTIFIER_HPP