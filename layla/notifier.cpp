/* layla/notifier.cpp
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

 #include "notifier.hpp"

class _CootLaylaNotifier {
    GObject base;
};

G_BEGIN_DECLS

G_DEFINE_TYPE(CootLaylaNotifier, coot_layla_notifier, G_TYPE_OBJECT)

static guint cif_file_generated_signal;

static void coot_layla_notifier_init(CootLaylaNotifier* self) {
    // This is the primary constructor
    
    // GObject doesn't run C++ constructors upon allocation
    // so we take care of this ourselves
    
}



static void coot_layla_notifier_dispose(GObject* _self) {
    CootLaylaNotifier* self = COOT_COOT_LAYLA_NOTIFIER(_self);
    // GObject doesn't run C++ destructors
    // so we take care of this ourselves
    
    
    G_OBJECT_CLASS(coot_layla_notifier_parent_class)->dispose(_self);
}

static void coot_layla_notifier_class_init(CootLaylaNotifierClass* klass) {
    // I think that this is a GObject class constructor that sets up the GObject class at runtime.

    cif_file_generated_signal = g_signal_new("cif-file-generated",
        G_TYPE_FROM_CLASS (klass),
        (GSignalFlags) (G_SIGNAL_RUN_LAST | G_SIGNAL_NO_RECURSE | G_SIGNAL_NO_HOOKS),
        0 /* class offset.Subclass cannot override the class handler (default handler). */,
        NULL /* accumulator */,
        NULL /* accumulator data */,
        NULL /* C marshaller. g_cclosure_marshal_generic() will be used */,
        G_TYPE_NONE /* return_type */,
        1     /* n_params */,
        G_TYPE_STRING
    );
    G_OBJECT_CLASS(klass)->dispose = coot_layla_notifier_dispose;

}

CootLaylaNotifier* 
coot_layla_notifier_new()
{
    return COOT_COOT_LAYLA_NOTIFIER(g_object_new (COOT_LAYLA_NOTIFIER_TYPE, NULL));
}

G_END_DECLS

void coot_layla_notifier_report_cif_file_generated(CootLaylaNotifier* self, const gchar* filepath) {
    g_signal_emit(self, cif_file_generated_signal, 0, filepath);
}
