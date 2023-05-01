#  Copyright 2012 by the University of Oxford
#  Copyright 2014 by Medical Research Council
#  Author: Paul Emsley
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or (at
#  your option) any later version.
#
#  This program is distributed in the hope that it will be useful, but
#  WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
#  02110-1301, USA

import gi
gi.require_version('Gtk', '4.0')
from gi.repository import Gtk

# import pygtk, gtk, pango
import xml.etree.ElementTree as et
import string
import math
import types
import gzip
import coot

class PDB_Entry:

    def __init__(self, entry_attrib, xml_file_name):
        self.xml_file_name = xml_file_name
        self.pdbid                 = False
        self.Fo_Fc_correlation     = False
        self.IoverSigma            = False
        self.numMillerIndices      = False
        self.percent_RSRZ_outliers = False
        self.absolute_percentile_percent_RSRZ_outliers = False
        self.relative_percentile_percent_RSRZ_outliers = False
        self.percent_RSRZ_outliers                     = False
        self.acentric_outliers     = False
        self.centric_outliers      = False
        self.EDS_R                 = False
        self.EDS_resolution        = False
        self.EDS_resolution_low    = False
        self.DataAnisotropy        = False
        self.DataCompleteness      = False
        self.Fo_Fc_correlation     = False
        self.WilsonBaniso          = False
        self.WilsonBestimate       = False
        self.relative_percentile_RNAsuiteness = False
        self.absolute_percentile_RNAsuiteness = False
        self.RNAsuiteness          = False
        self.TwinL                 = False
        self.TwinFraction          = False
        self.TwinL2                = False
        self.xtriage_input_columns = False
        self.TransNCS              = False
        self.DCC_R                 = False
        self.DCC_refinement_program = False
        self.DCC_Rfree              = False
        self.absolute_percentile_DCC_Rfree = False
        self.relative_percentile_DCC_Rfree = False
        self.CCP4version        = False
        self.RefmacVersion      = False
        self.bulk_solvent_k     = False
        self.bulk_solvent_b     = False
        self.absolute_percentile_percent_rama_outliers = False
        self.absolute_percentile_percent_rota_outliers = False
        self.relative_percentile_percent_rama_outliers = False
        self.relative_percentile_percent_rota_outliers = False
        self.percent_rama_outliers                     = False
        self.percent_rota_outliers                     = False
        self.RestypesNotcheckedForBondAngleGeometry    = False

        try:
            self.pdbid = entry_attrib['pdbid']
        except KeyError as e:
            pass

        try:
            self.absolute_percentile_percent_RSRZ_outliers = entry_attrib['absolute-percentile-percent-RSRZ-outliers']
            self.relative_percentile_percent_RSRZ_outliers = entry_attrib['relative-percentile-percent-RSRZ-outliers']
            self.percent_RSRZ_outliers                     = entry_attrib['percent-RSRZ-outliers']
        except KeyError as e:
            # it is OK if a Entry doesnt have these
            pass

        try:
            self.acentric_outliers                         = entry_attrib['acentric_outliers']
            self.centric_outliers                          = entry_attrib['centric_outliers']
        except KeyError as e:
            # it is OK if a Entry doesnt have these
            pass

        try:
            self.EDS_R              = entry_attrib['EDS_R']
            self.EDS_resolution     = entry_attrib['EDS_resolution']
            self.EDS_resolution_low = entry_attrib['EDS_resolution_low']
        except KeyError as e:
            # it is OK if a Entry doesnt have these
            pass

        try:
            self.DataAnisotropy                            = entry_attrib['DataAnisotropy']
            self.DataCompleteness                          = entry_attrib['DataCompleteness']
            self.Fo_Fc_correlation                         = entry_attrib['Fo_Fc_correlation']
        except KeyError as e:
            # it is OK if a Entry doesnt have these
            pass

        try:
            self.WilsonBaniso                              = entry_attrib['WilsonBaniso']
            self.WilsonBestimate                           = entry_attrib['WilsonBestimate']
        except KeyError as e:
            # it is OK if a Entry doesnt have these
            pass

        try:
            self.relative_percentile_RNAsuiteness          = entry_attrib['relative-percentile-RNAsuiteness']
            self.absolute_percentile_RNAsuiteness          = entry_attrib['absolute-percentile-RNAsuiteness']
            self.RNAsuiteness                              = entry_attrib['RNAsuiteness']
        except KeyError as e:
            # it is OK if a Entry doesnt have these
            pass

        try:
            self.TwinL              = entry_attrib['TwinL']
            self.TwinFraction       = entry_attrib['TwinFraction']
            self.TwinL2             = entry_attrib['TwinL2']
        except KeyError as e:
            # it is OK if a Entry doesnt have these
            pass

        try:
            self.xtriage_input_columns                     = entry_attrib['xtriage_input_columns']
        except KeyError as e:
            # it is OK if a Entry doesnt have one
            pass

        try:
            self.TransNCS                                  = entry_attrib['TransNCS']
        except KeyError as e:
            # it is OK if a Entry doesnt have one
            pass

        try:
            self.DCC_R                                     = entry_attrib['DCC_R']
            self.DCC_refinement_program                    = entry_attrib['DCC_refinement_program']
            self.DCC_Rfree                                 = entry_attrib['DCC_Rfree']
            self.absolute_percentile_DCC_Rfree             = entry_attrib['absolute-percentile-DCC_Rfree']
            self.relative_percentile_DCC_Rfree             = entry_attrib['relative-percentile-DCC_Rfree']
        except KeyError as e:
            # it is OK if a Entry doesnt have one
            pass

        try:
            self.CCP4version                               = entry_attrib['CCP4version']
        except KeyError as e:
            # it is OK if a Entry doesnt have one
            pass

        try:
            self.RefmacVersion      = entry_attrib['RefmacVersion']
        except KeyError as e:
            # it is OK if a Entry doesnt have one
            pass

        try:
            self.bulk_solvent_k     = entry_attrib['bulk_solvent_k']
            self.bulk_solvent_b     = entry_attrib['bulk_solvent_b']
        except KeyError as e:
            # it is OK if a Entry doesnt have one
            pass


class ModelledSubgroup:

    # <bond-outlier atom0="CD" atom1="CE" mean="1.508" obs="1.561" stdev="0.025" z="2.11"/>
    class bond_outlier:
        def __init__(self, bond_attrib):
            # print 'bond_attrib: ', bond_attrib
            self.atom0  = False
            self.atom1  = False
            self.mean   = False
            self.obs    = False
            self.stddev = False
            self.z      = False
            try:
                self.atom0  = bond_attrib['atom0']
                self.atom1  = bond_attrib['atom1']
                self.mean   = bond_attrib['mean']
                self.obs    = bond_attrib['obs']
                self.stdev  = bond_attrib['stdev']
                self.z      = bond_attrib['z']
            except KeyError as e:
                pass

    # <mog-bond-outlier Zscore="8.78" atoms="C2',N1'" mean="1.46" mindiff="0.10" numobs="26" obsval="1.32" stdev="0.02"/>
    class mog_bond_outlier:
        def __init__(self, mog_bond_attrib):
            # print 'mog_bond_attrib: ', mog_bond_attrib
            self.Zscore  = False
            self.atoms   = False
            self.mean    = False
            self.mindiff = False
            self.numobs  = False
            self.obsval  = False
            self.stdev   = False
            try:
                self.Zscore  = mog_bond_attrib['Zscore']
                self.atoms   = mog_bond_attrib['atoms']
                self.mean    = mog_bond_attrib['mean']
                self.mindiff = mog_bond_attrib['mindiff']
                self.numobs  = mog_bond_attrib['numobs']
                self.obsval  = mog_bond_attrib['obsval']
                self.stdev   = mog_bond_attrib['stdev']
            except KeyError as e:
                pass

    # <angle-outlier atom0="C4" atom1="C5" atom2="C7" mean="119.000" obs="121.304" stdev="0.600" z="3.84"/>
    class angle_outlier:
        def __init__(self, angle_attrib):
            # print 'angle_attrib: ', angle_attrib
            self.atom0  = False
            self.atom1  = False
            self.atom2  = False
            self.mean   = False
            self.obs    = False
            self.stddev = False
            self.z      = False
            try:
                self.atom0  = angle_attrib['atom0']
                self.atom1  = angle_attrib['atom1']
                self.atom2  = angle_attrib['atom2']
                self.mean   = angle_attrib['mean']
                self.obs    = angle_attrib['obs']
                self.stdev  = angle_attrib['stdev']
                self.z      = angle_attrib['z']
            except KeyError as e:
                pass

    # <mog-angle-outlier Zscore="4.09" atoms="O1',C1,C2" mean="107.14" mindiff="0.92" numobs="97"
    #                    obsval="113.79" stdev="1.63"/>
    class mog_angle_outlier:
        def __init__(self, mog_angle_attrib):
            # print 'mog_angle_attrib: ', mog_angle_attrib
            self.Zscore  = False
            self.atoms   = False
            self.mean    = False
            self.mindiff = False
            self.numobs  = False
            self.obsval  = False
            self.stdev   = False
            try:
                self.Zscore  = mog_angle_attrib['Zscore']
                self.atoms   = mog_angle_attrib['atoms']
                self.mean    = mog_angle_attrib['mean']
                self.mindiff = mog_angle_attrib['mindiff']
                self.numobs  = mog_angle_attrib['numobs']
                self.obsval  = mog_angle_attrib['obsval']
                self.stdev   = mog_angle_attrib['stdev']
            except KeyError as e:
                pass


    # <clash atom="O" cid="116" clashmag="0.52" dist="2.09"/>
    class clash:
        def __init__(self, clash_attrib):
            # print 'clash_attrib: ', clash_attrib
            self.atom     = False
            self.cid      = False
            self.clashmag = False
            self.dist     = False
            try:
                self.atom     = clash_attrib['atom']
                self.cid      = clash_attrib['cid']
                self.clashmag = clash_attrib['clashmag']
                self.dist     = clash_attrib['dist']
            except KeyError as e:
                pass


    def __init__(self, subgroup):

        # individual bits that may or may not be there:
        self.bond_outliers      = []
        self.angle_outliers     = []
        self.mog_bond_outliers  = []
        self.mog_angle_outliers = []
        self.clashes            = []
        self.rsrz     = False
        self.RNAscore = False
        self.suite    = False
        self.flippable_sidechain = None # i.e. unset ATM.
        self.rama = False
        self.num_H_reduce = False

        subgroup_attrib = subgroup.attrib
        for child in subgroup:
            if child.tag == 'bond-outlier':
                bo = self.bond_outlier(child.attrib)
                if bo:
                    self.bond_outliers.append(bo)
            if child.tag == 'mog-bond-outlier':
                bo = self.mog_bond_outlier(child.attrib)
                if bo:
                    self.mog_bond_outliers.append(bo)
            if child.tag == 'angle-outlier':
                ao = self.angle_outlier(child.attrib)
                if ao:
                    self.angle_outliers.append(ao)
            if child.tag == 'mog-angle-outlier':
                ao = self.mog_angle_outlier(child.attrib)
                if ao:
                    self.mog_angle_outliers.append(ao)
            if child.tag == 'clash':
                c = self.clash(child.attrib)
                if c:
                    self.clashes.append(c)

        try:
            self.altcode   = subgroup_attrib['altcode']
            self.resnum    = subgroup_attrib['resnum']
            self.chain     = subgroup_attrib['chain']
            self.icode     = subgroup_attrib['icode']
            self.resname   = subgroup_attrib['resname']
            self.model     = subgroup_attrib['model']
            self.rscc      = False
            self.NatomsEDS = False
            self.owab      = False
            self.avgoccu   = False
            self.rsr       = False
            self.rama      = False

            try:
                self.rsr   = subgroup_attrib['rsr']
            except KeyError as e:
                # it is OK if a ModelledSubgroup doesnt have a rsrz
                pass

            try:
                self.avgoccu   = subgroup_attrib['avgoccu']
            except KeyError as e:
                # it is OK if a ModelledSubgroup doesnt have a rsrz
                pass

            try:
                self.rama      = subgroup_attrib['rama']
            except KeyError as e:
                # it is OK if a ModelledSubgroup doesnt have a rama
                pass

            try:
                self.owab = subgroup_attrib['owab']
            except KeyError as e:
                # it is OK if a ModelledSubgroup doesnt have a owab
                pass

            try:
                self.NatomsEDS = subgroup_attrib['NatomsEDS']
            except KeyError as e:
                # it is OK if a ModelledSubgroup doesnt have a NatomsEDS
                pass

            try:
                self.rscc = subgroup_attrib['rscc']
            except KeyError as e:
                # it is OK if a ModelledSubgroup doesnt have a rscc
                pass

            try:
                self.rsrz = subgroup_attrib['rsrz']
            except KeyError as e:
                # it is OK if a ModelledSubgroup doesnt have a rsrz
                pass
            try:
                self.RNAscore = subgroup_attrib['RNAscore']
            except KeyError as e:
                # it is OK if a ModelledSubgroup doesnt have a RNAscore
                pass
            try:
                self.suite = subgroup_attrib['suite']
            except KeyError as e:
                # it is OK if a ModelledSubgroup doesnt have a suite
                pass

            return None

        except KeyError as e:
            print('ModelledSubgroup: KeyError with key', e)

    # Is this a residue for which we want to create a button?
    #
    # return a string (for the button, if yes) otherwise return False.
    #
    def is_problematic_p(self):

       is_bad = False
       problems = []
       for bo in self.bond_outliers:
           try:
               if abs(float(bo.z)) > 2:
                   is_bad = True
                   s = "Bond Outlier " + bo.atom0 + " " + bo.atom1 + " z = " + bo.z
                   problems.append(s)
           except TypeError as e:
               pass

       for ao in self.angle_outliers:
           try:
               if abs(float(ao.z)) > 3:
                   is_bad = True
                   s = "Angle Outlier " + ao.atom0 + " " + ao.atom1 + " " + ao.atom2 + " z = " + ao.z
                   problems.append(s)
           except TypeError as e:
               pass

       for bo in self.mog_bond_outliers:
           try:
               if abs(float(bo.Zscore)) > 2:
                   is_bad = True
                   s = "Mogul-based Bond Outlier " + bo.atoms + ", z = " + bo.Zscore
                   problems.append(s)
           except TypeError as e:
               pass

       for ao in self.mog_angle_outliers:
           try:
               if abs(float(ao.Zscore)) > 2:
                   is_bad = True
                   s = "Mogul-based Angle Outlier " + ao.atoms + ", z = " + ao.Zscore
                   problems.append(s)
           except TypeError as e:
               pass

       for clash in self.clashes:
           try:
               if (clash.clashmag):
                   if float(clash.clashmag) >= 0.4:
                       is_bad = True
                       s = "Clash atom " + clash.atom + " score: " + clash.clashmag
                       problems.append(s)
           except TypeError as e:
               pass

       if self.rsrz:
           if float(self.rsrz) > 0.40:
               is_bad = True
               s = "Bad RSRZ " + self.rsrz
               problems.append(s)

       if self.RNAscore:
           if float(self.RNAscore) < 0.5:
               is_bad = True
               s = "Bad RNAScore " + self.RNAscore
               problems.append(s)

       if self.rama:
           if self.rama == 'OUTLIER':
               is_bad = True
               s = "Ramachandran Outlier"
               problems.append(s)

       if is_bad:
           return problems
       else:
           return False

# return validation info, which wraps entry_validation_info
#
def parse_wwpdb_validation_xml(xml_string):

    try:
        tree = et.fromstring(xml_string)
        modelled_subgroups = []
        for child in tree:
            if child.tag == 'Entry':
                xml_file_name = 'Fake-XML-Filename'
                entry_info = PDB_Entry(child.attrib, xml_file_name)
            if child.tag == 'ModelledSubgroup':
                subgroup = ModelledSubgroup(child)
                if subgroup:
                    modelled_subgroups.append(subgroup)
        validation_info = [entry_info, modelled_subgroups]
        return validation_info
    except IOError as e:
        print("error:", e)
        return False


class validation_entry_to_canvas:

    def __init__(self, entry_validation_info_in):

        self.entry_validation_info = entry_validation_info_in;
        self.bar_length = 300
        self.bar_height = 10
        # self.x_bar_offset = 130
        self.x_bar_offset = 160
        self.y_bar_offset =  20
        self.y_pixels_bar_spacing = 30
        self.setup_colour_bar_buff()
        self.abs_bar_width = 6
        self.abs_bar_height = int(self.bar_height * 1.6)
        self.vbox = False

        if self.entry_validation_info != False:
            window = Gtk.Window()
            title = "PDB Validation Report for " # ...
            if self.entry_validation_info.pdbid:
                 title += self.entry_validation_info.pdbid
            else:
                 title += self.entry_validation_info.xml_file_name
            window.set_title(title)
            self.vbox = Gtk.VBox(False, 0)
            self.vbox.set_border_width(5)
            h_sep = Gtk.HSeparator()
            da = Gtk.DrawingArea()
            da.set_size_request(560,280)
            close_button = Gtk.Button("  Close  ")
            hbox = Gtk.HBox(False, 0)
            self.vbox.append(da)
            self.vbox.append(h_sep)
            hbox.append(close_button)
            self.vbox.pack_end(hbox, False, False, 6)
            window.add(self.vbox)
            window.show_all()
            # da.connect("expose-event", self.on_drawing_area_expose)
            da.connect("configure-event", self.on_drawing_area_expose)
            close_button.connect("clicked", lambda a : window.destroy())
            self.pangolayout = da.create_pango_layout("")

    def setup_colour_bar_buff(self):

        c = 3*self.bar_length*self.bar_height*['\0']
        for j in range(self.bar_length):

            f_j = float(j)/float(self.bar_length)
            # we need g to go 255:255:0
            r = int(255*(1-math.pow(f_j, 5)))
            # we need g to go 0:255:0
            g_1 = f_j                #  0 : 0.5 : 1
            g_2 = 2 * (g_1 - 0.5)    # -1 : 0.  : 1
            g_3 = g_2 * g_2          #  1 : 0.  : 1
            g_4 = 1 - g_3            #  0 : 1.  : 0
            g = int(240*g_4)
            b = int(255*math.pow(f_j, 0.2))
            for i in range(self.bar_height):
                idx = 3*(self.bar_length*i + j)
                c[idx  ] = chr(r)
                c[idx+1] = chr(g)
                c[idx+2] = chr(b)
        # self.colour_bar_buff = string.join(c, '') 20211021-PE 
        # is this what I want (now)?
        # print("debug:: c", c)
        self.colour_bar_buff = ""
        self.colour_bar_buff.join(c)
        print("debug:: self.colour_bar_buff", self.colour_bar_buff)

    def on_drawing_area_expose(self, da, event):

        style = da.get_style()
        gc = style.fg_gc[Gtk.STATE_NORMAL]

        n_sliders = self.draw_sliders(da, gc)
        self.draw_top_labels(da, gc)
        self.draw_bottom_labels(da, gc, n_sliders) # Worse, Better
        self.draw_key(da, gc, n_sliders) # percentile box descriptions
        self.draw_eds_resolution_text(da, gc, n_sliders)


    def draw_eds_resolution_text(self, da, gc, n_sliders):
        if (self.entry_validation_info.EDS_resolution):
            pangolayout = da.create_pango_layout("")
            font_desc = pango.FontDescription('Sans 8')
            pangolayout.set_font_description(font_desc)
            s = 'EDS Resolution: ' + self.entry_validation_info.EDS_resolution
            pangolayout.set_text(s)
            y_level  = self.y_bar_offset + self.y_pixels_bar_spacing * n_sliders + 90
            da.window.draw_layout(gc, 10, y_level, pangolayout)


    def draw_rgb_image(self, da, gc, x, y):

        da.window.draw_rgb_image(gc, x, y, self.bar_length, self.bar_height,
                                 Gtk.gdk.RGB_DITHER_NONE, self.colour_bar_buff, self.bar_length*3)

    def bar_for_abs(self, abs_percent, y_min, da, gc):

        x = int(self.x_bar_offset + 0.01 * abs_percent * self.bar_length - 0.5 * self.abs_bar_width)
        y = int(y_min - self.bar_height * 0.25)
        da.window.draw_rectangle(gc, True, x, y, self.abs_bar_width, self.abs_bar_height)

    def arcs_for_rel(self, rel_percent, y_min, da, gc):
        x = int(self.x_bar_offset + 0.01 * rel_percent * self.bar_length - 0.5 * self.abs_bar_width)
        y = int(y_min - self.bar_height * 0.25)
        da.window.draw_rectangle(gc, False, x, y, self.abs_bar_width, self.abs_bar_height)

    # Worse, Better
    def draw_bottom_labels(self, da, gc, n_sliders):
        x_worse  = self.x_bar_offset
        y_worse  = self.y_bar_offset + self.y_pixels_bar_spacing * n_sliders + 15
        x_better = self.x_bar_offset + self.bar_length - 32
        y_better = y_worse

        pl_wb = da.create_pango_layout("")
        pl_wb.set_text('Worse')
        font_desc = pango.FontDescription('Sans 8')
        font_desc.set_style(pango.STYLE_ITALIC)
        pl_wb.set_font_description(font_desc)
        da.window.draw_layout(gc, x_worse,  y_worse,  pl_wb)

        pl_wb.set_text('Better')
        da.window.draw_layout(gc, x_better, y_better, pl_wb)

    # percentile box descriptions
    def draw_key(self, da, gc, n_sliders):
        x_key_box_abs =  self.x_bar_offset
        y_key_box_abs = self.y_bar_offset + self.y_pixels_bar_spacing * (n_sliders + 1) + 10

        x_key_box_rel = x_key_box_abs
        y_key_box_rel = y_key_box_abs + 20

        x_key_1 =  x_key_box_abs + 10
        y_key_1  = y_key_box_abs

        x_key_2 =  x_key_1
        y_key_2  = y_key_1 + 20

        da.window.draw_rectangle(gc, True, x_key_box_abs, y_key_box_abs,
                                 self.abs_bar_width, self.abs_bar_height)
        da.window.draw_rectangle(gc, False, x_key_box_rel, y_key_box_rel,
                                 self.abs_bar_width, self.abs_bar_height)

        pl = da.create_pango_layout("")
        font_desc = pango.FontDescription('Sans 9')
        pl.set_font_description(font_desc)

        pl.set_text('Percentile relative to all x-ray structures')
        da.window.draw_layout(gc, x_key_1, y_key_1, pl)
        pl.set_text('Percentile relative to x-ray structures of similar resolution')
        da.window.draw_layout(gc, x_key_2, y_key_2, pl)


    # return True if the bar for the absolute percentile was drawn,
    # otherwise return False
    #
    def draw_slider(self, name, x_for_rj, abs_str, rel_str, value_str, slider_no, da, gc):

        y = self.y_bar_offset + self.y_pixels_bar_spacing*(slider_no+1)

        # colour bar
        self.draw_rgb_image(da, gc, self.x_bar_offset, y)

        local_gc = gc;
        local_gc.set_foreground(local_gc.get_colormap().alloc_color("#888888"))

        da.window.draw_rectangle(local_gc, False, self.x_bar_offset, y, self.bar_length, self.bar_height)
        local_gc.set_foreground(local_gc.get_colormap().alloc_color("#000000"))

        # Metric text
        pangolayout = da.create_pango_layout("")
        pangolayout.set_justify(1)
        pangolayout.set_alignment(pango.ALIGN_RIGHT)
        pangolayout.set_text(name)
        pangolayout.context_changed()
        da.window.draw_layout(gc, 4+x_for_rj, y-6, pangolayout)

        # Values text
        if isinstance(value_str, types.StringType):
            x_for_value = self.x_bar_offset + self.bar_length + 12
            pangolayout.set_text(value_str)
            # print "value", x_for_value, y
            da.window.draw_layout(gc, x_for_value, y-6, pangolayout)

        # Bars for percentile scores
        try:
            if rel_str != False:
                rel = float(rel_str)
                self.arcs_for_rel(rel, y, da, gc)
            if abs_str != False:
                abs = float(abs_str)
                self.bar_for_abs( abs, y, da, gc)
                return True

            rel = float(rel_str)
        except TypeError as e:
            print(e)

        # hopefully we don't get here.
        return False


    # pass abs values
    def draw_connecting_lines(self, pc_ranks, da, gc):

        if len(pc_ranks) > 1:
            for slider_no in range(len(pc_ranks)-1):
                y_min_1 = 20 + 30*(slider_no+1)
                y_min_2 = 20 + 30*(slider_no+2)
                abs_1 = float(pc_ranks[slider_no])
                abs_2 = float(pc_ranks[slider_no+1])

                x_1 = int(self.x_bar_offset + 0.01 * abs_1 * self.bar_length)
                y_1 = int(y_min_1 + self.bar_height + 1)

                x_2 = int(self.x_bar_offset + 0.01 * abs_2 * self.bar_length)
                y_2 = int(y_min_2 - self.bar_height * 0.25)

                da.window.draw_line(gc, x_1, y_1, x_2, y_2)


    def draw_sliders(self, da, gc):

        save_abs = []
        icount_slider = 0

        abs = self.entry_validation_info.absolute_percentile_DCC_Rfree
        rel = self.entry_validation_info.relative_percentile_DCC_Rfree
        value = self.entry_validation_info.DCC_Rfree
        if abs != False and rel != False:
            state = self.draw_slider("Rfree", 114, abs, rel, value, icount_slider, da, gc)
            if state:
                save_abs.append(abs)
            icount_slider += 1

        abs   = self.entry_validation_info.absolute_percentile_clashscore
        rel   = self.entry_validation_info.relative_percentile_clashscore
        value = self.entry_validation_info.clashscore
        if abs != False and rel != False:
            state = self.draw_slider("Clashscore", 78, abs, rel, value, icount_slider, da, gc)
            if state:
                save_abs.append(abs)
            icount_slider += 1

        abs   = self.entry_validation_info.absolute_percentile_percent_rama_outliers
        rel   = self.entry_validation_info.relative_percentile_percent_rama_outliers
        value = self.entry_validation_info.percent_rama_outliers
        if abs != False and rel != False:
            state = self.draw_slider("Ramachandran Outliers", 0, abs, rel, value, icount_slider, da, gc)
            if state:
                save_abs.append(abs)
            icount_slider += 1

        abs = self.entry_validation_info.absolute_percentile_percent_rota_outliers
        rel = self.entry_validation_info.relative_percentile_percent_rota_outliers
        value = self.entry_validation_info.percent_rota_outliers
        if abs != False and rel != False:
            state = self.draw_slider("Sidechain Outliers", 32, abs, rel, value, icount_slider, da, gc)
            if state:
                save_abs.append(abs)
            icount_slider += 1

        abs = self.entry_validation_info.relative_percentile_percent_RSRZ_outliers
        rel = self.entry_validation_info.absolute_percentile_percent_RSRZ_outliers
        value = self.entry_validation_info.percent_RSRZ_outliers
        if abs != False and rel != False:
            state = self.draw_slider("RSRZ Outliers", 56, abs, rel, value, icount_slider, da, gc)
            if state:
                save_abs.append(abs)
            icount_slider += 1

        abs = self.entry_validation_info.relative_percentile_RNAsuiteness
        rel = self.entry_validation_info.absolute_percentile_RNAsuiteness
        value = self.entry_validation_info.RNAsuiteness
        if abs != False and rel != False:
            state = self.draw_slider("RNASuiteness", 57, abs, rel, value, icount_slider, da, gc)
            if state:
                save_abs.append(abs)
            icount_slider += 1

        self.draw_connecting_lines(save_abs, da, gc)
        return icount_slider

    def draw_top_labels(self, da, gc):
        pangolayout = da.create_pango_layout("")

        pangolayout = da.create_pango_layout("")
        font_desc = pango.FontDescription('Sans 13')
        pangolayout.set_font_description(font_desc)

        pangolayout.set_text('Metric')
        y_level = 15
        da.window.draw_layout(gc, 100, y_level, pangolayout)
        pangolayout.set_text('Percentile Ranks')
        da.window.draw_layout(gc, 245, y_level, pangolayout)
        pangolayout.set_text('Value')
        da.window.draw_layout(gc, 470, y_level, pangolayout)

def add_residue_buttons(subgroups, vbox, imol):

    def go_to_residue(button, residue_spec):
        # print "go to imol ", imol, "residue", residue_spec
        set_go_to_atom_molecule(imol)
        set_go_to_atom_from_res_spec(residue_spec)

    if vbox:
        vbox_residue_buttons = Gtk.VBox(False, 0)
        scrolled_win = Gtk.ScrolledWindow()
        # scrolled_win.set_policy(Gtk.POLICY_AUTOMATIC, Gtk.POLICY_ALWAYS)
        scrolled_win.add_with_viewport(vbox_residue_buttons)
        scrolled_win.set_size_request(-1, 200)
        for group in subgroups:
            p = group.is_problematic_p()
            if (p):
                ri_string = 'Residue '
                ri_string += group.chain
                ri_string += ' '
                ri_string += group.resnum
                ri_string += ' '
                ri_string += group.resname
                ri_string += ':'
                is_first = True
                for p_i in p:
                    ri_string += '\n'
                    ri_string += '    '
                    ri_string += p_i
                residue_button = Gtk.Button(ri_string)
                if group.icode == ' ':
                    group.icode = ''
                try:
                    r_n = int(group.resnum)
                    residue_spec = [group.chain, r_n, group.icode ]
                    residue_button.connect("clicked", go_to_residue, residue_spec)
                    vbox_residue_buttons.append(residue_button)
                    residue_button.show()
                except ValueError:
                    print('problem parsing', group.chain, group.resnum, group.icode)
        vbox.append(scrolled_win)
        scrolled_win.show()
        vbox_residue_buttons.show()

def validation_to_gui(entry_validation_info, subgroups, imol):
    vi = validation_entry_to_canvas(entry_validation_info)
    add_residue_buttons(subgroups, vi.vbox, imol)

def sort_subgroups(subgroups):
    return subgroups

def pdb_validate(accession_code, imol):

    try:
        if len(accession_code) == 4:
            middle_letters = accession_code[1] + accession_code[2]
            gz_file_name = "pdb-validation-" + accession_code + ".xml.gz"
            url = 'http://ftp.ebi.ac.uk/pub/databases/pdb/validation_reports/'
            url += middle_letters
            url += '/'
            url += accession_code
            url += '/'
            url += accession_code
            url += '_validation.xml.gz'
            status = coot.coot_get_url(url, gz_file_name)
            # turn the gz_file_name into a string
            vi = False
            try:
                gz = gzip.open(gz_file_name)
                xml_string = gz.read()
                vi = parse_wwpdb_validation_xml(xml_string)
                gz.close()
            except IOError as e:
                print("Error", e)
            if vi:
                entry_validation_info = vi[0]
                subgroups = vi[1]
                ss = sort_subgroups(subgroups)
                validation_to_gui(entry_validation_info, ss, imol)
            else:
                print("something went wrong when getting", url)
                s = "Something went wrong"
                add_status_bar_text(s)

        else:
            print('WARNING:: invalid accession code', accession_code)
    except KeyError as e:
        print("Failure to get validation for ", accession_code)

# xml_list = ["3NPQ.xml", "1FV2.xml", "2PE5.xml", "1HAK.xml", "2FGG.xml", "1CBS.xml"]

# xml_list = ["1FV2.xml"]
# xml_list = ["2PE5.xml"]
# xml_list = ["1HAK.xml"]
# xml_list = ["1CBS.xml"]
# xml_list = ["1CBS-1.xml"]
# xml_list = ["2FGG.xml"]
# # xml_list = ["3NPQ.xml", "1FV2.xml", "2PE5.xml", "1HAK.xml", "2FGG.xml", "1CBS.xml"]
# # xml_list = ["3NPQ.xml"] # DNA only?
# # xml_list = ["3NPQ.xml", "1FV2.xml", "2PE5.xml", "1HAK.xml", "2FGG.xml", "1CBS.xml"]

# for xml_file in xml_list:
#     vi = parse_wwpdb_validataion_xml(xml_file)
#     if vi:
#         entry_validation_info = vi[0]
#         subgroups = vi[1]
#         ss = sort_subgroups(subgroups)
#         validation_to_gui(entry_validation_info, ss)

# # coot_real_exit(0)
# gtk.main()

# pdb_validate('1cbs')


