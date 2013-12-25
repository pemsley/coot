
import pygtk, gtk, pango
import xml.etree.ElementTree as et
import string
import math

class PDB_Entry:

    def setup(self):
	self.other = 'Done'
	
    def __init__(self, entry_attrib):
	self.setup()
	self.Fo_Fc_correlation     = False
	self.IoverSigma            = False
	self.numMillerIndices      = False
	self.percent_RSRZ_outliers = False
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
        try: 
	    self.PDB_revision_number                       = entry_attrib['PDB-revision-number']
	    self.absolute_percentile_clashscore            = entry_attrib['absolute-percentile-clashscore']
	    self.relative_percentile_clashscore            = entry_attrib['relative-percentile-clashscore']
	    self.RestypesNotcheckedForBondAngleGeometry    = entry_attrib['RestypesNotcheckedForBondAngleGeometry']
	    self.attemptedValidationSteps                  = entry_attrib['attemptedValidationSteps']
	    self.clashscore                                = entry_attrib['clashscore']
	    self.num_H_reduce                              = entry_attrib['num-H-reduce']
	except KeyError as e:
	    print 'PDB_Entry: KeyError with key', e
	    return None

	try:
	    self.absolute_percentile_percent_rama_outliers = entry_attrib['absolute-percentile-percent-rama-outliers']
	    self.absolute_percentile_percent_rota_outliers = entry_attrib['absolute-percentile-percent-rota-outliers']
	    self.relative_percentile_percent_rama_outliers = entry_attrib['relative-percentile-percent-rama-outliers']
	    self.relative_percentile_percent_rota_outliers = entry_attrib['relative-percentile-percent-rota-outliers']
	    self.percent_rama_outliers                     = entry_attrib['percent-rama-outliers']
	    self.percent_rota_outliers                     = entry_attrib['percent-rota-outliers']
	except KeyError as e:
	    # it is OK if a Entry doesnt have these
	    pass

	# need data for these of course
	try:
	    self.Fo_Fc_correlation = entry_attrib['Fo_Fc_correlation']
	    self.IoverSigma        = entry_attrib['IoverSigma']
	    self.numMillerIndices                          = entry_attrib['numMillerIndices']
	except KeyError as e:
	    # it is OK if a Entry doesnt have these
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
	    self.rscc      = subgroup_attrib['rscc']
	    self.resnum    = subgroup_attrib['resnum']
	    self.chain     = subgroup_attrib['chain']
	    self.NatomsEDS = subgroup_attrib['NatomsEDS']
	    self.owab      = subgroup_attrib['owab']
	    self.avgoccu   = subgroup_attrib['avgoccu']
	    self.icode     = subgroup_attrib['icode']
	    self.resname   = subgroup_attrib['resname']
	    self.model     = subgroup_attrib['model']
	    self.rsr       = subgroup_attrib['rsr']
	    
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
	    print 'ModelledSubgroup: KeyError with key', e

def parse_wwpdb_validataion_xml(xml_file_name):

    tree = et.parse(xml_file_name)
    root = tree.getroot()
    modelled_subgroups = []
    for child in root:
	if child.tag == 'Entry':
	    entry_info = PDB_Entry(child.attrib)
	if child.tag == 'ModelledSubgroup':
	    subgroup = ModelledSubgroup(child)
	    if subgroup:
		modelled_subgroups.append(subgroup)
    return entry_info


def test_parse_xml(xml_file_name):

    print xml_file_name
    tree = et.parse(xml_file_name)
    root = tree.getroot()
    print root
    for child in root:
	# print child.tag, child.attrib
	if child.tag == 'Entry':
	    print '==========================================='
	    entry_info = PDB_Entry(child.attrib)
    print 'entry_info.other: ', entry_info.other
    print 'entry_info.RefmacVersion: ', entry_info.RefmacVersion


def test_parse_xml_aside_1(xml_file_name):

    tree = et.parse(xml_file_name)
    root = tree.getroot()
    print root
    # print dir(root)
    print root.tag
    for child in root:
	print child.tag, child.attrib

def draw_rgb_image(da, gc, x, y):
    b = 80*3*80*['\0']
    for i in range(80):
	for j in range(80):
	    b[3*80*i+3*j]   = chr(255)
	    b[3*80*i+3*j+1] = chr(int(255*(j/40-math.trunc(j/40))))
	    b[3*80*i+3*j+2] = chr(j*3)
    buff = string.join(b, '')
    da.window.draw_rgb_image(gc, x, y, 200, 20,
			     gtk.gdk.RGB_DITHER_NONE, buff, 80*3)
    return
    
def on_drawing_area_expose(da, event):
    style = da.get_style()
    gc = style.fg_gc[gtk.STATE_NORMAL]
    black_color = gtk.gdk.Color(red=0, green=0, blue=0)
    draw_rgb_image(da, gc, 60, 20)
    da.window.draw_rectangle(gc, False, 60, 20, 200, 20)


def validation_entry_to_canvas(entry_validation_info):
    
    window = gtk.Window(gtk.WINDOW_TOPLEVEL)
    window.set_title("Validation Report")
    vbox = gtk.VBox(False, 0)
    vbox.set_border_width(5)
    h_sep = gtk.HSeparator()
    da = gtk.DrawingArea()
    da.set_size_request(300,200)
    close_button = gtk.Button("  Close  ")
    vbox.pack_start(da,           False, 6)
    vbox.pack_start(h_sep,        False, 6)
    vbox.pack_start(close_button, False, 6)
    window.add(vbox)
    window.show_all()
    da.connect("expose-event", on_drawing_area_expose)

    style = da.get_style()
    gc = style.fg_gc[gtk.STATE_NORMAL]
    pangolayout = da.create_pango_layout("")

    black_color = gtk.gdk.Color(red=0, green=0, blue=0)
    da.window.draw_rectangle(gc, True, 10, 20, 30, 50)
    da.window.draw_point(gc, 30, 40)
    da.window.draw_layout(gc, 102, 100, pangolayout)


# parse_wwpdb_validataion_xml("3NPQ.xml")
# test_parse_xml("test.xml")
# c = TestComplex(1,2)
# c = TestComplex(1)

vi = parse_wwpdb_validataion_xml("3NPQ.xml")
validation_entry_to_canvas(vi)

# coot_real_exit(0)
gtk.main()

