

camera.setFovy(5)
camera.setZClipFront(100)
camera.setFogFront(-100)
camera.setFogDepthRange(200)
camera.setFovy(3)
sceneSetup.getLight(0).setTranslation(CXXCoord_float(40.,40.,40,0.))

colorScheme = ColorScheme.colorByElementScheme()
colorScheme.addRule(SolidColorRule.colorRuleForSelectionStringAndName("/*/*/*.*/*:*","WHITE"))
colorScheme.addRule(SolidColorRule.colorRuleForSelectionStringAndName('("A") & [C]',"CORAL"))
colorScheme.addRule(SolidColorRule.colorRuleForSelectionStringAndName('("U") & [C]',"RED"))
colorScheme.addRule(SolidColorRule.colorRuleForSelectionStringAndName('("G") & [C]',"CYAN"))
colorScheme.addRule(SolidColorRule.colorRuleForSelectionStringAndName('("C") & [C]',"LAWNGREEN"))

import pypdb
pdbCodes = ["1ffk"]
for pdbCode in pdbCodes:
    pdbContents = pypdb.get_pdb_file(pdbCode)
    molecule = MyMolecule.createFromString(pdbContents)
    sceneSetup.setTranslation(molecule.getCentre()*-1.)
    molinst3=MolecularRepresentationInstance.create(molecule,ColorScheme.colorChainsScheme(),"ALL","DishyBases")
    sceneSetup.addRepresentationInstance(molinst3)
    molinst4=MolecularRepresentationInstance.create(molecule,ColorScheme.colorChainsScheme(),"ALL","Ribbon")
    molinst4.getRepresentation().updateFloatParameter("ribbonStyleCoilThickness",1.)
    molinst4.getRepresentation().updateFloatParameter("ribbonStyleDNARNAWidth",2.)
    sceneSetup.addRepresentationInstance(molinst4)
    self.wfWidget.repaint()
    
