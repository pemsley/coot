/*
 *  MolecularRepresentation.h
 *  Aesop
 *
 *  Created by Martin Noble on 17/02/2009.
 *  Copyright 2009 Dept. of Biochemistry, Oxford University. All rights reserved.
 *
 */
#ifndef MolecularRepresentation_h
#define MolecularRepresentation_h

#include <string>
#include <vector>
#include <memory>
#include <list>
#include <map>
#include "MyMolecule.h"
#include "ColorScheme.h"
#include "DisplayPrimitive.h"
#include "CompoundSelection.h"
#include "Representation.h"
#include "DiscreteSegment.h"

#include "mmdb2/mmdb_manager.h"

class Renderer;
class RepresentationInstance;

class MolecularRepresentation : public Representation {
private:
	//Modelled
    std::shared_ptr<MyMolecule> myMolecule;
    std::shared_ptr<CompoundSelection> selection;
    std::shared_ptr<ColorScheme> colorScheme;
    std::string renderStyle;
	//not modelled    
    int selHnd;
    std::vector<DiscreteSegment *> segments;
public:
    static std::string renderStyles[];
	//Also modelled 
    
    MolecularRepresentation() : Representation(), myMolecule(0), selection(0), colorScheme(0), selHnd(-999) {
        redrawNeeded = true;
        inRedraw = false;
        installDefaults();
    };
    MolecularRepresentation(std::shared_ptr<MyMolecule>_myMolecule, std::shared_ptr<ColorScheme> _colorScheme, std::shared_ptr<CompoundSelection> _compoundSelection, std::string _renderStyle) : Representation(), myMolecule(_myMolecule), selection(_compoundSelection), colorScheme(_colorScheme), renderStyle(_renderStyle), selHnd(-999) {
        redrawNeeded = true;
        inRedraw = false;
        installDefaults();
        //std::cout << "At inception " << getDoDraw();
    };
    void installDefaults(){
        updateBoolParameter("doDraw",true);
        //Here install render style defaults
        updateFloatParameter("ribbonStyleCoilThickness",0.3);
        updateFloatParameter("ribbonStyleHelixWidth",1.2);
        updateFloatParameter("ribbonStyleStrandWidth",1.2);
        updateFloatParameter("ribbonStyleArrowWidth",1.5);
        updateFloatParameter("ribbonStyleDNARNAWidth",1.5);
        updateIntParameter  ("ribbonStyleAxialSampling",6);
        updateIntParameter  ("cylindersStyleAngularSampling",6);
        updateFloatParameter("cylindersStyleCylinderRadius",0.2);
        updateFloatParameter("cylindersStyleBallRadius",0.2);        
        updateFloatParameter("surfaceStyleProbeRadius",1.4);
        updateFloatParameter("ballsStyleRadiusMultiplier",1.);
        updateBoolParameter ("smoothBetas",true);
    };
    virtual ~MolecularRepresentation(){
        deletePrimitives();
        std::vector<DiscreteSegment *>::iterator iter = segments.begin();
        std::vector<DiscreteSegment *>::iterator end = segments.end();
        for (; iter != end; ++iter){
            delete *iter;
        }
    };
    int drawCalphas();
    int drawRibbon();
	int drawSpheres();
	int drawBondsAsSticks();
	int drawBondsAsNewSticks();
	int drawBondsAsCylinders();
	int drawSurfaceOfKind(int surfaceKind);
	int drawMolecularSurface();
	int drawVdWSurface();
	int drawAccessibleSurface();
	int drawHydrogenBonds();
    int drawDishyBases();
    
    void renderWithRenderer(Renderer *renderer);
	
    void deletePrimitives(){
        /*
         for (unsigned long i=0; i< displayPrimitives.size(); i++){
            delete displayPrimitives[i];
            displayPrimitives[i] = 0;
        }
         */
        displayPrimitives.clear();
    };    
    std::shared_ptr<CompoundSelection> getCompoundSelection(){
        return selection;
    };
    void setCompoundSelection(std::shared_ptr<CompoundSelection> _selection){
        selection = _selection;
        redrawNeeded = true;
    };
    void setColorScheme(std::shared_ptr<ColorScheme> _scheme){
        colorScheme = _scheme;
        redrawNeeded = true;
    };
    std::shared_ptr<ColorScheme> getColorScheme(){
        return colorScheme;
    };
    virtual void redraw(){
		if (selection != 0 && colorScheme != 0 && renderStyle != "" && myMolecule != 0){
			deletePrimitives();
			if (renderStyle == "Ribbon"){
				drawRibbon();
			}
			else if (renderStyle == "Calpha"){
				drawCalphas();
			}
			else if (renderStyle == "Sticks"){
				drawBondsAsNewSticks();
			}
			else if (renderStyle == "Cylinders"){
				drawBondsAsCylinders();
			}
			else if (renderStyle == "Spheres"){
				drawSpheres();
			}
			else if (renderStyle == "MolecularSurface"){
				drawMolecularSurface();
			}
			else if (renderStyle == "VdWSurface"){
				drawVdWSurface();
			}
            else if (renderStyle == "DishyBases"){
                drawDishyBases();
            }
			else if (renderStyle == "AccessibleSurface"){
				drawAccessibleSurface();
			}
			else if (renderStyle == "HydrogenBonds"){
				drawHydrogenBonds();
			}
		}
        redrawNeeded = false;
    };
    void setRenderStyle(std::string _renderStyle){
        renderStyle = _renderStyle;
        redrawNeeded = true;
    };
    std::string getRenderStyle() const{
        return renderStyle;
    };
    void setMolecule(std::shared_ptr<MyMolecule> _molecule){
        myMolecule = _molecule;
        redrawNeeded = true;
    };
	std::shared_ptr<MyMolecule> getMolecule(){
		return myMolecule;
	};
    static int nStyles();
    static std::string styleName(int iStyle) {
        return renderStyles[iStyle];
    };
    virtual bool getDoDraw() {
        bool result = boolParameters["doDraw"];
        if (result) {
            if (getMolecule() == 0) result = false;
            else result = getMolecule()->getDoDraw();
        }
        if (inRedraw) result = false;
        return result;
    };
    void colorByOwnPotential();
    void colorByPotential(std::string compoundSelectionString, std::shared_ptr<MyMolecule> theMolecule);
};

#endif
