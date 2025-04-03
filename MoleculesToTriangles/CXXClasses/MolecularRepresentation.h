/*
 * MoleculesToTriangles/CXXClasses/MolecularRepresentation.h
 *
 * Copyright 2009 by Martin Noble, University of Oxford
 * Author: Martin Noble
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
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
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

#include "compat/coot-sysdep.h"
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
        updateIntParameter  ("dishStyleAngularSampling",32);
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
    int drawStickBases();

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
       // std::cout << "++++++++++++++++++++++++++++++++++++++++++++ redraw with renderStyle " << renderStyle << std::endl;
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
                        else if (renderStyle == "StickBases"){
                                drawStickBases();
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
