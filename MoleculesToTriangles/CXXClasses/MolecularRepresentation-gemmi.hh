/*
 * MoleculesToTriangles/CXXClasses/MolecularRepresentation-gemmi.hh
 *
 * gemmi-native twin of MolecularRepresentation.h (mmdb->gemmi migration).
 * Reuses the (mmdb-free) Representation base and DisplayPrimitive hierarchy;
 * consumes the coot::m2t MyMolecule / ColorScheme / DiscreteSegment twins and
 * gemmi-selection. Type is coot::m2t::MolecularRepresentation.
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 */

#ifndef MolecularRepresentation_gemmi_hh
#define MolecularRepresentation_gemmi_hh

#include <string>
#include <vector>
#include <memory>

#include "Representation.h"
#include "DisplayPrimitive.h"
#include "MyMolecule-gemmi.hh"
#include "ColorScheme-gemmi.hh"
#include "DiscreteSegment-gemmi.hh"
#include "MoleculesToTriangles/CXXClasses/gemmi-selection.hh"

class Renderer;

namespace coot {
   namespace m2t {

      class MolecularRepresentation : public Representation {
      private:
         std::shared_ptr<MyMolecule> myMolecule;
         std::shared_ptr<compound_selection_t> selection;
         std::shared_ptr<ColorScheme> colorScheme;
         std::string renderStyle;
         std::vector<DiscreteSegment *> segments;
      public:
         static std::string renderStyles[];

         MolecularRepresentation() : Representation(), myMolecule(nullptr), selection(nullptr), colorScheme(nullptr) {
            redrawNeeded = true;
            inRedraw = false;
            installDefaults();
         }
         MolecularRepresentation(std::shared_ptr<MyMolecule> _myMolecule,
                                 std::shared_ptr<ColorScheme> _colorScheme,
                                 std::shared_ptr<compound_selection_t> _compoundSelection,
                                 std::string _renderStyle)
            : Representation(), myMolecule(_myMolecule), selection(_compoundSelection),
              colorScheme(_colorScheme), renderStyle(_renderStyle) {
            redrawNeeded = true;
            inRedraw = false;
            installDefaults();
         }
         void installDefaults() {
            updateBoolParameter("doDraw", true);
            updateFloatParameter("ribbonStyleCoilThickness", 0.3);
            updateFloatParameter("ribbonStyleHelixWidth", 1.2);
            updateFloatParameter("ribbonStyleStrandWidth", 1.2);
            updateFloatParameter("ribbonStyleArrowWidth", 1.5);
            updateFloatParameter("ribbonStyleDNARNAWidth", 1.5);
            updateIntParameter  ("ribbonStyleAxialSampling", 6);
            updateIntParameter  ("cylindersStyleAngularSampling", 6);
            updateFloatParameter("cylindersStyleCylinderRadius", 0.2);
            updateIntParameter  ("dishStyleAngularSampling", 32);
            updateFloatParameter("cylindersStyleBallRadius", 0.2);
            updateFloatParameter("surfaceStyleProbeRadius", 1.4);
            updateFloatParameter("ballsStyleRadiusMultiplier", 1.);
            updateBoolParameter ("smoothBetas", true);
         }
         virtual ~MolecularRepresentation() {
            deletePrimitives();
            for (DiscreteSegment *s : segments) delete s;
         }

         int drawCalphas();
         int drawRibbon();
         int drawSpheres();
         int drawBondsAsSticks();
         int drawBondsAsNewSticks();
         int drawBondsAsCylinders();
         int drawSurfaceOfKind(int surfaceKind);   // TODO GEMMI Phase 5 (CXXSurface)
         int drawMolecularSurface();               // TODO GEMMI Phase 5
         int drawVdWSurface();                     // TODO GEMMI Phase 5
         int drawAccessibleSurface();              // TODO GEMMI Phase 5
         int drawHydrogenBonds();
         int drawDishyBases();
         int drawStickBases();

         void renderWithRenderer(Renderer *renderer);

         void deletePrimitives() { displayPrimitives.clear(); }
         std::shared_ptr<compound_selection_t> getCompoundSelection() { return selection; }
         void setCompoundSelection(std::shared_ptr<compound_selection_t> _selection) {
            selection = _selection; redrawNeeded = true;
         }
         void setColorScheme(std::shared_ptr<ColorScheme> _scheme) { colorScheme = _scheme; redrawNeeded = true; }
         std::shared_ptr<ColorScheme> getColorScheme() { return colorScheme; }

         // For ResidueNumber ramp rules, narrow each rule's start/end value to the
         // actual residue-number span of its selection (the gemmi equivalent of the
         // original's prepareForSelectionInMMDB), so the colour ramp spans the chain.
         void prepareRampRules();

         virtual void redraw() {
            if (selection != nullptr && colorScheme != nullptr && renderStyle != "" && myMolecule != nullptr) {
               deletePrimitives();
               prepareRampRules();
               if      (renderStyle == "Ribbon")           drawRibbon();
               else if (renderStyle == "Calpha")           drawCalphas();
               else if (renderStyle == "Sticks")           drawBondsAsNewSticks();
               else if (renderStyle == "Cylinders")        drawBondsAsCylinders();
               else if (renderStyle == "Spheres")          drawSpheres();
               else if (renderStyle == "MolecularSurface") drawMolecularSurface();
               else if (renderStyle == "VdWSurface")       drawVdWSurface();
               else if (renderStyle == "DishyBases")       drawDishyBases();
               else if (renderStyle == "StickBases")       drawStickBases();
               else if (renderStyle == "AccessibleSurface")drawAccessibleSurface();
               else if (renderStyle == "HydrogenBonds")    drawHydrogenBonds();
            }
            redrawNeeded = false;
         }
         void setRenderStyle(std::string _renderStyle) { renderStyle = _renderStyle; redrawNeeded = true; }
         std::string getRenderStyle() const { return renderStyle; }
         void setMolecule(std::shared_ptr<MyMolecule> _molecule) { myMolecule = _molecule; redrawNeeded = true; }
         std::shared_ptr<MyMolecule> getMolecule() { return myMolecule; }
         static int nStyles();
         static std::string styleName(int iStyle) { return renderStyles[iStyle]; }
         virtual bool getDoDraw() {
            bool result = boolParameters["doDraw"];
            if (result) {
               if (getMolecule() == nullptr) result = false;
               else result = getMolecule()->getDoDraw();
            }
            if (inRedraw) result = false;
            return result;
         }
         // TODO GEMMI Phase 6 (electrostatics): colorByOwnPotential / colorByPotential
      };
   }
}

#endif // MolecularRepresentation_gemmi_hh
