//
//  MyMolecule_BoostWrapper.cpp
//  MoleculesToTriangles
//
//  Created by Martin Noble on 24/01/2017.
//  Copyright © 2017 MartinNobleSoftware. All rights reserved.
//

#include <stdio.h>
#include <boost/python.hpp>
using namespace boost::python;
#include "MyMolecule.h"
#include "CompoundSelection.h"
#include "ColorScheme.h"
#include "ColorRule.h"
#include "SolidColorRule.h"
#include "RepresentationInstance.h"
#include "MolecularRepresentationInstance.h"
#include "Light.h"
#include "SceneSetup.h"
#include "Camera.h"
//#include "CameraPort.h"
#include "oglPolyhedron.h"
//#include "LinesPrimitive.h"
#include "BondsPrimitive.h"
#include "Renderer.h"
#include "RendererGL.h"
#include "RendererGLSL.hpp"
#include "RotatedTranslatedScaledEntity.h"
#include "MoleculesToTriangles/CXXSurface/CXXCoord.h"
#include "mmdb2/mmdb_manager.h"


struct RepresentationWrap : Representation, wrapper<Representation>
{
    void updateFloatParameter(std::string parameterName,const float parameterValue){
        this->get_override("updateFloatParameter")(parameterName, parameterValue);
    };
    void updateIntParameter(std::string parameterName, int parameterValue){
        this->get_override("updateIntParameter")(parameterName, parameterValue);
    };
    void updateBoolParameter(std::string parameterName, bool parameterValue){
        this->get_override("updateBoolParameter")(parameterName, parameterValue);
    };
    bool getDoDraw(){
        return this->get_override("getDoDraw")();
    };
    bool setDoDraw(const bool yesOrNo){
        return this->get_override("getDoDraw")(yesOrNo);
    };
    void redraw(){
        this->get_override("redraw")();
    };
};

struct ColorRuleWrap : ColorRule, wrapper<ColorRule>
{
    virtual FCXXCoord colorForAtom (const mmdb::Atom *atom) {
        std::cout <<"Here we are\n";
        return this->get_override("colorForAtom")(atom);
    };
};


struct RendererWrap : Renderer, wrapper<Renderer>
{
    virtual void drawTestTriangle(const FCXXCoord &coord) const {
        this->get_override("drawTestTriangle")(coord);
    };
    virtual void vboRenderVCNFixed(VertexColorNormalPrimitive *prim) {
        this->get_override("vboRenderVCNFixed")(prim);
    };
    virtual void vboRenderVCN(VertexColorNormalPrimitive *prim) {
        this->get_override("vboRenderVCNFixed")(prim);
    };
    virtual void vboRenderVN(VertexColorNormalPrimitive *prim) {
        this->get_override("vboRenderVN")(prim);
    };
    virtual void vboRenderVC(VertexColorNormalPrimitive *prim) {
        this->get_override("vboRenderVC")(prim);
    };
    virtual void setupCamera(Camera *anArg) {
        this->get_override("setupCamera")(anArg);
    };
    virtual void setupScene(SceneSetup *anArg) {
        this->get_override("setupScene")(anArg);
    };
    virtual void init() {
        this->get_override("init")();
    };
    virtual void resize(int w, int h) {
        this->get_override("resize")(w, h);
    };
    virtual void setupLightAsIndexFromViewpointWithScale(Light *light, int asIndex, FCXXCoord  fromViewpoint, float scale) {
        this->get_override("setupLightAsIndexFromViewpointWithScale")(light, asIndex, fromViewpoint, scale);
    };
    virtual void setupRepresentationInstance(RepresentationInstance *instance) {
        this->get_override("setupRepresentationInstance")(instance);
    };
    virtual void restoreModelviewMatrix() {
        this->get_override("restoreModelviewMatrix")();
    };
    virtual void renderBondsPrimitive(BondsPrimitive *anArg) {
        this->get_override("renderBondsPrimitive")(anArg);
    };
    /*
     virtual void renderLinesPrimitive(LinesPrimitive *anArg) {
     this->get_override("renderLinesPrimitive")(anArg);
     };
     */
    virtual void renderPolyhedron(oglPolyhedron *anArg) {
        this->get_override("renderPolyhedron")(anArg);
    };
    virtual void renderVertexColorNormalPrimitive(VertexColorNormalPrimitive *anArg) {
        this->get_override("renderVertexColorNormalPrimitive")(anArg);
    };
    virtual void renderVertexColorPrimitive(VertexColorNormalPrimitive *anArg) {
        this->get_override("renderVertexColorPrimitive")(anArg);
    };
    virtual FCXXCoord  angstromsForPixels(float x, float y) {
        return this->get_override("angstromsForPixels")(x,y);
    };
    virtual void render(Camera *camera) {
        this->get_override("render")(camera);
    };
    virtual void liberateHandlesForDisplayPrimitive(DisplayPrimitive *prim) {
        this->get_override("liberateHandlesForDisplayPrimitive")(prim);
    };
    virtual void clearCameraCanvas(Camera *camera) {
        this->get_override("clearCameraCanvas")(camera);
    };
};

BOOST_PYTHON_MODULE(libBoostedMoleculesToTrianglesCXXClasses)
{
    
    class_<mmdb::Manager>("Manager")
    ;
    
    class_<MyMolecule, std::shared_ptr<MyMolecule> >("MyMolecule",init<std::string>())
    .def("setDoDraw", &MyMolecule::setDoDraw)
    .def("getCentre",&MyMolecule::getCentre)
    .def("centreOfSelectionString",&MyMolecule::centreOfSelectionString)
    .def("writePDB",&MyMolecule::writePDB)
    .def("getMmdb",&MyMolecule::getMmdb, return_value_policy<return_opaque_pointer>())
    .def(self_ns::str(self))                     // __str__
    .def("create",&MyMolecule::create)
    .staticmethod("create")
    .def("createFromString",&MyMolecule::createFromString)
    .staticmethod("createFromString")
    ;
    
    class_<CompoundSelection, std::shared_ptr<CompoundSelection> >("CompoundSelection",init<std::string,std::string>())
    .def("create", &CompoundSelection::create)
    .def("describe", &CompoundSelection::describe)
    .def("deleteInMMDB", &CompoundSelection::deleteInMMDB)
    .staticmethod("create")
    ;
    
    class_<ColorScheme, std::shared_ptr<ColorScheme> >("ColorScheme")
    .def("colorByElementScheme", &ColorScheme::colorByElementScheme)
    .staticmethod("colorByElementScheme")
    .def("colorBySecondaryScheme", &ColorScheme::colorBySecondaryScheme)
    .staticmethod("colorBySecondaryScheme")
    .def("colorBFactorScheme", &ColorScheme::colorBFactorScheme)
    .staticmethod("colorBFactorScheme")
    .def("colorRampChainsScheme", &ColorScheme::colorRampChainsScheme)
    .staticmethod("colorRampChainsScheme")
    .def("colorChainsScheme", &ColorScheme::colorChainsScheme)
    .staticmethod("colorChainsScheme")
    .def("colorSchemeForColorName", &ColorScheme::colorSchemeForColorName)
    .staticmethod("colorSchemeForColorName")
    .def("addRule", &ColorScheme::addRule)
    ;

    class_<SolidColorRule, std::shared_ptr<SolidColorRule> >("SolidColorRule")
    .def("colorRuleForSelectionStringAndName", &SolidColorRule::colorRuleForSelectionStringAndName)
    .staticmethod("colorRuleForSelectionStringAndName")
    .def("colorRuleForSelectionStringAndColor", &SolidColorRule::colorRuleForSelectionStringAndColor)
    .staticmethod("colorRuleForSelectionStringAndColor")
    ;
    
    class_<RotatedTranslatedScaledEntity, std::shared_ptr<RotatedTranslatedScaledEntity> >("RotatedTranslatedScaledEntity")
    .def("setTranslation", &RotatedTranslatedScaledEntity::setTranslation )
    .def("getTranslation", &RotatedTranslatedScaledEntity::getTranslation )
    .def("setRotation", &RotatedTranslatedScaledEntity::setRotation )
    .def("getRotation", &RotatedTranslatedScaledEntity::getRotation )
    .def("setScale", &RotatedTranslatedScaledEntity::setScale )
    .def("getScale", &RotatedTranslatedScaledEntity::getScale )
    .def("translateBy", &RotatedTranslatedScaledEntity::translateBy )
    .def("rotateBy", &RotatedTranslatedScaledEntity::rotateBy )
    ;
    
    class_<RepresentationWrap, std::shared_ptr<RepresentationWrap>, boost::noncopyable >("Representation")
    .def("updateFloatParameter", &RepresentationWrap::updateFloatParameter)
    .def("updateIntParameter", &RepresentationWrap::updateIntParameter)
    .def("updateBoolParameter", &RepresentationWrap::updateBoolParameter)
    .def("getDoDraw", pure_virtual(&RepresentationWrap::getDoDraw))
    ;
    
    class_<MolecularRepresentation, std::shared_ptr<MolecularRepresentation>, bases<RepresentationWrap> >("MolecularRepresentation")
    .def("setRenderStyle", &MolecularRepresentation::setRenderStyle)
    .def("setColorScheme", &MolecularRepresentation::setColorScheme)
    .def("setMolecule", &MolecularRepresentation::setMolecule)
    .def("setCompoundSelection", &MolecularRepresentation::setCompoundSelection)
    .def("colorByPotential", &MolecularRepresentation::colorByPotential)
    .def("colorByOwnPotential", &MolecularRepresentation::colorByOwnPotential)
    .def("updateFloatParameter", &MolecularRepresentation::updateFloatParameter)
    .def("updateIntParameter", &MolecularRepresentation::updateIntParameter)
    .def("updateBoolParameter", &MolecularRepresentation::updateBoolParameter)
    ;
    
    class_<RepresentationInstance, std::shared_ptr<RepresentationInstance>, bases< RotatedTranslatedScaledEntity> >("RepresentationInstance")
    .def("getRepresentation", &RepresentationInstance::getRepresentation)
    .def("setDoDraw", &RepresentationInstance::setDoDraw)
    .def("setPostTranslation", &RepresentationInstance::setPostTranslation)
    .def("setPostRotation", &RepresentationInstance::setPostRotation)
    ;
    
    class_<MolecularRepresentationInstance, std::shared_ptr<MolecularRepresentationInstance>, bases< RepresentationInstance> >("MolecularRepresentationInstance",init<std::shared_ptr<MyMolecule>, std::shared_ptr<ColorScheme>, std::shared_ptr<CompoundSelection>, std::string>())
    .def("getRepresentation", &MolecularRepresentationInstance::getRepresentation)
    .def("setDoDraw", &MolecularRepresentationInstance::setDoDraw)
    .def("create", &MolecularRepresentationInstance::create)
    .staticmethod("create")
    ;
    
    class_<Light, std::shared_ptr<Light>, bases<RotatedTranslatedScaledEntity> >("Light")
    .def("defaultLight", &Light::defaultLight)
    .def("setDrawLight", &Light::setDrawLight)
    .def("setIntensity", &Light::setIntensity)
    .def("getIntensity", &Light::getIntensity)
    .def("getAmbient", &Light::getAmbient)
    .def("getSpecular", &Light::getSpecular)
    .def("getDiffuse", &Light::getDiffuse)
    .def("setAmbient", &Light::setAmbient)
    .def("setSpecular", &Light::setSpecular)
    .def("setDiffuse", &Light::setDiffuse)
    .def("setLightType", &Light::setLightType)
    .def("getLightType", &Light::getLightType)
    .staticmethod("defaultLight")
    ;
    
    class_<SceneSetup, std::shared_ptr<SceneSetup>, bases<RotatedTranslatedScaledEntity> >("SceneSetup")
    .def("defaultSceneSetup", &SceneSetup::defaultSceneSetup)
    .def("addCamera", &SceneSetup::addCamera)
    .def("setBackgroundColor", &SceneSetup::setBackgroundColor)
    .def("addLight", &SceneSetup::addLight)
    .def("getLight", &SceneSetup::getLight)
    .def("getRepresentationInstance", &SceneSetup::getRepresentationInstance)
    .def("addRepresentationInstance", &SceneSetup::addRepresentationInstance)
    .def("removeRepresentationInstance", &SceneSetup::removeRepresentationInstance)
    .staticmethod("defaultSceneSetup")
    ;
    
    class_<Camera, std::shared_ptr<Camera>, bases<RotatedTranslatedScaledEntity> >("Camera")
    .def("renderWithRenderer", &Camera::renderWithRenderer)
    .def("setSceneSetup", &Camera::setSceneSetup)
    .def("getSceneSetup", &Camera::getSceneSetup)
    .def("defaultCamera", &Camera::defaultCamera)
    .def("getFovy", &Camera::getFovy)
    .def("setFovy", &Camera::setFovy)
    .def("getFogFront", &Camera::getFogFront)
    .def("setFogFront", &Camera::setFogFront)
    .def("getFogDepthRange", &Camera::getFogDepthRange)
    .def("setFogDepthRange", &Camera::setFogDepthRange)
    .def("getZClipFront", &Camera::getZClipFront)
    .def("setZClipFront", &Camera::setZClipFront)
    .staticmethod("defaultCamera")
    ;
    
    /*
     class_<CameraPort, std::shared_ptr<CameraPort> >("CameraPort")
     .def("setCamera", &CameraPort::setCamera)
     .def("setRendererGL", &CameraPort::setRendererGL)
     .def("runLoop", &CameraPort::runLoop)
     .def("runThreaded", &CameraPort::runThreaded)
     .def("create", &CameraPort::create)
     .staticmethod("create")
     ;
     */
    
    class_<RendererWrap, std::shared_ptr<RendererWrap>, boost::noncopyable >("Renderer")
    ;
    
    class_<RendererGL, std::shared_ptr<RendererGL>, bases<Renderer> >("RendererGL")
    .def("create", &RendererGL::create)
    .def("init", &RendererGL::init)
    .def("resize", &RendererGL::resize)
    .staticmethod("create")
    ;
    
    class_<FCXXCoord  >("CXXCoord_float")
    .def(init<CXXCoord_ftype,CXXCoord_ftype,CXXCoord_ftype>())
    .def(init<CXXCoord_ftype,CXXCoord_ftype,CXXCoord_ftype,CXXCoord_ftype>())
    .def("dump", &FCXXCoord ::dump)
    .def(self *= CXXCoord_ftype())
    .def(self * CXXCoord_ftype())
    .def(self / CXXCoord_ftype())
    .def(self += FCXXCoord ())
    .def(self + FCXXCoord ())
    .def(self - FCXXCoord ())
    ;
    
    class_<RendererGLSL, std::shared_ptr<RendererGLSL>, bases<RendererGL> >("RendererGLSL",init<>())
    .def("test", &RendererGLSL::test)
    .def("init", &RendererGLSL::init)
    .def("setProgram", &RendererGLSL::setProgram)
    ;
    
    implicitly_convertible<std::shared_ptr<RendererGL>,std::shared_ptr<Renderer> >();
    implicitly_convertible<std::shared_ptr<RendererGLSL>,std::shared_ptr<RendererGL> >();
    implicitly_convertible<std::shared_ptr<MolecularRepresentationInstance>,std::shared_ptr<RepresentationInstance> >();
    implicitly_convertible<std::shared_ptr<MolecularRepresentation>,std::shared_ptr<Representation> >();
    implicitly_convertible<std::shared_ptr<SolidColorRule>,std::shared_ptr<ColorRule> >();
}

