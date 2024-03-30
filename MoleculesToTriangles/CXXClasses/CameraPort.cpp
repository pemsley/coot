//
//  CameraPort.cpp
//  MoleculesToTriangles
//
//  Created by Martin Noble on 27/01/2017.
//  Copyright © 2017 MartinNobleSoftware. All rights reserved.
//

#include "CameraPort.h"
#include "SceneSetup.h"
#include <thread>
#include <set>
#include "ColorScheme.h"
#include "MolecularRepresentation.h"
#include "Light.h"

void CameraPort::initialiseMoleculesToTriangles(){
    
    //Initialize your renderer...this does things like compiling and linking shaders and should be done only once
    renderer = RendererGLSL::create();
    
    //Generate a sceneSetup...this stores lights and representaiton instances
    sceneSetup = SceneSetup::defaultSceneSetup();
    //Because I am using my "Camera" class to define things like filed of view, perspective etc, I have to do this, but
    //applicaitons which take control of their own model view matrix don't
    camera->setSceneSetup(sceneSetup);
    
    //Add a simple light..set some parameters and move it
    auto simpleLight = Light::defaultLight();
    simpleLight->setIntensity(0.75);
    simpleLight->setDrawLight(false);
    sceneSetup->addLight(simpleLight);
    //Can retrieve reference to the light if so preferred
    sceneSetup->getLight(0)->setTranslation(FCXXCoord(400.,400.,0,0.));
    
    //Add another simple light
    auto simpleLight2 = Light::defaultLight();
    sceneSetup->addLight(simpleLight2);
    simpleLight2->setIntensity(0.5);
    simpleLight2->setDrawLight(false);
    simpleLight2->setTranslation(FCXXCoord(0.,0.,400,0.));
    
}

std::shared_ptr<MolecularRepresentationInstance> CameraPort::addRepresentationInstance(mmdb::Manager *ofMmdb, std::shared_ptr<ColorScheme> colorScheme, std::string selectionString, std::string renderStyle)
{

    MyMolecule *mm = new MyMolecule(ofMmdb);
    // auto myMolecule = std::shared_ptr<MyMolecule>(new MyMolecule::MyMolecule(ofMmdb));
    auto myMolecule = std::shared_ptr<MyMolecule>(mm);
    auto newRepresentationInstance = MolecularRepresentationInstance::create(myMolecule,
                                                                             colorScheme,
                                                                             selectionString,
                                                                             renderStyle);
    camera->getSceneSetup()->addRepresentationInstance(newRepresentationInstance);
    return newRepresentationInstance;
}

void CameraPort::renderSceneSetup()
{
    sceneSetup->renderWithRendererFromViewpointDontChangeMVM(renderer, FCXXCoord (camera->getTranslationX(), camera->getTranslationY(), camera->getTranslationZ()));
}




const std::vector<std::shared_ptr<ColorScheme> > CameraPort::colorSchemes = CameraPort::createColorSchemes();
std::vector<std::shared_ptr<ColorScheme> > CameraPort::createColorSchemes()
{
    //Create a list of color schemes t
    std::vector<std::shared_ptr<ColorScheme> > result;
    //Here the default ones available
    result.push_back(ColorScheme::colorByElementScheme());
    result.push_back(ColorScheme::colorBySecondaryScheme());
    result.push_back(ColorScheme::colorChainsScheme());
    result.push_back(ColorScheme::colorRampChainsScheme());
    
    //And here something a bit more bespoke
    auto bespokeColorScheme = std::shared_ptr<ColorScheme>(new ColorScheme());
    result.push_back(bespokeColorScheme);
    //Known color names are enumerated in SolidColorRule.h
    bespokeColorScheme->addRule(SolidColorRule::colorRuleForSelectionStringAndName("A/", "LAWNGREEN"));
    //Otherwise can use RGBA values
    bespokeColorScheme->addRule(SolidColorRule::colorRuleForSelectionStringAndColor("A/30-60", FCXXCoord(1.f,0.f,0.f,1.f)));
    
    return result;
}

const std::vector<std::string> CameraPort::displayStyles = CameraPort::createDisplayStyles();
std::vector<std::string> CameraPort::createDisplayStyles()
{
    std::vector<std::string> result;
    result.push_back("Ribbon");
    result.push_back("Calpha");
    result.push_back("Sticks");
    result.push_back("Cylinders");
    result.push_back("Spheres");
    result.push_back("MolecularSurface");
    result.push_back("VdWSurface");
    result.push_back("AccessibleSurface");
    return result;
}

CameraPort::CameraPort() : iCurrentColorScheme(0), iCurrentDisplayStyle(0), oldX(-1e30), oldY(-1e30){

    camera = Camera::defaultCamera();
    //Setting the following parameters unnecessary if you are handling the modelview matrix
    camera->setFovy(5);
    camera->setZClipFront(20);
    camera->setFovy(3);

    std::cout << glfwGetVersionString();
    if (!glfwInit()){
        window = 0;
    }
    else {
        glfwWindowHint(GLFW_SAMPLES, 4);
        window = glfwCreateWindow(640, 480, "Hello World", NULL, NULL);
    }
};

void CameraPort::runLoop(std::string message){
    glfwMakeContextCurrent(window);
    renderer->init();
    glfwSetWindowUserPointer(window, (void *)this);
    glfwSetCursorPosCallback(window, CameraPort::cursorMoved);
    glfwSetKeyCallback(window, CameraPort::key_callback);
    
    /* Loop until the user closes the window */
    while (!glfwWindowShouldClose(window))
    {
        /* Render here */
        glClear(GL_COLOR_BUFFER_BIT);
        
        renderer->setupCamera(camera.get());
        renderSceneSetup();
        /* Swap front and back buffers */
        glfwSwapBuffers(window);
        /* Poll for and process events */
        glfwWaitEventsTimeout(1.0);
    }
    
    glfwTerminate();
}

std::shared_ptr<CameraPort> CameraPort::create(){
    std::shared_ptr<CameraPort> result(new CameraPort());
    if (!result->getWindow()){
        glfwTerminate();
        result = 0;
    }
    return result;
};

void CameraPort::zoom(float factor)
{
    camera->setFovy(camera->getFovy()/factor);
}

void CameraPort::cycleColorScheme()
{
    iCurrentColorScheme++;
    if (iCurrentColorScheme>=colorSchemes.size()) iCurrentColorScheme = 0;
    std::shared_ptr<ColorScheme> nextColorScheme = colorSchemes[iCurrentColorScheme];
    
    std::set<std::shared_ptr<MolecularRepresentation> > reps = camera->getSceneSetup()->molecularRepresentations();
    for (auto rep=reps.begin(); rep!= reps.end(); ++rep){
        (*rep)->setColorScheme(nextColorScheme);
    }
}

void CameraPort::cycleDisplayStyle()
{
    iCurrentDisplayStyle++;
    if (iCurrentDisplayStyle>=displayStyles.size()) iCurrentDisplayStyle = 0;
    std::string nextDisplayStyle = displayStyles[iCurrentDisplayStyle];
    
    std::set<std::shared_ptr<MolecularRepresentation> > reps = camera->getSceneSetup()->molecularRepresentations();
    for (auto rep=reps.begin(); rep!= reps.end(); ++rep){
        (*rep)->setRenderStyle(nextDisplayStyle);
    }
}




void CameraPort::cursorMoved(GLFWwindow *win, double newX, double newY){
    auto cameraPort = static_cast<CameraPort*>(glfwGetWindowUserPointer(win));
    float deltaX = float(newX - cameraPort->getOldX());
    float deltaY = float(newY - cameraPort->getOldY());
    float magnitude = sqrtf((deltaX*deltaX) + (deltaY*deltaY));
    int state = glfwGetMouseButton(win, GLFW_MOUSE_BUTTON_LEFT);
    if (state == GLFW_PRESS){
        if (magnitude > 0 && magnitude < 1e15){
            float xComponent =  deltaY / magnitude;
            float yComponent =  deltaX / magnitude;
            FCXXCoord rotation(magnitude/4., xComponent, yComponent, 0.);
            cameraPort->getCamera()->getSceneSetup()->rotateBy(rotation);
        }
    }
    cameraPort->setOldX(newX);
    cameraPort->setOldY(newY);
}

void CameraPort::key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    auto cameraPort = static_cast<CameraPort*>(glfwGetWindowUserPointer(window));
    if (key == GLFW_KEY_N && action == GLFW_PRESS)
        cameraPort->zoom(1.05);
    if (key == GLFW_KEY_M && action == GLFW_PRESS)
        cameraPort->zoom(0.95);
    if (key == GLFW_KEY_C && action == GLFW_PRESS)
        cameraPort->cycleColorScheme();
    if (key == GLFW_KEY_S && action == GLFW_PRESS)
        cameraPort->cycleDisplayStyle();
}
