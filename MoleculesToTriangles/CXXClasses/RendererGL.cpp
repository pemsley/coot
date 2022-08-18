/*
 *  RendererGL.mm
 *  Aesop
 *
 *  Created by Martin Noble on 03/06/2010.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include "RendererGL.h"

//#include "GL/glew.h"

#ifdef HAVE_CONFIG_H
#include "config.h"
#else
#define HAVE_OPENGL_GL_H
#endif

#if defined(HAVE_WINDOWS_H) && defined(_WIN32)
  # include <windows.h>
  #endif
  #ifdef HAVE_GL_GL_H
  # warning Have_GL_GL_H
  # define GL_GLEXT_LEGACY
  # define GL_GLEXT_PROTOTYPES
  # include <GL/gl.h>
  # include <GL/glu.h>
  # include <GL/glext.h>
  #elif defined(HAVE_OPENGL_GL_H)
  # include <OpenGL/gl.h>
  # include <OpenGL/glu.h>
  # include <OpenGL/glext.h>
  # define glGenVertexArrays glGenVertexArraysAPPLE
  # define glDeleteVertexArrays glDeleteVertexArraysAPPLE
  # define glBindVertexArray glBindVertexArrayAPPLE
  #else
  # error no gl.h
  #endif


#include "Camera.h"
#include "SceneSetup.h"
#include "Light.h"
#include "DisplayPrimitive.h"
#include "RepresentationInstance.h"
//#include "LinesPrimitive.h"
#include "BondsPrimitive.h"
#include "SpherePrimitive.h"
#include "oglPolyhedron.h"
#include "VertexColorNormalPrimitive.h"
#include "CXXMatrix.h"
#include "MolecularRepresentation.h"

#ifndef RendererHandles_type
#define RendererHandles_type
typedef struct RendererHandles_ {
    unsigned int vertexHandle;
    unsigned int colorHandle;
    unsigned int indexHandle;
    unsigned int arrayObjectHandle;
} RendererHandles;
#endif

void RendererGL::init(){
    /*
    GLenum err = glewInit();
    if (GLEW_OK != err)
    {
         //Problem: glewInit failed, something is seriously wrong. 
        fprintf(stderr, "Error: %s\n", glewGetErrorString(err));
    }
    fprintf(stdout, "Status: Using GLEW %s\n", glewGetString(GLEW_VERSION));
*/
}

RendererGL::~RendererGL(){
    auto allocatedHandle = allocatedHandles.begin();
    for (; allocatedHandle != allocatedHandles.end(); ++allocatedHandle){
        allocatedHandle->first->removeRenderer(this);
        glDeleteBuffers(1, &(allocatedHandle->second.vertexHandle));
        glDeleteBuffers(1, &(allocatedHandle->second.indexHandle));
        glDeleteVertexArrays(1, &(allocatedHandle->second.arrayObjectHandle));
    }
    allocatedHandles.clear();
}
void RendererGL::myglEnable(int flag){
#ifdef DEBUG_MINE
    std::cout << "In RendererGL::myglEnable\n";
#endif
    ::glEnable(flag);
}
void RendererGL::myglDisable(int flag){
    ::glDisable(flag);
}
void RendererGL::myglClearDepth(float value){
    ::glClearDepth(value);
}
void RendererGL::myglGetFloatv(int propertyToGet, float *dest){
    ::glGetFloatv(propertyToGet, dest);
}
void RendererGL::myglGetIntegerv(int propertyToGet, int *dest)
{
    ::glGetIntegerv(propertyToGet, dest);
}

void RendererGL::drawTestTriangle(const FCXXCoord & at) const
{
    GLuint VertexArrayID;
    glGenVertexArrays(1, &VertexArrayID);
    glBindVertexArray(VertexArrayID);
    FCXXCoord  v1(-1.,-1.,0.);
    FCXXCoord  v2(1.,-1.,0.);
    FCXXCoord  v3(0, 1.,0.);
    // An array of 3 vectors which represents 3 vertices
    static const GLfloat g_vertex_buffer_data[] = {
        static_cast<GLfloat>((v1+at)[0]),static_cast<GLfloat>((v1+at)[1]),static_cast<GLfloat>((v1+at)[2]),
        static_cast<GLfloat>((v2+at)[0]),static_cast<GLfloat>((v2+at)[1]),static_cast<GLfloat>((v2+at)[2]),
        static_cast<GLfloat>((v3+at)[0]),static_cast<GLfloat>((v3+at)[1]),static_cast<GLfloat>((v3+at)[2])
    };
    // This will identify our vertex buffer
    GLuint vertexbuffer;
    // Generate 1 buffer, put the resulting identifier in vertexbuffer
    glGenBuffers(1, &vertexbuffer);
    // The following commands will talk about our 'vertexbuffer' buffer
    glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
    // Give our vertices to OpenGL.
    glBufferData(GL_ARRAY_BUFFER, sizeof(g_vertex_buffer_data), g_vertex_buffer_data, GL_STATIC_DRAW);
    // 1rst attribute buffer : vertices
    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
    glVertexAttribPointer(
                          0,                  // attribute 0. No particular reason for 0, but must match the layout in the shader.
                          3,                  // size
                          GL_FLOAT,           // type
                          GL_FALSE,           // normalized?
                          0,                  // stride
                          (void*)0            // array buffer offset
                          );
    // Draw the triangle !
    glDrawArrays(GL_TRIANGLES, 0, 3); // Starting from vertex 0; 3 vertices total -> 1 triangle
    glDisableVertexAttribArray(0);
}

void RendererGL::resize(int w, int h)
{
    glViewport(0, 0, w, h);
}

std::shared_ptr<Renderer> RendererGL::create()
{
    return std::shared_ptr<Renderer>(new RendererGL());
}



#define RendererType RendererGL

/*
 *  RendererImplementations.h
 *  AesopCD
 *
 *  Created by Martin Noble on 08/11/2010.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

void RendererType::setupCamera(Camera *camera){
#ifdef DEBUG_MINE
    std::cout << "In setupCamera\n";
#endif
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    CXXMatrix lookAtMatrix(CXXMatrix::lookAtMatrix(camera->getRotation(), camera->getTranslation()));
    GLfloat lookAtFloats[16];
    lookAtMatrix.gLOperator(lookAtFloats);
    glMultMatrixf(lookAtFloats);
    
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    
    double fovy = camera->getFovy();//(180./M_PI) * 20. * 2. * atan2(camera->getSceneSetup()->getScale(), eyeDistance);
    double cameraNear = camera->getZClipFront();
    double cameraFar = cameraNear - camera->getZClipDepth();
    double eyeDistance = camera->getTranslation().get3DLength();
    double zNear = eyeDistance - cameraNear*camera->getSceneSetup()->getScale();
    zNear = (zNear>0.01?zNear:0.01);
    double zFar = eyeDistance - cameraFar*camera->getSceneSetup()->getScale();
    zFar = (zFar>0.02?zFar:0.02);
    GLint viewportParams[6];
    glGetIntegerv(GL_VIEWPORT, viewportParams);
    double aspectRatio = (float)(viewportParams[2])/(float)(viewportParams[3]);
    
    gluPerspective(fovy, aspectRatio, zNear, zFar);
    //glOrtho(-aspectRatio, aspectRatio, -1.0f, 1.0f, -100.f, 100.f);
    
    if (!camera->getSceneSetup()){
        throw "Camera trying to render with out having an assigned scene setup";
    }
    
    
    //Scene fog
    FCXXCoord backgroundColor(camera->getSceneSetup()->getBackgroundColor());
    glClearColor(backgroundColor[0], backgroundColor[1], backgroundColor[2], backgroundColor[3]);
    myglClearDepth(1.);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    GLfloat fogColor[] = {static_cast<GLfloat>(backgroundColor[0]), static_cast<GLfloat>(backgroundColor[1]), static_cast<GLfloat>(backgroundColor[2]), static_cast<GLfloat>(backgroundColor[3])};
    
    glFogfv(GL_FOG_COLOR, fogColor);
    glFogi(GL_FOG_MODE, GL_LINEAR);
    float cameraDistance = camera->getTranslation().get3DLength();
    float fogStart = cameraDistance + (camera->getFogFront()  * camera->getSceneSetup()->getScale());
    float fogEnd = fogStart + (camera->getFogDepthRange() * camera->getSceneSetup()->getScale());
    if (fogStart < 0.) fogStart = 0.;
    if (fogEnd < 1.) fogEnd = 1.;
    glFogf(GL_FOG_START, fogStart);
    glFogf(GL_FOG_END, fogEnd);
    
}


void RendererType::setupScene(SceneSetup *sceneSetup){
#ifdef DEBUG_MINE
    std::cout << "In setupScene\n";
#endif
    myglEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);
    
    //Make ready for lighting
    myglEnable(GL_LIGHTING);
    
    // Overall modelView matrix
    FCXXCoord cxxRotation(sceneSetup->getRotation());
    glMatrixMode(GL_MODELVIEW);
    glRotatef(cxxRotation[0], cxxRotation[1], cxxRotation[2], cxxRotation[3]);
    float cScale = sceneSetup->getScale();
    glScalef (cScale, cScale, cScale);
    FCXXCoord cxxTranslation(sceneSetup->getTranslation());
    glTranslatef(cxxTranslation[0], cxxTranslation[1], cxxTranslation[2]);
    drawTestTriangle(cxxTranslation*-1.);
    /*
    std::cout << "Camera is " << camera << std::endl;
    FCXXCoord cameraTranslation = camera->getTranslation();
    float cameraDistance = cameraTranslation.get3DLength();
    float fogStart = cameraDistance + (camera->getFogFront()  * camera->getSceneSetup()->getScale());
    float fogEnd = fogStart + (camera->getFogDepthRange() * camera->getSceneSetup()->getScale());
    if (fogStart < 0.) fogStart = 0.;
    if (fogEnd < 1.) fogEnd = 1.;
 
    //glFogi(GL_FOG_COORD_SRC, GL_FOG_COORDINATE);
    GLfloat fogColor[] = {1.,0.,0.,0.};
    glFogfv(GL_FOG_COLOR, fogColor);
    glFogi(GL_FOG_MODE, GL_LINEAR);
    glFogf(GL_FOG_START, -20);
    glFogf(GL_FOG_END, 20);
*/
}

void RendererType::setupLightAsIndexFromViewpointWithScale(Light *light, int asIndex, FCXXCoord fromViewpoint, float scale){
    glMatrixMode(GL_MODELVIEW);
    /*
     GLfloat params[16];
    glGetFloatv(GL_MODELVIEW_MATRIX, params);
    for (int i=0; i<16; i++) std::cout << params[i] << " ";
    std::cout << std::endl;
     */
    glPushMatrix();
    glLoadIdentity();
    
    //Apply "look At" and scale matrices
    gluLookAt(fromViewpoint[0], fromViewpoint[1], fromViewpoint[2],//0., 0., dist,
              0., 0., 0.,
              0., 1., 0.);
    glScalef(scale, scale, scale);
    
    GLenum lightEnum = GL_LIGHT0 + asIndex;
    GLfloat lightPosition[] = {
        static_cast<GLfloat>(light->getTranslation()[0]),
        static_cast<GLfloat>(light->getTranslation()[1]),
        static_cast<GLfloat>(light->getTranslation()[2]),
        0.
    };
    
    glLightf(lightEnum, GL_LINEAR_ATTENUATION, 0.0);
    if (light->getLightType() == Light::Positional){
        glLightf(lightEnum, GL_QUADRATIC_ATTENUATION, 1.0f/(scale*scale));
        lightPosition[3] = 1.;
    }
    else {
        glLightf(lightEnum, GL_QUADRATIC_ATTENUATION, 0.0);
        lightPosition[3] = 0.;
    }
    glLightfv (lightEnum, GL_POSITION, lightPosition);
    //std::cout << "Light " << lightEnum << " at " << lightPosition[0]<<" " << lightPosition[1]<<" " << lightPosition[2]<<std::endl;
    {
        FCXXCoord intensityAdjustedColor(light->getAmbient() * light->getIntensity());
        GLfloat glLightColor[] = {
            static_cast<GLfloat>(intensityAdjustedColor[0]),
            static_cast<GLfloat>(intensityAdjustedColor[1]),
            static_cast<GLfloat>(intensityAdjustedColor[2]),
            1.};
        glLightfv (lightEnum, GL_AMBIENT, glLightColor);
    }
    {
        FCXXCoord intensityAdjustedColor(light->getDiffuse() * light->getIntensity());
        GLfloat glLightColor[] = {
            static_cast<GLfloat>(intensityAdjustedColor[0]),
            static_cast<GLfloat>(intensityAdjustedColor[1]),
            static_cast<GLfloat>(intensityAdjustedColor[2]),
            1.};
        glLightfv (lightEnum, GL_DIFFUSE, glLightColor);
    }
    {
        FCXXCoord intensityAdjustedColor(light->getSpecular() * light->getIntensity());
        GLfloat glLightColor[] = {
            static_cast<GLfloat>(intensityAdjustedColor[0]),
            static_cast<GLfloat>(intensityAdjustedColor[1]),
            static_cast<GLfloat>(intensityAdjustedColor[2]),
            1.};
        glLightfv (lightEnum, GL_SPECULAR, glLightColor);
    }
    myglEnable(lightEnum);
    
    if (light->getDrawLight()){
        FCXXCoord lightColor(1.,1.,1.,1.);
        SpherePrimitive lightBall(FCXXCoord (lightPosition[0], lightPosition[1],lightPosition[2], 1.), 0.2, lightColor);
        float glowColor[] = {0.8, 0.8, 0.65, 1.0};
        float blackColor[] = {1.0, 1.0, 1.0, 1.0};
        glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, glowColor);
        lightBall.renderWithRenderer(std::shared_ptr<Renderer>(this));
        glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, blackColor);
    }
    
    glPopMatrix();
}

void RendererType::setupRepresentationInstance(RepresentationInstance *instance){
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    GLfloat matrix[16];
    glGetFloatv (GL_MODELVIEW_MATRIX, matrix);
    glLoadIdentity();
    
    //Instance post transformation
    {
        auto postTransformation = instance->getPostTransformation();
        if (postTransformation) {
            FCXXCoord cxxTranslation(instance->getPostTransformation()->getTranslation());
            glTranslatef(cxxTranslation[0], cxxTranslation[1], cxxTranslation[2]);
            float cScale = instance->getPostTransformation()->getScale();
            glScalef (cScale, cScale, cScale);
            FCXXCoord cxxRotation(instance->getPostTransformation()->getRotation());
            glRotatef(cxxRotation[0], cxxRotation[1], cxxRotation[2], cxxRotation[3]);
        }
    }
    glMultMatrixf(matrix);
    
    /*
     NSArray *myActions = [self getOrderedActions];
     Action *action;
     for (action in myActions){
     if ([action isMemberOfClass:[RotationAction class]]){
     RotationAction *rotationAction = (RotationAction *)action;
     GLfloat localRotation[4];
     [rotationAction glRotation:localRotation atTime:[NSDate timeIntervalSinceReferenceDate]];
     glRotatef(localRotation[0], localRotation[1], localRotation[2], localRotation[3]);
     }
     }
     */
    FCXXCoord cxxRotation(instance->getRotation());
    glRotatef(cxxRotation[0], cxxRotation[1], cxxRotation[2], cxxRotation[3]);
    
    float cScale = instance->getScale();
    glScalef (cScale, cScale, cScale);
    
    FCXXCoord cxxTranslation(instance->getTranslation());
    glTranslatef(cxxTranslation[0], cxxTranslation[1], cxxTranslation[2]);
}

void RendererType::restoreModelviewMatrix() {
    //Note here we need to restore previous modelview matrix:  may need thought for
    //other renderers
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();
}

/*void RendererType::renderLinesPrimitive(LinesPrimitive *lines) {
 glLineWidth (1.0);
 myglDisable(GL_LIGHTING);
 myglDisable(GL_TEXTURE_2D);
 
 glEnableClientState(GL_VERTEX_ARRAY);
 glEnableClientState(GL_COLOR_ARRAY);
 
 glColorPointer(4, GL_FLOAT, 0, lines->getColorArray());
 glVertexPointer(3, GL_FLOAT, 0, lines->getVertexArray());
 glDrawArrays(GL_LINES, 0, lines->getNBonds());
 
 glDisableClientState(GL_VERTEX_ARRAY);
 glDisableClientState(GL_COLOR_ARRAY);
 
 myglEnable(GL_LIGHTING);
 }*/

void RendererType::renderPolyhedron(oglPolyhedron *prim) {
    GLint oldMatrixMode;
    myglGetIntegerv(GL_MATRIX_MODE, &oldMatrixMode);
    if (oldMatrixMode != GL_MODELVIEW){
        glMatrixMode(GL_MODELVIEW);
    }
    glPushMatrix();
    float compositeScale = prim->getRadius();
    glTranslatef(prim->getCentre()[0], prim->getCentre()[1], prim->getCentre()[2]);
    glScalef(compositeScale, compositeScale, compositeScale);
    myglEnable(GL_NORMALIZE);
    GLfloat matrix[16];
    myglGetFloatv(GL_MODELVIEW_MATRIX, matrix);
    
    //Specify material properties that should be taken from the color
    //glColorMaterial not available in OpenGL ES
    //glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
    
    //Set material properties that are not per-vertex
    GLfloat specularColor[] = {1., 1., 1., 1.};
    GLfloat blackColor[] = {0., 0., 0., 1.};
    GLfloat polyhedronColor[] = {static_cast<GLfloat>(prim->getColor()[0]), static_cast<GLfloat>(prim->getColor()[1]),
        static_cast<GLfloat>(prim->getColor()[2]),static_cast<GLfloat>(prim->getColor()[3])};
    
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specularColor);
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, blackColor);
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, polyhedronColor);
    glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 32.0);
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, blackColor);
    
    vboRenderVN(prim);
    /*
     glEnableClientState(GL_VERTEX_ARRAY);
     glEnableClientState(GL_NORMAL_ARRAY);
     glDisableClientState(GL_COLOR_ARRAY);
     
     GLvoid *vertexStart = (GLvoid *) &(prim->getVertexColorNormalArray()[0].vertex[0]);
     glVertexPointer(3, GL_FLOAT, sizeof(VertexColorNormal), vertexStart);
     GLvoid *colorStart = (GLvoid *) &(prim->getVertexColorNormalArray()[0].color[0]);
     glColorPointer(4, GL_FLOAT, sizeof(VertexColorNormal), colorStart);
     GLvoid *normalStart = (GLvoid *) &(prim->getVertexColorNormalArray()[0].normal[0]);
     glNormalPointer(GL_FLOAT, sizeof(VertexColorNormal), normalStart);
     
     glDrawElements(GL_TRIANGLES, 3*prim->nTriangles(), kGLIndexType, prim->getIndexArray());
     
     glDisableClientState(GL_VERTEX_ARRAY);
     glDisableClientState(GL_NORMAL_ARRAY);
     */
    
    myglDisable(GL_COLOR_MATERIAL);
    
    glPopMatrix();
    //Restore modelview matrix mode and modelview matrx
    if (oldMatrixMode != GL_MODELVIEW){
        glMatrixMode(oldMatrixMode);
    }
}

void RendererType::renderVertexColorNormalPrimitive(VertexColorNormalPrimitive *prim) {
#ifdef DEBUG_MINE
    std::cout << "In renderVertexColorNormalPrimitiveA\n";
#endif
    
    myglEnable(GL_NORMALIZE);
    
    //Specify material properties that should be taken from the color
    //glColorMaterial not available in OpenGL ES
    //glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
    myglEnable(GL_COLOR_MATERIAL);
    
    //Set material properties that are not per-vertex
    GLfloat specularColor[] = {1., 1., 1., 1.};
    GLfloat blackColor[] = {0., 0., 0., 1.};
    
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specularColor);
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, blackColor);
    glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, blackColor);
    glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 32.0);
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, blackColor);
    vboRenderVCN(prim);
    /*
     glEnableClientState(GL_VERTEX_ARRAY);
     glEnableClientState(GL_NORMAL_ARRAY);
     glEnableClientState(GL_COLOR_ARRAY);
     
     GLvoid *vertexStart = (GLvoid *) &(prim->getVertexColorNormalArray()[0].vertex);
     glVertexPointer(3, GL_FLOAT, sizeof(VertexColorNormalPrimitive::VertexColorNormal), vertexStart);
     GLvoid *colorStart = (GLvoid *) &(prim->getVertexColorNormalArray()[0].color);
     glColorPointer(4, GL_FLOAT, sizeof(VertexColorNormalPrimitive::VertexColorNormal), colorStart);
     GLvoid *normalStart = (GLvoid *) &(prim->getVertexColorNormalArray()[0].normal);
     glNormalPointer(GL_FLOAT, sizeof(VertexColorNormalPrimitive::VertexColorNormal), normalStart);
     glDrawElements(GL_TRIANGLES, 3*prim->nTriangles(), kGLIndexType, prim->getIndexArray());
     
     glDisableClientState(GL_VERTEX_ARRAY);
     glDisableClientState(GL_NORMAL_ARRAY);
     glDisableClientState(GL_COLOR_ARRAY);
     
     */
    myglDisable(GL_COLOR_MATERIAL);
}

void RendererType::renderVertexColorPrimitive(VertexColorNormalPrimitive *prim) {
    
    glLineWidth (1.0);
    myglDisable(GL_LIGHTING);
    myglDisable(GL_TEXTURE_2D);
    
    vboRenderVC(prim);
    
    myglEnable(GL_LIGHTING);
    /*
     glEnableClientState(GL_VERTEX_ARRAY);
     glEnableClientState(GL_NORMAL_ARRAY);
     glEnableClientState(GL_COLOR_ARRAY);
     
     GLvoid *vertexStart = (GLvoid *) &(prim->getVertexColorNormalArray()[0].vertex);
     glVertexPointer(3, GL_FLOAT, sizeof(VertexColorNormal), vertexStart);
     GLvoid *colorStart = (GLvoid *) &(prim->getVertexColorNormalArray()[0].color);
     glColorPointer(4, GL_FLOAT, sizeof(VertexColorNormal), colorStart);
     GLvoid *normalStart = (GLvoid *) &(prim->getVertexColorNormalArray()[0].normal);
     glNormalPointer(GL_FLOAT, sizeof(VertexColorNormal), normalStart);
     glDrawElements(GL_TRIANGLES, 3*prim->nTriangles(), kGLIndexType, prim->getIndexArray());
     
     glDisableClientState(GL_VERTEX_ARRAY);
     glDisableClientState(GL_NORMAL_ARRAY);
     glDisableClientState(GL_COLOR_ARRAY);
     */
    
}

FCXXCoord RendererType::angstromsForPixels(float x, float y){
    GLfloat modelviewMatrix[16];
    myglGetFloatv(GL_MODELVIEW_MATRIX, modelviewMatrix);
    CXXMatrix totalTransformation(modelviewMatrix);
    CXXMatrix invertedTotalTransformation = totalTransformation.inverseGLOperator();
    
    //Remove translation from the resulting matrix
    for (int i=0; i<3; i++) invertedTotalTransformation[3][i] = 0.;
    invertedTotalTransformation[3][3] = 1.;
    
    GLint viewportParams[6];
    glGetIntegerv(GL_VIEWPORT, viewportParams);
    FCXXCoord oldVector((2.*(double)x) /(double)viewportParams[2], -(2.*(double)y)/(double)viewportParams[3], 0., 1.);
    
    FCXXCoord additionalTranslation = invertedTotalTransformation * oldVector;
    return additionalTranslation;
}

void RendererType::vboRenderVCN(VertexColorNormalPrimitive *prim)
{
#ifdef DEBUG_MINE
    std::cout << "In RendererGL::vboRencerVCN render " << prim->nVertices();
#endif
    RendererHandles rendererHandles;
    std::map<DisplayPrimitive *, RendererHandles>::iterator mapIter = allocatedHandles.find(prim);
    if (mapIter != allocatedHandles.end()){
        rendererHandles = mapIter->second;
    }
    else {
        //std::cout << "Allocating handles for primitive" << prim << std::endl;
        GLuint vertexHandle, indexHandle, arrayObjectHandle;
        
        glGenVertexArrays(1, &arrayObjectHandle);
        glBindVertexArray(arrayObjectHandle);
        rendererHandles.arrayObjectHandle = arrayObjectHandle;
        
        glGenBuffers(1, &vertexHandle);
        rendererHandles.vertexHandle = vertexHandle;
        
        glBindBuffer(GL_ARRAY_BUFFER, vertexHandle);
        glBufferData(GL_ARRAY_BUFFER, prim->nVertices()*sizeof(VertexColorNormalPrimitive::VertexColorNormal),
                     prim->getVertexColorNormalArray(), GL_STATIC_DRAW);
        
        glEnableClientState(GL_VERTEX_ARRAY);
        glVertexPointer(3, GL_FLOAT, sizeof(VertexColorNormalPrimitive::VertexColorNormal),
                        (void*)offsetof(VertexColorNormalPrimitive::VertexColorNormal,vertex));
        
        glEnableClientState(GL_COLOR_ARRAY);
        glColorPointer(4, GL_FLOAT, sizeof(VertexColorNormalPrimitive::VertexColorNormal),
                       (void*)offsetof(VertexColorNormalPrimitive::VertexColorNormal,color));
        
        glEnableClientState(GL_NORMAL_ARRAY);
        glNormalPointer(GL_FLOAT, sizeof(VertexColorNormalPrimitive::VertexColorNormal),
                        (void*)offsetof(VertexColorNormalPrimitive::VertexColorNormal,normal));
        
        glGenBuffers(1, &indexHandle);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexHandle);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, prim->nTriangles() * 3 * sizeof (GLIndexType),
                     prim->getIndexArray(), GL_STATIC_DRAW);
        
        rendererHandles.vertexHandle = vertexHandle;
        rendererHandles.indexHandle = indexHandle;
        
        prim->addRenderer(this);
        allocatedHandles[prim] = rendererHandles;
        
        glBindVertexArray(0);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
        
        glDisableClientState(GL_NORMAL_ARRAY);
        glDisableClientState(GL_COLOR_ARRAY);
        glDisableClientState(GL_VERTEX_ARRAY);
    }
    
    glBindVertexArray(rendererHandles.arrayObjectHandle);
    int nTriangleVertices = 3*prim->nTriangles();
    glDrawElements(GL_TRIANGLES, nTriangleVertices, kGLIndexType, (void*)0);
    glBindVertexArray(0);
}

void RendererType::vboRenderVCNFixed(VertexColorNormalPrimitive *prim)
{
    RendererHandles rendererHandles;
    std::map<DisplayPrimitive *, RendererHandles>::iterator mapIter = allocatedHandles.find(prim);
    if (mapIter != allocatedHandles.end()){
        rendererHandles = mapIter->second;
    }
    else {
        GLuint vertexHandle, indexHandle, arrayObjectHandle;
        
        glGenBuffers(1, &vertexHandle);
        rendererHandles.vertexHandle = vertexHandle;
        
        glBindBuffer(GL_ARRAY_BUFFER, vertexHandle);
        glBufferData(GL_ARRAY_BUFFER, prim->nVertices()*sizeof(VertexColorNormalPrimitive::VertexColorNormal),
                     prim->getVertexColorNormalArray(), GL_STATIC_DRAW);
        
        glGenVertexArrays(1, &arrayObjectHandle);
        rendererHandles.arrayObjectHandle = arrayObjectHandle;
        
        glBindVertexArray(arrayObjectHandle);
        
        glBindBuffer(GL_ARRAY_BUFFER, vertexHandle);
        
        glEnableClientState(GL_VERTEX_ARRAY);
        glVertexPointer(3, GL_FLOAT, sizeof(VertexColorNormalPrimitive::VertexColorNormal),
                        (void*)offsetof(VertexColorNormalPrimitive::VertexColorNormal,vertex));
        
        glEnableClientState(GL_COLOR_ARRAY);
        glColorPointer(4, GL_UNSIGNED_BYTE, sizeof(VertexColorNormalPrimitive::VertexColorNormal),
                       (void*)offsetof(VertexColorNormalPrimitive::VertexColorNormal,color));
        
        glEnableClientState(GL_NORMAL_ARRAY);
        glNormalPointer(GL_FLOAT, sizeof(VertexColorNormalPrimitive::VertexColorNormal),
                        (void*)offsetof(VertexColorNormalPrimitive::VertexColorNormal,normal));
        
        glGenBuffers(1, &indexHandle);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexHandle);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, prim->nTriangles() * 3 * sizeof (GLIndexType),
                     prim->getIndexArray(), GL_STATIC_DRAW);
        
        rendererHandles.vertexHandle = vertexHandle;
        rendererHandles.indexHandle = indexHandle;
        
        prim->addRenderer(this);
        allocatedHandles[prim] = rendererHandles;
        
        glBindVertexArray(0);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
        
        glDisableClientState(GL_NORMAL_ARRAY);
        glDisableClientState(GL_COLOR_ARRAY);
        glDisableClientState(GL_VERTEX_ARRAY);
    }
    
    glBindVertexArray(rendererHandles.arrayObjectHandle);
    int nTriangleVertices = 3*prim->nTriangles();
    glDrawElements(GL_TRIANGLES, nTriangleVertices, kGLIndexType, (void*)0);
    
    glBindVertexArray(0);
}

void RendererType::vboRenderVC(VertexColorNormalPrimitive *prim)
{
    RendererHandles rendererHandles;
    std::map<DisplayPrimitive *, RendererHandles>::iterator mapIter = allocatedHandles.find(prim);
    if (mapIter != allocatedHandles.end()){
        rendererHandles = mapIter->second;
    }
    else {
        GLuint vertexHandle, indexHandle;
        
        glGenBuffers(1, &vertexHandle);
        glGenBuffers(1, &indexHandle);
        
        glBindBuffer(GL_ARRAY_BUFFER, vertexHandle);
        glBufferData(GL_ARRAY_BUFFER, prim->nVertices()*sizeof(VertexColorNormalPrimitive::VertexColor), prim->getVertexColorArray(), GL_STATIC_DRAW);
        
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexHandle);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, prim->nLines() * 2 * sizeof (GLIndexType),
                     prim->getIndexArray(), GL_STATIC_DRAW);
        
        rendererHandles.vertexHandle = vertexHandle;
        rendererHandles.indexHandle = indexHandle;
        
        prim->addRenderer(this);
        allocatedHandles[prim] = rendererHandles;
    }
    
    glBindBuffer(GL_ARRAY_BUFFER, rendererHandles.vertexHandle);
    glEnableClientState(GL_VERTEX_ARRAY);
    glVertexPointer(3, GL_FLOAT, sizeof(VertexColorNormalPrimitive::VertexColor),
                    (void*)offsetof(VertexColorNormalPrimitive::VertexColor,vertex));
    glEnableClientState(GL_COLOR_ARRAY);
    glColorPointer(4, GL_FLOAT, sizeof(VertexColorNormalPrimitive::VertexColor),
                   (void*)offsetof(VertexColorNormalPrimitive::VertexColor,color));
    
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, rendererHandles.indexHandle);
    int a = prim->nLines();
    glDrawElements(GL_LINES, 2*prim->nLines(), kGLIndexType, (void*)0);
    
    
    glBindBuffer(GL_ARRAY_BUFFER,0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);
    
    glDisableClientState(GL_COLOR_ARRAY);
    glDisableClientState(GL_VERTEX_ARRAY);
    
}


void RendererType::vboRenderVN(VertexColorNormalPrimitive *prim)
{
    RendererHandles rendererHandles;
    std::map<DisplayPrimitive *, RendererHandles>::iterator mapIter = allocatedHandles.find(prim);
    if (mapIter != allocatedHandles.end()){
        rendererHandles = mapIter->second;
    }
    else {
        GLuint vertexHandle, indexHandle;
        
        glGenBuffers(1, &vertexHandle);
        glGenBuffers(1, &indexHandle);
        
        glBindBuffer(GL_ARRAY_BUFFER, vertexHandle);
        glBufferData(GL_ARRAY_BUFFER, prim->nVertices()*sizeof(VertexColorNormalPrimitive::VertexNormal),
                     prim->getVertexNormalArray(), GL_STATIC_DRAW);
        
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexHandle);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, prim->nTriangles() * 3 * sizeof (GLIndexType),
                     prim->getIndexArray(), GL_STATIC_DRAW);
        
        rendererHandles.vertexHandle = vertexHandle;
        rendererHandles.indexHandle = indexHandle;
        
        prim->addRenderer(this);
        allocatedHandles[prim] = rendererHandles;
    }
    
    glBindBuffer(GL_ARRAY_BUFFER, rendererHandles.vertexHandle);
    glEnableClientState(GL_VERTEX_ARRAY);
    glVertexPointer(3, GL_FLOAT, sizeof(VertexColorNormalPrimitive::VertexNormal),
                    (void*)offsetof(VertexColorNormalPrimitive::VertexNormal,vertex));
    glEnableClientState(GL_NORMAL_ARRAY);
    glNormalPointer(GL_FLOAT, sizeof(VertexColorNormalPrimitive::VertexNormal),
                    (void*)offsetof(VertexColorNormalPrimitive::VertexNormal,normal));
    
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, rendererHandles.indexHandle);
    glDrawElements(GL_TRIANGLES, 3*prim->nTriangles(), kGLIndexType, (void*)0);
    
    glBindBuffer(GL_ARRAY_BUFFER,0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);
    
    glDisableClientState(GL_NORMAL_ARRAY);
    glDisableClientState(GL_VERTEX_ARRAY);
}

void RendererType::liberateHandlesForDisplayPrimitive(DisplayPrimitive *prim)  {
    //std::cout << "In liberateHandlesForDisplayPrimitive\n";
    std::map<DisplayPrimitive *, RendererHandles>::iterator rendererHandles = allocatedHandles.find(prim);
    if (rendererHandles != allocatedHandles.end()){
        //std::cout << "Doing the erase\n";
        glDeleteBuffers(1, &rendererHandles->second.vertexHandle);
        glDeleteBuffers(1, &rendererHandles->second.indexHandle);
        glDeleteVertexArrays(1, &(rendererHandles->second.arrayObjectHandle));
        allocatedHandles.erase(rendererHandles);
    }
    rendererHandles = allocatedHandles.find(prim);
    //std::cout << "In list is now " << (rendererHandles == allocatedHandles.end()) << std::endl;
}

void RendererType::renderBondsPrimitive(BondsPrimitive *prim)
{
    glLineWidth (1.0);
    myglDisable(GL_LIGHTING);
    myglDisable(GL_TEXTURE_2D);
    
    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_COLOR_ARRAY);
    
    BondsPrimitive::VertexColor *vertexColorArray = prim->getVertexColorArray();
    glColorPointer(4, GL_FLOAT, sizeof(BondsPrimitive::VertexColor), &(vertexColorArray[0].color));
    glVertexPointer(3, GL_FLOAT, sizeof(BondsPrimitive::VertexColor), &(vertexColorArray[0].vertex));
    glDrawElements(GL_LINES, 2*prim->getNBonds(), kGLIndexType, prim->getIndexArray());
    glDisableClientState(GL_VERTEX_ARRAY);
    glDisableClientState(GL_COLOR_ARRAY);
    
    myglEnable(GL_LIGHTING);
}

void RendererType::render(Camera *camera) {
    FCXXCoord backgroundColor(camera->getSceneSetup()->getBackgroundColor());
    // Clear the Canvas
    glClearColor(backgroundColor[0], backgroundColor[1], backgroundColor[2], backgroundColor[3]);
    myglClearDepth(1.);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    myglEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);
    //Make ready for lighting
    myglEnable(GL_LIGHTING);
    
    //gluLookAt stuff
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    CXXMatrix lookAtMatrix(CXXMatrix::lookAtMatrix(camera->getRotation(), camera->getTranslation()));
    GLfloat lookAtFloats[16];
    lookAtMatrix.gLOperator(lookAtFloats);
    glMultMatrixf(lookAtFloats);
    //projection
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    
    GLint viewportParams[6];
    glGetIntegerv(GL_VIEWPORT, viewportParams);
    double aspectRatio = (float)(viewportParams[2])/(float)(viewportParams[3]);
    double eyeDistance = camera->getTranslation().get3DLength();
    double fovy = camera->getFovy();//(180./M_PI) * 2. * atan2(camera->getSceneSetup()->getScale(), eyeDistance);
    double cameraNear = camera->getZClipFront();
    double cameraFar = cameraNear - camera->getZClipDepth();
    double zNear = eyeDistance - cameraNear*camera->getSceneSetup()->getScale();
    zNear = (zNear>0.01?zNear:0.01);
    double zFar = eyeDistance - cameraFar*camera->getSceneSetup()->getScale();
    zFar = (zFar>0.02?zFar:0.02);
    gluPerspective(fovy, aspectRatio, zNear, zFar);
    
    //Make ready for lighting
    myglEnable(GL_LIGHTING);
    GLint a = GL_TRUE;
    glLightModeliv(GL_LIGHT_MODEL_TWO_SIDE, &a);
    std::shared_ptr<SceneSetup> sceneSetup = camera->getSceneSetup();
    //Lights
    auto lightsPntr = sceneSetup->lightsBegin();
    auto lightsEnd = sceneSetup->lightsEnd();
    for (int i=0; i<3; i++){
        GLenum lightEnum = GL_LIGHT0 + i;
        myglDisable(lightEnum);
    }
    int asIndex = 0;
    for (; lightsPntr != lightsEnd; ++lightsPntr){
        auto light(*lightsPntr);
        float scale = camera->getSceneSetup()->getScale();
        glMatrixMode(GL_MODELVIEW);
        glPushMatrix();
        glLoadIdentity();
        
        glMultMatrixf(lookAtFloats);
        
        glScalef(scale, scale, scale);
        
        GLenum lightEnum = GL_LIGHT0 + asIndex;
        GLfloat lightPosition[] = {
            static_cast<GLfloat>(light->getTranslation()[0]),
            static_cast<GLfloat>(light->getTranslation()[1]),
            static_cast<GLfloat>(light->getTranslation()[2]),
            0.
        };
        
        glLightf(lightEnum, GL_LINEAR_ATTENUATION, 0.0);
        if (light->getLightType() == Light::Positional){
            glLightf(lightEnum, GL_QUADRATIC_ATTENUATION, 1.0f);
            lightPosition[3] = 1.;
        }
        else {
            glLightf(lightEnum, GL_QUADRATIC_ATTENUATION, 0.0);
            lightPosition[3] = 0.;
        }
        glLightfv (lightEnum, GL_POSITION, lightPosition);
        {
            FCXXCoord intensityAdjustedColor(light->getAmbient() * light->getIntensity());
            GLfloat glLightColor[] = {
                static_cast<GLfloat>(intensityAdjustedColor[0]),
                static_cast<GLfloat>(intensityAdjustedColor[1]),
                static_cast<GLfloat>(intensityAdjustedColor[2]),
                1.};
            glLightfv (lightEnum, GL_AMBIENT, glLightColor);
        }
        {
            FCXXCoord intensityAdjustedColor(light->getDiffuse() * light->getIntensity());
            GLfloat glLightColor[] = {
                static_cast<GLfloat>(intensityAdjustedColor[0]),
                static_cast<GLfloat>(intensityAdjustedColor[1]),
                static_cast<GLfloat>(intensityAdjustedColor[2]),
                1.};
            glLightfv (lightEnum, GL_DIFFUSE, glLightColor);
        }
        {
            FCXXCoord intensityAdjustedColor(light->getSpecular() * light->getIntensity());
            GLfloat glLightColor[] = {
                static_cast<GLfloat>(intensityAdjustedColor[0]),
                static_cast<GLfloat>(intensityAdjustedColor[1]),
                static_cast<GLfloat>(intensityAdjustedColor[2]),
                1.};
            glLightfv (lightEnum, GL_SPECULAR, glLightColor);
        }
        myglEnable(lightEnum);
        
        if (light->getDrawLight()){
            FCXXCoord lightColor(1.,1.,1.,1.);
            SpherePrimitive lightBall(FCXXCoord (lightPosition[0], lightPosition[1],lightPosition[2], 1.), 0.2, lightColor);
            float glowColor[] = {0.8, 0.1, 0.1, 1.0};
            float blackColor[] = {1.0, 1.0, 1.0, 1.0};
            glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, glowColor);
            lightBall.renderWithRenderer(std::shared_ptr<RendererGL>(this));
            glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, blackColor);
        }
        
        glPopMatrix();
        asIndex++;
    }
    // Representation instances
    auto representationInstancesPntr = sceneSetup->representationInstancesBegin();
    auto representationInstancesEnd = sceneSetup->representationInstancesEnd();
    for (; representationInstancesPntr != representationInstancesEnd; ++representationInstancesPntr){
        auto instance(*representationInstancesPntr);
        if (instance->getDoDraw()){
            std::shared_ptr<Representation>representation = instance->getRepresentation();
            if (representation->getDoDraw()){
                glMatrixMode(GL_MODELVIEW);
                glPushMatrix();
                
                glLoadIdentity();
                glMultMatrixf(lookAtFloats);
                
                //Instance post transformation
                {
                    FCXXCoord cxxRotation(instance->getPostTransformation()->getRotation());
                    glRotatef(cxxRotation[0], cxxRotation[1], cxxRotation[2], cxxRotation[3]);
                    
                    float cScale = instance->getPostTransformation()->getScale();
                    glScalef (cScale, cScale, cScale);
                    
                    FCXXCoord cxxTranslation(instance->getPostTransformation()->getTranslation());
                    glTranslatef(cxxTranslation[0], cxxTranslation[1], cxxTranslation[2]);
                    
                }
                // Overall modelView matrix
                {
                    FCXXCoord cxxRotation(sceneSetup->getRotation());
                    glRotatef(cxxRotation[0], cxxRotation[1], cxxRotation[2], cxxRotation[3]);
                    float cScale = sceneSetup->getScale();
                    glScalef (cScale, cScale, cScale);
                    FCXXCoord cxxTranslation(sceneSetup->getTranslation());
                    glTranslatef(cxxTranslation[0], cxxTranslation[1], cxxTranslation[2]);
                }
                //Instance pre transformation
                {
                    FCXXCoord cxxRotation(instance->getRotation());
                    glRotatef(cxxRotation[0], cxxRotation[1], cxxRotation[2], cxxRotation[3]);
                    
                    float cScale = instance->getScale();
                    glScalef (cScale, cScale, cScale);
                    
                    FCXXCoord cxxTranslation(instance->getTranslation());
                    glTranslatef(cxxTranslation[0], cxxTranslation[1], cxxTranslation[2]);
                }
                instance->getRepresentation()->renderWithRenderer(std::shared_ptr<RendererGL>(this));
                
                glPopMatrix();
            }
        }
    }
    //Lights
    lightsPntr = sceneSetup->lightsBegin();
    asIndex = 0;
    for (; lightsPntr != lightsEnd; ++lightsPntr){
        GLenum lightEnum = GL_LIGHT0 + asIndex;
        myglDisable(lightEnum);
    }
}

void RendererType::clearCameraCanvas(Camera *camera)
{
    //Scene setup
    std::shared_ptr<SceneSetup> sceneSetup = camera->getSceneSetup();
    FCXXCoord backgroundColor(sceneSetup->getBackgroundColor());
    
    // Clear the Canvas
    glClearColor(backgroundColor[0], backgroundColor[1], backgroundColor[2], backgroundColor[3]);
    myglClearDepth(1000.);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    myglEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);
    
    //Scene fog
    myglEnable(GL_FOG);
    
    GLfloat fogColor[] = {1.,0.,0.,1.};//{static_cast<GLfloat>(backgroundColor[0]), static_cast<GLfloat>(backgroundColor[1]), static_cast<GLfloat>(backgroundColor[2]), static_cast<GLfloat>(backgroundColor[3])};
    FCXXCoord cameraTranslation = camera->getTranslation();
    float cameraDistance = cameraTranslation.get3DLength();
    float fogStart = cameraDistance + (camera->getFogFront()  * camera->getSceneSetup()->getScale());
    float fogEnd = fogStart + (camera->getFogDepthRange() * camera->getSceneSetup()->getScale());
    if (fogStart < 0.) fogStart = 0.;
    if (fogEnd < 1.) fogEnd = 1.;
    glFogi(GL_FOG_COORD_SRC, GL_FOG_COORDINATE);
    glFogfv(GL_FOG_COLOR, fogColor);
    glFogi(GL_FOG_MODE, GL_LINEAR);
    glFogf(GL_FOG_START, fogStart);
    glFogf(GL_FOG_END, fogEnd);
    
}

void RendererType::drawTestSquare()
{
    static const GLfloat squareVertices[] = {
        -0.5f, -0.5f, 1,
        0.5f, -0.5f, 1,
        -0.5f,  0.5f, 1,
        0.5f,  0.5f, 1
    };
    
    static const GLubyte squareColors[] = {
        255, 255,   0, 255,
        0,   255, 255, 255,
        0,     0,   0,   0,
        255,   0, 255, 255,
    };
    glClearColor(0.0f, 0.0f, 1.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);
    
    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_COLOR_ARRAY);
    
    glVertexPointer(3, GL_FLOAT, 0, squareVertices);
    glColorPointer(4, GL_UNSIGNED_BYTE, 0, squareColors);
    
    glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
    
    glDisableClientState(GL_VERTEX_ARRAY);
    glDisableClientState(GL_COLOR_ARRAY);
    
}
