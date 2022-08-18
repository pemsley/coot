
#include "RendererGLSL.hpp"
// #include "PointLight.hpp" // now we read the shaders from files, not the header
// #include "GL/glew.h"

// init statics that use to be in PointLight.hpp
std::string RendererGLSL::PointLightFragmentShaderText;
std::string RendererGLSL::PointLightVertexShaderText;


#ifdef HAVE_CONFIG_H
#include "config.h"
#else
#define HAVE_OPENGL_GL_H
#endif

#if 0 // -------------------------------------------
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
#endif // 0 -------------------------------------

#include <epoxy/gl.h>


#include <stdio.h>
#include <iostream>
#include <fstream>
#include "DisplayPrimitive.h"
#include "VertexColorNormalPrimitive.h"


void RendererGLSL::setProgram(int _program)
{
#ifdef DEBUG_MINE
    std::cout << "In setProgram";
#endif
    program=_program;
};

void RendererGLSL::vboRenderVCN(VertexColorNormalPrimitive *prim)
{
#ifdef DEBUG_MINE
    std::cout << "In GLSL::vboRenderVCN render " << prim->nVertices();
#endif
    RendererHandles rendererHandles;
    std::map<DisplayPrimitive *, RendererHandles>::iterator mapIter = allocatedHandles.find(prim);
    if (mapIter != allocatedHandles.end()){
        rendererHandles = mapIter->second;
    }
    else {
#ifdef DEBUG_MINE
        std::cout << "Allocating new handles for "<<prim << std::endl;
#endif
        GLuint vertexHandle, indexHandle, arrayObjectHandle;

        glGenVertexArrays(1, &arrayObjectHandle);
        glBindVertexArray(arrayObjectHandle);
        rendererHandles.arrayObjectHandle = arrayObjectHandle;

        glGenBuffers(1, &vertexHandle);
        rendererHandles.vertexHandle = vertexHandle;

        glBindBuffer(GL_ARRAY_BUFFER, vertexHandle);
        glBufferData(GL_ARRAY_BUFFER, prim->nVertices()*sizeof(VertexColorNormalPrimitive::VertexColorNormal),
                     prim->getVertexColorNormalArray(), GL_DYNAMIC_DRAW);

        myglEnableClientState(GL_VERTEX_ARRAY);
        myglVertexPointer(3, GL_FLOAT, sizeof(VertexColorNormalPrimitive::VertexColorNormal),
                          (void*)offsetof(VertexColorNormalPrimitive::VertexColorNormal,vertex));

        myglEnableClientState(GL_COLOR_ARRAY);
        myglColorPointer(4, GL_FLOAT, sizeof(VertexColorNormalPrimitive::VertexColorNormal),
                         (void*)offsetof(VertexColorNormalPrimitive::VertexColorNormal,color));

        myglEnableClientState(GL_NORMAL_ARRAY);
        myglNormalPointer(GL_FLOAT, sizeof(VertexColorNormalPrimitive::VertexColorNormal),
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
        
        glDisableVertexAttribArray(0);
        glDisableVertexAttribArray(1);
        glDisableVertexAttribArray(2);
    }
    
    GLenum err = glGetError(); std::cout << "   Here in vboRenderVCN() pre  glBindVertexArray() err " << err << std::endl;
    glBindVertexArray(rendererHandles.arrayObjectHandle);
    err = glGetError(); std::cout << "   Here in vboRenderVCN() post glBindVertexArray() err " << err << std::endl;
    glEnableVertexAttribArray(0);
    glEnableVertexAttribArray(1);
    glEnableVertexAttribArray(2);
    int nTriangleVertices = 3*prim->nTriangles();
    err = glGetError(); std::cout << "   Here in vboRenderVCN() pre glDrawElements() err " << err << std::endl;
    err = glGetError(); std::cout << "   Here in vboRenderVCN() pre glDrawElements() err " << err << std::endl;
    glDrawElements(GL_TRIANGLES, nTriangleVertices, kGLIndexType, (void*)0);
    err = glGetError(); std::cout << "   Here in vboRenderVCN() post glDrawElements() err " << err << std::endl;

 
    glBindVertexArray(0);
}

void RendererGLSL::renderVertexColorNormalPrimitive(VertexColorNormalPrimitive *prim)
{
    GLenum err = glGetError(); std::cout << "      Here in RendererGLSL::renderVertexColorNormalPrimitive() --start-- about to glUseProgram() "
                                         << program << " err " << err << std::endl;
    glUseProgram(program);
    err = glGetError(); std::cout << "      Here A in RendererGLSL::renderVertexColorNormalPrimitive() err " << err << std::endl;
    std::cout << "In RendererGLSL::renderVertexColorNormalPrimitive() A\n";

    err = glGetError(); std::cout << "      Here B in RendererGLSL::renderVertexColorNormalPrimitive() calling vboRenderVCN err " << err << std::endl;
    vboRenderVCN(prim);
    err = glGetError(); std::cout << "      Here C in RendererGLSL::renderVertexColorNormalPrimitive() calling myglDisable() err " << err << std::endl;
    // myglDisable(GL_COLOR_MATERIAL);
    err = glGetError(); std::cout << "      Here D in RendererGLSL::renderVertexColorNormalPrimitive() calling glUseProgram(0) err " << err << std::endl;
    glUseProgram(0);
    err = glGetError(); std::cout << "      Here E in RendererGLSL::renderVertexColorNormalPrimitive() --end-- err " << err << std::endl;
}

void RendererGLSL::init()
{
    std::cout << "RendererGLSL::init(): Off to load shaders...\n\n";
    loadShaders();
    GLenum err = glGetError(); std::cout << "\n      Here in RendererGLSL::init() after loadShaders() err " << err << std::endl;
    int n;

    // auto a = glGetString(GL_VERSION);
    // std::cout << "OpenGL version: " << a << std::endl;

    glGetIntegerv(GL_MAX_VERTEX_UNIFORM_COMPONENTS, &n);
    err = glGetError(); std::cout << "      Here in RendererGLSL::init() B err " << err << std::endl;
    std::cout << "INFO:: GL_MAX_VERTEX_UNIFORM_COMPONENTS " << n << std::endl;
    glUseProgram(program);
    err = glGetError(); std::cout << "      Here in RendererGLSL::init() C0 err " << err << std::endl;

    const char *names[] = {
        "mygl_ModelViewMatrix",
        "mygl_ProjectionMatrix",
        "mygl_NormalMatrix",
        "mygl_UseLight0",
        "mygl_LightSource[0].ambient",
        "mygl_LightSource[0].position",
        "mygl_LightSource[0].diffuse",
        "mygl_LightSource[0].specular",
        "mygl_LightSource[0].constantAttenuation",
        "mygl_LightSource[0].linearAttenuation",
        "mygl_LightSource[0].quadraticAttenuation",
        "mygl_UseLight1",
        "mygl_LightSource[1].position",
        "mygl_LightSource[1].ambient",
        "mygl_LightSource[1].diffuse",
        "mygl_LightSource[1].specular",
        "mygl_LightSource[1].constantAttenuation",
        "mygl_LightSource[1].linearAttenuation",
        "mygl_LightSource[1].quadraticAttenuation",
        "mygl_UseColorArray",
        "mygl_FrontMaterial.emission",
        "mygl_FrontMaterial.ambient",
        "mygl_FrontMaterial.diffuse",
        "mygl_FrontMaterial.shininess",
        "mygl_FrontMaterial.specular",
        "mygl_BackMaterial.emission",
        "mygl_BackMaterial.ambient",
        "mygl_BackMaterial.diffuse",
        "mygl_BackMaterial.shininess",
        "mygl_BackMaterial.specular",
        "mygl_FrontLightModelProduct.sceneColor"
    };
    for (int i=0; i<(sizeof(names) / sizeof(char *)); i++){
        int location = glGetUniformLocation(program, names[i]);
        std::cout << "location of "  << names[i] << " is " << location << std::endl;
        if (location == -1)
          std::cout << "      Some problem with glGetUniformLocation() for " << names[i] << std::endl;
        err = glGetError(); std::cout << "      Here in RendererGLSL::init() loop D1 i " << i << " err "
                                      << err << std::endl;
        uniforms[std::string(names[i])] = location;
        std::cout << "Loc of "  << names[i] << " is " << location << std::endl;
    }
    err = glGetError(); std::cout << "      Here in RendererGLSL::init() D2 err " << err << std::endl;
    GLfloat nullColor[] = {0.,0.,0.,1.};
    glUniform4fv(uniforms["mygl_FrontMaterial.emission"], 4, nullColor);
    err = glGetError(); std::cout << "      Here in RendererGLSL::init() E1a err " << err << std::endl;
    err = glGetError(); std::cout << "      Here in RendererGLSL::init() E1b err " << err << std::endl;
    glUniform4fv(uniforms["mygl_FrontMaterial.ambient"], 4, nullColor);
    err = glGetError(); std::cout << "      Here in RendererGLSL::init() E2a err " << err << std::endl;
    err = glGetError(); std::cout << "      Here in RendererGLSL::init() E2b err " << err << std::endl;
    glUniform4fv(uniforms["mygl_FrontMaterial.specular"], 4, nullColor);
    err = glGetError(); std::cout << "      Here in RendererGLSL::init() E3a err " << err << std::endl;
    err = glGetError(); std::cout << "      Here in RendererGLSL::init() E3b err " << err << std::endl;
    glUniform4fv(uniforms["mygl_FrontLightModelProduct.sceneColor"], 4, nullColor);
    glUniform1i(uniforms["mygl_UseLight0"], 1);
    glUniform1i(uniforms["mygl_UseLight1"], 1);
    glUniform1i(uniforms["mygl_UseColorArray"], 1);
    err = glGetError(); std::cout << "      Here in RendererGLSL::init() Ef err " << err << std::endl;
    glUseProgram(0);
    err = glGetError(); std::cout << "      Here in RendererGLSL::init() F err " << err << std::endl;
}


void RendererGLSL::myglMaterialfv(int faceEnum, int propertyEnum, float *values){
    std::string frontProperty("mygl_FrontMaterial.");
    std::string backProperty("mygl_BackMaterial.");
    switch (propertyEnum){
        case GL_AMBIENT:
            frontProperty.append("ambient");
            backProperty.append("ambient");
            break;
        case GL_DIFFUSE:
            frontProperty.append("diffuse");
            backProperty.append("diffuse");
            break;
        case GL_SPECULAR:
            frontProperty.append("specular");
            backProperty.append("specular");
            break;
        case GL_EMISSION:
            frontProperty.append("emission");
            backProperty.append("emission");
            break;
    }
    switch (faceEnum){
        case GL_FRONT:
            glUniform4fv(uniforms[frontProperty], 4, values);
            //std::cout << uniforms[frontProperty] << "("<<frontProperty<<") set to  " << values[0] << "  " << values[1] << "  " << values[2] << "  " << values[3] << "  " << std::endl;
            break;
        case GL_BACK:
            glUniform4fv(uniforms[backProperty], 4, values);
            //std::cout << uniforms[backProperty] << "("<<backProperty<<") set to  " << values[0] << std::endl;
            break;
        case GL_FRONT_AND_BACK:
            glUniform4fv(uniforms[frontProperty], 4, values);
            //std::cout << uniforms[frontProperty] << "("<<frontProperty<<") set to  " << values[0] << "  " << values[1] << "  " << values[2] << "  " << values[3] << "  " << std::endl;
            glUniform4fv(uniforms[backProperty], 4, values);
            //std::cout << uniforms[backProperty] << "("<<backProperty<<") set to  " << values[0] << std::endl;
            break;
        default:
            break;
    }
}

void RendererGLSL::myglMaterialf(int faceEnum, int propertyEnum, float value){
    std::string frontProperty("mygl_FrontMaterial.");
    std::string backProperty("mygl_BackMaterial.");
    switch (propertyEnum){
        case GL_SHININESS:
            frontProperty.append("shininess");
            backProperty.append("shininess");
            break;
        default:
            break;
    }
    switch (faceEnum){
        case GL_FRONT:
            glUniform1f(uniforms[frontProperty], value);
            break;
        case GL_BACK:
            glUniform1f(uniforms[backProperty], value);
            break;
        case GL_FRONT_AND_BACK:
            glUniform1f(uniforms[frontProperty], value);
            glUniform1f(uniforms[backProperty], value);
            break;
        default:
            break;
    }
}

void RendererGLSL::myglEnableClientState(int clientStateEnum){
    if (clientStateEnum == GL_VERTEX_ARRAY){
        glEnableVertexAttribArray(0);
    }
    else if (clientStateEnum == GL_COLOR_ARRAY){
        glUniform1i(uniforms["mygl_UseColorArray"], 1);
        glEnableVertexAttribArray(1);
    }
    else if (clientStateEnum == GL_NORMAL_ARRAY){
        glEnableVertexAttribArray(2);
    }
}

void RendererGLSL::myglVertexPointer(int size, int type, int stride, const void *pointer){
    glVertexAttribPointer(0, size, type, 1, stride, pointer);
}

void RendererGLSL::myglColorPointer(int size, int type, int stride, const void *pointer){
    glVertexAttribPointer(1, size, type, 1, stride, pointer);
}

void RendererGLSL::myglNormalPointer(int type, int stride, const void *pointer){
    glVertexAttribPointer(2, 3, type, 1, stride, pointer);
}

void RendererGLSL::loadShaders()
{
    std::cout << "In loadShaders...\n";
    GLuint vs = glCreateShader(GL_VERTEX_SHADER);
    //loadShaderFile("/Users/martin/Dropbox/Programming/MoleculesToTriangles/CXXClasses/PointLight.vert", vs);

    std::pair<bool, std::string> shader_code = readShaderFile("PointLight.vert.glsl");
    if (shader_code.first)
       PointLightVertexShaderText = shader_code.second;

    const char *vertexSource = PointLightVertexShaderText.c_str();
    GLint vertexShaderSize = PointLightVertexShaderText.size();
    glShaderSource(vs, 1, &vertexSource, NULL);
     
    glCompileShader(vs);
    int status;
    glGetShaderiv(vs, GL_COMPILE_STATUS, &status);
    if(status == GL_FALSE)
    {
        char infoLog[1024];
        glGetShaderInfoLog(vs, 1024, NULL, infoLog);
        std::cout << "The vertex shader failed to compile with the following errors:" << std::endl
                  << infoLog << std::endl;
        std::cout << PointLightVertexShaderText << std::endl;
        glDeleteShader(vs);
    }
    else std::cout << "The vertex shader compiled without errors\n";
    
    GLuint fs = glCreateShader(GL_FRAGMENT_SHADER);
    //loadShaderFile("/Users/martin/Dropbox/Programming/MoleculesToTriangles/CXXClasses/PointLight.frag", fs);

    shader_code = readShaderFile("PointLight.frag.glsl");
    if (shader_code.first)
       PointLightFragmentShaderText = shader_code.second;

    const char *fragmentSource = PointLightFragmentShaderText.c_str();
    GLint fragmentShaderSize = PointLightFragmentShaderText.size();
    glShaderSource(fs, 1, &fragmentSource, &fragmentShaderSize);
 
    glCompileShader(fs);
    glGetShaderiv(fs, GL_COMPILE_STATUS, &status);
    if(status == GL_FALSE)
    {
        char infoLog[1024];
        glGetShaderInfoLog(fs, 1024, NULL, infoLog);
        std::cout << "The fragment shader failed to compile with the following errors:" << std::endl
                  << infoLog << std::endl;
        std::cout << PointLightFragmentShaderText << std::endl;
        glDeleteShader(fs);
    }
    else std::cout << "The fragment shader compiled without errors\n";

    program = glCreateProgram();
    glAttachShader(program, fs);
    glAttachShader(program, vs);
    
    glBindAttribLocation(program, 0, "mygl_Vertex");
    glBindAttribLocation(program, 1, "mygl_Color");
    glBindAttribLocation(program, 2, "mygl_Normal");
    
    glLinkProgram(program);
    GLint infologLength = 0;
    glGetProgramiv(program, GL_INFO_LOG_LENGTH, &infologLength);
    std::cerr<<"Link Log Length "<<infologLength<<"\n";
    if(infologLength > 1)
    {
        char *infoLog = new char[infologLength];
        GLint charsWritten  = 0;
        
        glGetProgramInfoLog(program, infologLength, &charsWritten, infoLog);
        
        std::cerr<<infoLog<<std::endl;
        delete [] infoLog;
        glGetProgramiv(program, GL_LINK_STATUS, &infologLength);
        if(infologLength == GL_FALSE)
        {
            std::cout<<"Program link failed exiting \n";
        }
    }
}

std::shared_ptr<Renderer> RendererGLSL::create()
{
    auto result = std::shared_ptr<Renderer>(new RendererGLSL());
    return result;
}

std::pair<bool, std::string> RendererGLSL::readShaderFile(std::string strFilename) const {

    std::string s;
	 bool status = false;
    std::ifstream shaderSource(strFilename);
	 if (shaderSource.is_open()) {
       s = std::string((std::istreambuf_iterator<char>(shaderSource)), std::istreambuf_iterator<char>());
		 shaderSource.close();
		 status = true;
	 } else {
        std::cerr << " File not found " << strFilename << std::endl;
	 }
	 return std::pair<bool, std::string>(status, s);
}

bool RendererGLSL::loadShaderFile(std::string strFilename, GLuint iHandle)
{
    std::ifstream shaderSource(strFilename);
    if (!shaderSource.is_open())
    {
        std::cerr<< " File not found "<< strFilename.c_str()<< std::endl;
        return false;
    }
    // now read in the data
    std::string strSource = std::string((std::istreambuf_iterator<char>(shaderSource)), std::istreambuf_iterator<char>());
    shaderSource.close();
    strSource+="\0";
    //pass the code to OGL
    const char* data=strSource.c_str();
    glShaderSource(iHandle, 1, &data, NULL);
    return true;
}
