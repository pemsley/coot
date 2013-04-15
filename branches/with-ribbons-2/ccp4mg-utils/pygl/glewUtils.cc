#ifdef _WIN32
#include <windows.h>
#include "glew.h"
#include "wglew.h"
#include <stdio.h>
#endif

void BlendFuncSeparate(GLenum, GLenum, GLenum, GLenum){
#ifdef _WIN32
  static bool doneglewInit=false;
  if(!doneglewInit){
    int err = glewInit();
    if (GLEW_OK != err) {
      fprintf(stderr, "Error [main]: glewInit failed: %s\r\n", glewGetErrorString(err));
      fflush(stderr);
      return;
    }
  }
  doneglewInit=true;
#endif
  if(GLEW_EXT_blend_func_separate){
    glBlendFuncSeparateEXT(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, GL_ONE, GL_ONE_MINUS_SRC_ALPHA);
  }
}
