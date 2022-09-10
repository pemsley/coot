#!/usr/bin/python
import sys

# top_srcdir is a different directory to top_builddir. We want to make PointLight.hpp
# in top_builddir based on sources (and py file) in top_srcdir.
# argv[1] : "PointLight" (output filename prefix)
# argv[2] : top_srddir (directory where PointLight.{vert,frag} can be found
#
with open(sys.argv[1]+".hpp","w") as hppFile:
    
    hppFile.write('std::string RendererGLSL::'+sys.argv[1]+'FragmentShaderText="')
    with open(sys.argv[2] + '/' + sys.argv[1]+".frag","r") as fragmentFile:
        hppFile.write("#version 120\\n#define highp\\n#define lowp\\n#define mediump\\n"+fragmentFile.read().replace('\n','\\n'))
    hppFile.write('\\n\\0";\n')
    hppFile.write('std::string RendererGLSL::'+sys.argv[1]+'VertexShaderText="')
    with open(sys.argv[2] + '/' + sys.argv[1]+".vert","r") as vertexFile:
        hppFile.write("#version 120\\n#define highp\\n#define lowp\\n#define mediump\\n"+vertexFile.read().replace('\n','\\n'))
    hppFile.write('\\n\\0";\n')


    
