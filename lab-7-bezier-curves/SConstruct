import os
env = Environment(ENV = os.environ)
env.Append(CXXFLAGS=["-g", "-std=c++11", "-Wall", "-Werror"])
env.Append(LIBS=["GL","GLU","glut","GLEW"])
# If you are not on the department machines (such as ti-05), you will probably need to comment these out.
env.Append(CPPPATH=["/usr/csshare/pkgs/glew-1.13.0/include"])
env.Append(LIBPATH=["/usr/csshare/pkgs/glew-1.13.0/lib64/"])
env.Append(RPATH=["/usr/csshare/pkgs/glew-1.13.0/lib64/"])
# End optional block
env.Program("bezier",["main.cpp"])
