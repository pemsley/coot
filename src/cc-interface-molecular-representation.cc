

#ifdef USE_MOLECULES_TO_TRIANGLES

#ifdef USE_PYTHON
#include <Python.h>  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#endif

#include <cstddef>

// #include "globjects.h" //includes gtk/gtk.h
#include "graphics-info.h"
#include "c-interface.h"
#include "cc-interface.hh"
#include "cc-interface-molecular-representation.hh"

#ifdef USE_PYTHON

// Martin's Triangles

// e.g. 0, "//C", "RampChainsScheme", "Ribbon"
int add_molecular_representation_py(int imol, PyObject *atom_selection_py, PyObject *ColorScheme_py, PyObject *style_py) {

   int status = -1;
   if (is_valid_model_molecule(imol)) {
#ifdef USE_MOLECULES_TO_TRIANGLES
      // check that these are strings
      std::string atom_selection = PyBytes_AS_STRING(PyUnicode_AsUTF8String(atom_selection_py));
      std::string ColorScheme    = PyBytes_AS_STRING(PyUnicode_AsUTF8String(ColorScheme_py));
      std::string style          = PyBytes_AS_STRING(PyUnicode_AsUTF8String(style_py));
      status = graphics_info_t::molecules[imol].add_molecular_representation(atom_selection, ColorScheme, style);
      graphics_draw();
#endif
   }
   return status;

}
#endif // USE_PYTHON

#ifdef USE_GUILE
// e.g. 0, "//C", "RampChainsScheme", "Ribbon"
int add_molecular_representation_scm(int imol, SCM atom_selection_scm, SCM ColorScheme_scm, SCM style_scm) {

   int status = -1;
   if (is_valid_model_molecule(imol)) {
#ifdef USE_MOLECULES_TO_TRIANGLES
      // check that these are strings
      std::string atom_selection = scm_to_locale_string(atom_selection_scm);
      std::string ColorScheme    = scm_to_locale_string(ColorScheme_scm);
      std::string style          = scm_to_locale_string(style_scm);
      std::cout << "Calling add_molecular_representation with " << atom_selection << " " << ColorScheme << "  " << style << std::endl;
      status = graphics_info_t::molecules[imol].add_molecular_representation(atom_selection, ColorScheme, style);
      graphics_draw();
#endif
   }
   return status;
}
#endif // USE_GUILE

void remove_molecular_representation(int imol, int rep_no) {

   if (is_valid_model_molecule(imol)) {
      graphics_info_t::molecules[imol].remove_molecular_representation(rep_no);
      graphics_draw();
   }
}

extern "C" void add_molecular_representation_test() {
   int status = -1;
   std::pair<bool, std::pair<int, coot::atom_spec_t> > active_atom = active_atom_spec();
   if (active_atom.first) {
      int imol = active_atom.second.first;
      std::cout << "Ribbons on molecule " << imol << std::endl;
      if (is_valid_model_molecule(imol)) {
         std::string atom_selection = "//A";
         std::string ColorScheme = "colorRampChainsScheme";
         std::string style = "Ribbon";
         status = graphics_info_t::molecules[imol].add_molecular_representation(atom_selection, ColorScheme, style);
         graphics_info_t::graphics_draw();
      }
   }
}

#else

#endif // USE_MOLECULES_TO_TRIANGLES



void import_bild(const std::string &file_name) {

   class c_info_t {
   public:
      glm::vec3 start_point;
      glm::vec3 end_point;
      coot::colour_holder col;
      float radius;
      c_info_t(const float &x1, const float &y1, const float &z1,
               const float &x2, const float &y2, const float &z2,
               const float &w, const coot::colour_holder &col_in) {
         start_point = glm::vec3(x1, y1, z1);
         end_point   = glm::vec3(x2, y2, z2);
         col = col_in;
         radius = w;
      }
   };

   Shader *shader_p = &graphics_info_t::shader_for_moleculestotriangles;

   auto show_cylinders = [shader_p] (const std::vector<c_info_t> &cv) {
                            meshed_generic_display_object m;
                            for (auto ci : cv) {
                               std::pair<glm::vec3, glm::vec3> pp(ci.start_point, ci.end_point);
                               float h = glm::distance(ci.start_point, ci.end_point);
                               m.add_cylinder(pp, ci.col, ci.radius, 16, true, true,
                                              meshed_generic_display_object::FLAT_CAP,
                                              meshed_generic_display_object::FLAT_CAP);
                            }
                            Material material;
                            gtk_gl_area_attach_buffers(GTK_GL_AREA(graphics_info_t::glareas[0]));
                            m.mesh.setup(shader_p, material);
                            graphics_info_t::generic_display_objects.push_back(m);
                         };

   if (coot::file_exists(file_name)) {
      std::ifstream f(file_name.c_str());
      if (f) {
         std::vector<std::string> lines;
         lines.reserve(4000);
	 std::string line;
	 while (std::getline(f, line)) {
            lines.push_back(line);
	 }
         if (lines.size() > 1) {

            std::vector<c_info_t> cylinder_infos;
            coot::colour_holder current_colour;

            for (auto line : lines) {
               std::vector<std::string> parts = coot::util::split_string_no_blanks(line);
               if (parts.size() == 4) {
                  if (parts[0] == ".color") {
                     try {
                        float r = std::stof(parts[1]);
                        float g = std::stof(parts[2]);
                        float b = std::stof(parts[3]);
                        current_colour = coot::colour_holder(r,g,b);
                     }
                     catch (const std::runtime_error &rte) {
                        std::cout << "WARNING:: failed to read " << rte.what() << std::endl;
                     }
                  }
               }
               if (parts.size() == 8) {
                  if (parts[0] == ".cylinder") {
                     try {
                        float x1 = std::stof(parts[1]);
                        float y1 = std::stof(parts[2]);
                        float z1 = std::stof(parts[3]);
                        float x2 = std::stof(parts[4]);
                        float y2 = std::stof(parts[5]);
                        float z2 = std::stof(parts[6]);
                        float w  = std::stof(parts[7]);
                        cylinder_infos.push_back(c_info_t(x1, y1, z1, x2, y2, z2, w, current_colour));
                     }
                     catch (const std::runtime_error &rte) {
                        std::cout << "WARNING:: failed to read " << rte.what() << std::endl;
                     }
                  }
               }
            }

            if (! cylinder_infos.empty())
               show_cylinders(cylinder_infos);

         } else {
            std::cout << "WARNING:: problematic bild file " << file_name << std::endl;
         }
      }
   } else {
      std::cout << "WARNING:: file not found " << file_name << std::endl;
   }


}
