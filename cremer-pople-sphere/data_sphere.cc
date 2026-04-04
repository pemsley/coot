//
// Created by Jordan Dialpuri on 04/04/2026.
//

#include "data_sphere.hh"
#include <cmath>

float get_energy(float theta, float phi) {
   float base = 0.5f + 0.5f * pow(sin(theta), 2.0f); // low at poles

   float perturb =
         0.10f * sin(2.0f * theta) * cos(3.0f * phi)
       + 0.05f * cos(4.0f * theta)
       + 0.05f * sin(theta) * sin(2.0f * phi);

   return base + perturb;
}
glm::vec4 energy_to_colour(float t) {
   t = glm::clamp(t, 0.0f, 1.0f);

   if (t < 0.5f)
      return glm::mix(glm::vec4(0.1f,0.1f,0.9f,1),
                      glm::vec4(1,1,1,1), t*2.0f);

   return glm::mix(glm::vec4(1,1,1,1),
                   glm::vec4(0.9f,0.1f,0.1f,1),
                   (t-0.5f)*2.0f);
}

coot::simple_mesh_t make_data_sphere(float r, int lat, int lon) {
   coot::simple_mesh_t mesh;
   const float pi = (float)M_PI;

   for (int i = 0; i <= lat; i++) {
      float theta = pi * i / lat;

      for (int j = 0; j <= lon; j++) {
         float phi = 2*pi * j / lon;

         glm::vec3 n(sin(theta)*cos(phi),
                     sin(theta)*sin(phi),
                     cos(theta));

         mesh.vertices.emplace_back(
             r*n,
             n,
             energy_to_colour(get_energy(theta, phi))
         );
      }
   }

   for (int i = 0; i < lat; i++)
      for (int j = 0; j < lon; j++) {
         unsigned a=i*(lon+1)+j;
         unsigned b=(i+1)*(lon+1)+j;
         unsigned c=(i+1)*(lon+1)+j+1;
         unsigned d=i*(lon+1)+j+1;

         mesh.triangles.emplace_back(a,b,c);
         mesh.triangles.emplace_back(a,c,d);
      }

   return mesh;
}