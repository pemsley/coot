//
// Created by Jordan Dialpuri on 04/04/2026.
//

#include "data_sphere.hh"
#include <algorithm>
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

glm::vec4 density_to_colour(float t) {
   t = glm::clamp(t, 0.0f, 1.0f);

   // 5-stop heat map: dark-navy → blue → teal → yellow → red
   const glm::vec4 stops[5] = {
      { 0.05f, 0.05f, 0.20f, 1.0f },  // 0.00 – no data  (dark navy)
      { 0.10f, 0.20f, 0.90f, 1.0f },  // 0.25 – sparse   (blue)
      { 0.10f, 0.80f, 0.70f, 1.0f },  // 0.50 – moderate (teal)
      { 0.90f, 0.85f, 0.10f, 1.0f },  // 0.75 – dense    (yellow)
      { 0.90f, 0.10f, 0.10f, 1.0f },  // 1.00 – hotspot  (red)
   };

   const float scaled = t * 4.0f;
   const int   idx    = std::min(static_cast<int>(scaled), 3);
   const float frac   = scaled - static_cast<float>(idx);
   return glm::mix(stops[idx], stops[idx + 1], frac);
}

coot::simple_mesh_t make_data_sphere(float r, int lat, int lon,
                                     const std::vector<float>& density_grid) {
   coot::simple_mesh_t mesh;
   const float pi = static_cast<float>(M_PI);

   for (int i = 0; i <= lat; i++) {
      const float theta = pi * static_cast<float>(i) / static_cast<float>(lat);

      for (int j = 0; j <= lon; j++) {
         const float phi = 2.0f * pi * static_cast<float>(j) / static_cast<float>(lon);

         const glm::vec3 n(sin(theta) * cos(phi),
                           sin(theta) * sin(phi),
                           cos(theta));

         // Nearest-neighbour lookup into the lat x lon density grid
         const int bi = std::min(i, lat - 1);
         const int bj = j % lon;
         const float d = density_grid[bi * lon + bj];

         mesh.vertices.emplace_back(r * n, n, density_to_colour(d));
      }
   }

   for (int i = 0; i < lat; i++)
      for (int j = 0; j < lon; j++) {
         unsigned a =  i      * (lon + 1) + j;
         unsigned b = (i + 1) * (lon + 1) + j;
         unsigned c = (i + 1) * (lon + 1) + j + 1;
         unsigned d =  i      * (lon + 1) + j + 1;

         mesh.triangles.emplace_back(a, b, c);
         mesh.triangles.emplace_back(a, c, d);
      }

   return mesh;
}