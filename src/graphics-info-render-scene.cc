/*
 * src/graphics-info-render-scene.cc
 *
 * Copyright 2022 by Medical Research Council
 * Author: Paul Emsley
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copies of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 * See http://www.gnu.org/licenses/
 *
 */

#include "analysis/stats.hh"
#define GLM_ENABLE_EXPERIMENTAL // # for norm things
// #include <glm/ext.hpp>
#include <glm/gtx/string_cast.hpp>  // to_string()

#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <gtk/gtk.h>
#include <epoxy/gl.h>


#include "graphics-info.h"

// renderCube() renders a 1x1 3D cube in NDC.
// -------------------------------------------------

unsigned int quadVAO = 0;
unsigned int quadVBO;

// static
void
graphics_info_t::render_scene_sans_depth_blur(stereo_eye_t eye,
                                              Shader *shader_for_tmeshes_p, Shader *shader_for_meshes_p,
                                              Shader *shader_for_tmeshes_with_shadows_p,
                                              Shader *shader_for_meshes_with_shadows_p,
                                              int width, int height) {

   auto renderQuad = [] () {
                        if (quadVAO == 0) {
                           float quadVertices[] = {
                                                   // positions        // texture Coords
                                                   -1.0f,  1.0f, 0.0f, 0.0f, 1.0f,
                                                   -1.0f, -1.0f, 0.0f, 0.0f, 0.0f,
                                                    1.0f,  1.0f, 0.0f, 1.0f, 1.0f,
                                                    1.0f, -1.0f, 0.0f, 1.0f, 0.0f,
                           };
                           // setup plane VAO
                           GLenum err = glGetError();
                           if (err)
                              std::cout << "GL ERROR:: lambda renderQuad() quadVAO == 0 --- start --- " << err << std::endl;
                           glGenVertexArrays(1, &quadVAO);
                           glGenBuffers(1, &quadVBO);
                           glBindVertexArray(quadVAO);
                           glBindBuffer(GL_ARRAY_BUFFER, quadVBO);
                           glBufferData(GL_ARRAY_BUFFER, sizeof(quadVertices), &quadVertices, GL_STATIC_DRAW);
                           glEnableVertexAttribArray(0);
                           glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void*)0);
                           glEnableVertexAttribArray(1);
                           glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void*)(3 * sizeof(float)));
                           err = glGetError();
                           if (err)
                              std::cout << "GL ERROR:: render_scene_sans_depth_blur() renderQuad() quadVAO == 0 --- end setup --- " << err << std::endl;
                        }
                        GLenum err = glGetError();
                        if (err)
                           std::cout << "GL ERROR:: render_scene_sans_depth_blur() lambda renderQuad() A " << err << std::endl;
                        glBindVertexArray(quadVAO);
                        err = glGetError();
                        if (err)
                           std::cout << "GL ERROR:: render_scene_sans_depth_blur() lambda renderQuad() B " << err << std::endl;
                        glEnableVertexAttribArray(0);
                        err = glGetError();
                        if (err)
                           std::cout << "GL ERROR:: render_scene_sans_depth_blur() lambda renderQuad() C1 " << err << std::endl;
                        glEnableVertexAttribArray(1);
                        err = glGetError();
                        if (err)
                           std::cout << "GL ERROR:: render_scene_sans_depth_blur() lambda renderQuad() C2 " << err << std::endl;
                        glBindVertexArray(quadVAO);
                        err = glGetError();
                        if (err)
                           std::cout << "GL ERROR:: render_scene_sans_depth_blur() lambda renderQuad() D " << err << std::endl;
                        glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
                        err = glGetError();
                        if (err)
                           std::cout << "GL ERROR:: render_scene_sans_depth_blur() lambda renderQuad() E " << err << std::endl;
                        glBindVertexArray(0);
                        err = glGetError();
                        if (err)
                           std::cout << "GL ERROR:: render_scene_sans_depth_blur() lambda renderQuad() F " << err << std::endl;
                     };

   auto render_to_shadow_map = [] () {

      // std::cout << "---------------------------------------- render_to_shadow_map() ------------------------------- "
      // << shadow_strength << std::endl;

                                 if (shadow_strength == 0.0f) return;
                                 GLenum err = glGetError();
                                 if (err)
                                    std::cout << "GL ERROR:: lambda render_to_shadow_map() --- start --- " << err << std::endl;

                                  graphics_info_t g;
                                  glViewport(0, 0, shadow_texture_width, shadow_texture_height);
                                  glBindFramebuffer(GL_FRAMEBUFFER, shadow_depthMap_framebuffer);

                                  glClear(GL_DEPTH_BUFFER_BIT);
                                  unsigned int light_index = 0;
                                  // make these static at some stage
                                  g.draw_Models_for_shadow_map(light_index);
                                  g.draw_molecules_for_shadow_map(light_index); // for coot non-Model objects

                               };

   GLenum err = glGetError();
   if (err)
      std::cout << "GL ERROR:: render_scene_sans_depth_blur() --- start --- " << err << std::endl;

   graphics_info_t di;
   GtkWidget *gl_area = graphics_info_t::glareas[0];

   int ssao_blur_size = graphics_info_t::ssao_blur_size;

   // std::cout << "render_scene_sans_depth_blur() " << std::endl;

   if (di.displayed_image_type == SHOW_AO_SCENE) {

      {

         int w = width;
         int h = height;

         // render
         // ------

         // 1. geometry pass: render scene's geometry/color data into gbuffer
         // -----------------------------------------------------------------
         // glBindFramebuffer(GL_FRAMEBUFFER, di.gBufferFBO);

         di.framebuffer_for_ssao_gbuffer.bind();
         glClearColor(0.0, 0.0, 0.0, 1.0);
         glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

         bool do_orthographic_projection = ! perspective_projection_flag;
         glm::mat4 projection_matrix = di.get_projection_matrix(do_orthographic_projection, w, h);
         glm::mat4 view_matrix       = di.get_view_matrix();
         glm::mat4 model_matrix      = di.get_model_matrix();
         model_matrix = glm::translate(model_matrix, glm::vec3(1.1f, 4.0f, 2.2f));
         model_matrix = glm::scale(model_matrix, glm::vec3(0.5f, 0.5f, 0.5f));

         // std::cout << "We are using shaderGeometryPass " << shaderGeometryPass.name << std::endl;
         di.shaderGeometryPass.Use();
         di.shaderGeometryPass.set_mat4_for_uniform("model", model_matrix);
         di.shaderGeometryPass.set_mat4_for_uniform("view", view_matrix);
         di.shaderGeometryPass.set_mat4_for_uniform("projection", projection_matrix);

         if (false) {
            std::cout << "sending model      " << glm::to_string(model_matrix) << std::endl;
            std::cout << "sending view       " << glm::to_string(view_matrix) << std::endl;
            std::cout << "sending projection " << glm::to_string(projection_matrix) << std::endl;
         }

         shaderGeometryPass.set_int_for_uniform("invertedNormals", 0);

         model_matrix = di.get_model_matrix();
         model_matrix = glm::translate(model_matrix, glm::vec3(0.0f, 0.5f, 0.0));
         model_matrix = glm::rotate(model_matrix, glm::radians(-90.0f), glm::vec3(1.0, 0.0, 0.0));
         model_matrix = glm::scale(model_matrix, glm::vec3(1.0f));
         di.shaderGeometryPass.set_mat4_for_uniform("model", model_matrix);

         // render_3d_scene_for_ssao(); // 20220130-PE render_3d_scene() was coot before Models/crows
         draw_models_for_ssao();
         draw_molecules_for_ssao(eye);

         glBindFramebuffer(GL_FRAMEBUFFER, 0);


         // 2. generate SSAO texture
         // ------------------------
         glBindFramebuffer(GL_FRAMEBUFFER, di.ssaoFBO);
         // std::cout << "Here in render_scene_sans_depth_blur() bound ssaFBO frame buffer " << di.ssaoFBO << std::endl;

         glClear(GL_COLOR_BUFFER_BIT);
         // std::cout << "We are using shaderSSAO " << shaderSSAO.name << std::endl;
         di.shaderSSAO.Use();
         di.shaderSSAO.set_int_for_uniform("gPosition", 0);
         di.shaderSSAO.set_int_for_uniform("gNormal",   1);
         di.shaderSSAO.set_int_for_uniform("texNoise",  2);

         // Send kernel + rotation
         di.shaderSSAO.set_int_for_uniform("n_ssao_kernel_samples", di.n_ssao_kernel_samples);
         for (unsigned int i = 0; i < di.n_ssao_kernel_samples; ++i) {
            // std::cout << "kernel sample " << i << " " << glm::to_string(di.ssaoKernel[i]) << std::endl;
            di.shaderSSAO.set_vec3_for_uniform("samples[" + std::to_string(i) + "]", di.ssaoKernel[i]);
         }
         di.shaderSSAO.set_mat4_for_uniform("projection", projection_matrix);

         if (false)
            std::cout << "sending ssao radius " << SSAO_radius << " bias " << SSAO_bias
                      << " n-kernel-samples "  << n_ssao_kernel_samples << std::endl;
         shaderSSAO.set_float_for_uniform("radius", SSAO_radius);
         shaderSSAO.set_float_for_uniform("bias",   SSAO_bias);
         glActiveTexture(GL_TEXTURE0);
         glBindTexture(GL_TEXTURE_2D, di.framebuffer_for_ssao_gbuffer.gPosition);
         glActiveTexture(GL_TEXTURE1);
         glBindTexture(GL_TEXTURE_2D, di.framebuffer_for_ssao_gbuffer.gNormal);
         glActiveTexture(GL_TEXTURE2);
         glBindTexture(GL_TEXTURE_2D, di.noiseTexture);
         if (false)
            std::cout << "debug:: in render_scene() gPosition is " << di.framebuffer_for_ssao_gbuffer.gPosition
                      << " and gNormal is " << di.framebuffer_for_ssao_gbuffer.gNormal << std::endl;
         renderQuad();
         glBindFramebuffer(GL_FRAMEBUFFER, 0);

         // std::cout << "debug gPosition gNormal noiseTexture " << di.gPosition << " "
         //           << di.gNormal << " " << di.noiseTexture << std::endl;

         err = glGetError();
         if (err)
            std::cout << "GL ERROR:: render_scene_sans_depth_blur() post noisetexture " << err << std::endl;

         // 3. blur SSAO texture to remove noise
         // ------------------------------------
         glBindFramebuffer(GL_FRAMEBUFFER, di.ssaoBlurFBO);
         // std::cout << "Here in render_scene_sans_depth_blur() bound ssaFBOblur frame buffer " << di.ssaoBlurFBO
         // << std::endl;
         glClear(GL_COLOR_BUFFER_BIT);
         di.shaderSSAOBlur.Use();
         di.shaderSSAOBlur.set_int_for_uniform("ssaoInput", 0);
         di.shaderSSAOBlur.set_int_for_uniform("blur_size", ssao_blur_size);
         glActiveTexture(GL_TEXTURE0);
         glBindTexture(GL_TEXTURE_2D, di.ssaoColorBuffer);
         err = glGetError();
         if (err)
            std::cout << "GL ERROR:: render_scene_sans_depth_blur() post bind ssaoColorBuffer " << err << std::endl;
         renderQuad();
         err = glGetError();
         if (err)
            std::cout << "GL ERROR:: render_scene_sans_depth_blur() post SSAO renderQuad() " << err << std::endl;
         glBindFramebuffer(GL_FRAMEBUFFER, 0);
         // std::cout << "render_scene_sans_depth_blur() here 4 " << std::endl;

         err = glGetError();
         if (err)
            std::cout << "GL ERROR:: render_scene_sans_depth_blur() post ssaoColorBuffer " << err << std::endl;

         bool show_shadow_map = false;

         if (show_shadow_map) { // debugging shadow map

            render_to_shadow_map();

            // now we have a texture in a framebuffer. I want to see what that looks like
            di.attach_buffers(); // switch back to the defaut framebuffer
            glViewport(0, 0, graphics_x_size, graphics_y_size);
            float gamma_pe = 1.2f;
            glClearColor(pow(0.07f, gamma_pe), pow(0.13f, gamma_pe), pow(0.17f, gamma_pe), 1.0f);
            // std::cout << "clearing with gamma colour" << std::endl;
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

            glActiveTexture(GL_TEXTURE0 + 0);
            glBindTexture(GL_TEXTURE_2D, shadow_depthMap_texture);
            di.tmesh_for_shadow_map.draw(&di.shader_for_shadow_map_image_texture_mesh);

         } else {

            render_to_shadow_map(); // needed? yes.

            {

               // std::cout << "   framebuffer_scale: " << framebuffer_scale << std::endl;
               // glViewport(0,0, width * framebuffer_scale, height * framebuffer_scale);
               glViewport(0,0, width, height);
               di.framebuffer_for_effects.bind();

               // are these needed if the background image is drawn?
               glClearColor(background_colour.r, background_colour.g, background_colour.b, 1.0);
               glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

               if (draw_background_image_flag) {
                  // 20220929-PE meh. I don't understand where this draw call should go in "Fancy" mode.

                  texture_for_background_image.Bind(0);
                  tmesh_for_background_image.draw(&shader_for_background_image, HUDTextureMesh::TOP_LEFT);
               }

               glEnable(GL_DEPTH_TEST);
               glDepthFunc(GL_LESS);
               glDisable(GL_BLEND);
               bool draw_shadows = true;
               // std::cout << "show_just_shadows " << show_just_shadows << std::endl;

               // maybe doing the shadows separately is the right approach, in fact.
               // because currently *every* draw call of each of the tmeshes
               // samples the shadow map (for every pixel). And the shadow map
               // samples a 3x3 patch for each pixel.
               // Maybe it's better to draw (like SHOW_BASIC_SCENE) to a framebuffer
               // colour texture without shadows and then use a new shader to combine
               // that image with the shadow map.
               //
               // I will write a shader for the meshes (shader_for_meshes_with_shadows_p)
               // that doesn't draw shadows.

               // std::cout << "render_scene_sans_depth_blur() here 5 " << std::endl;

               // std::cout << "draw_models_with_shadows() with shader " << shader_for_meshes_with_shadows_p->name << std::endl;
               di.draw_models_with_shadows(eye,
                                           shader_for_tmeshes_with_shadows_p,
                                           shader_for_meshes_with_shadows_p,
                                           width, height,
                                           draw_shadows,
                                           shadow_strength,
                                           show_just_shadows);

               render_3d_scene_with_shadows(eye); // no longer does a glClear()

            }


            // if I enable di.attach_buffers() above and comment out below, then it renders

            di.attach_buffers();
            glViewport(0, 0, graphics_x_size, graphics_y_size);
            // useful for debugging:
            // GLint local_fbo;
            // glGetIntegerv(GL_FRAMEBUFFER_BINDING, &local_fbo);

            {

               glClearColor(background_colour.r, background_colour.g, background_colour.b, 1.0);
               glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // needed?

               glBindVertexArray(di.screen_AO_quad_vertex_array_id); // maybe make this a HUDTexture?
               glm::vec4 bg_col(background_colour, 1.0f);

               glActiveTexture(GL_TEXTURE0 + 0);
               glBindTexture(GL_TEXTURE_2D, di.framebuffer_for_effects.get_texture_colour());
               glActiveTexture(GL_TEXTURE0 + 1);
               glBindTexture(GL_TEXTURE_2D, di.framebuffer_for_effects.get_texture_depth());
               glActiveTexture(GL_TEXTURE0 + 2);
               glBindTexture(GL_TEXTURE_2D, di.ssaoColorBufferBlur);
               di.shader_for_effects.Use();
               // std::cout << "using shader " << shader_for_effects.name << std::endl;
               shader_for_effects.set_int_for_uniform("screenTexture", 0);
               shader_for_effects.set_int_for_uniform("screenDepth",   1);
               shader_for_effects.set_int_for_uniform("ssao",          2); // sampler2D
               if (graphics_info_t::graphics_is_gl_es)
                  std::cout << "INFO:: no effects shader in OpenGL-ES mode" << std::endl;
               else
                  shader_for_effects.set_int_for_uniform("effects_output_type", effects_shader_output_type);
               shader_for_effects.set_bool_for_uniform("use_ssao", di.use_ssao);
               shader_for_effects.set_bool_for_uniform("show_ssao", di.show_just_ssao);
               shader_for_effects.set_float_for_uniform("ssao_strength", di.ssao_strength);
               shader_for_effects.set_float_for_uniform("brightness", effects_brightness);
               shader_for_effects.set_float_for_uniform("gamma",      effects_gamma);

               err = glGetError();
               if (err)
                  std::cout << "GL ERROR:: render_scene_sans_depth_blur() D err " << err << std::endl;

               glEnableVertexAttribArray(0);
               glEnableVertexAttribArray(1);
               glDrawArrays(GL_TRIANGLES, 0, 6);
               err = glGetError(); if (err) std::cout << "GL ERROR:: render_scene_sans_depth_blur() E err "
                                                      << err << std::endl;

               // std::cout << "DEBUG:: render_scene_sans_depth_blur() -- done -- " << std::endl;
            }

            di.draw_particles(eye);

            // 20220504-PE Either way round there is text texture masking problems over angry diego.
            //             Is the atom label shader not discarding the pixels correctly?
            draw_bad_nbc_atom_pair_markers(eye, PASS_TYPE_STANDARD);
            draw_molecules_atom_labels(eye);
            draw_rotation_centre_crosshairs(eye, GTK_GL_AREA(gl_area), PASS_TYPE_STANDARD);
            if (show_fps_flag)
               draw_hud_fps();

         }
      }
   }
}

void
graphics_info_t::render_scene_with_depth_blur(stereo_eye_t eye,
                                              Shader *shader_for_tmeshes_p, Shader *shader_for_meshes_p,
                                              Shader *shader_for_tmeshes_with_shadows_p,
                                              Shader *shader_for_meshes_with_shadows_p,
                                              int width, int height) {

   auto renderQuad = [] () {
                        if (quadVAO == 0) {
                           float quadVertices[] = {
                                                   // positions        // texture Coords
                                                   -1.0f,  1.0f, 0.0f, 0.0f, 1.0f,
                                                   -1.0f, -1.0f, 0.0f, 0.0f, 0.0f,
                                                    1.0f,  1.0f, 0.0f, 1.0f, 1.0f,
                                                    1.0f, -1.0f, 0.0f, 1.0f, 0.0f,
                           };
                           // setup plane VAO
                           glGenVertexArrays(1, &quadVAO);
                           glGenBuffers(1, &quadVBO);
                           glBindVertexArray(quadVAO);
                           glBindBuffer(GL_ARRAY_BUFFER, quadVBO);
                           glBufferData(GL_ARRAY_BUFFER, sizeof(quadVertices), &quadVertices, GL_STATIC_DRAW);
                           glEnableVertexAttribArray(0);
                           glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void*)0);
                           glEnableVertexAttribArray(1);
                           glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void*)(3 * sizeof(float)));
                        }
                        GLenum err = glGetError();
                        if (err)
                           std::cout << "GL ERROR:: render_scene_with_depth_blur() lambda renderQuad() A " << err << std::endl;
                        glBindVertexArray(quadVAO);
                        err = glGetError();
                        if (err)
                           std::cout << "GL ERROR:: render_scene_with_depth_blur() lambda renderQuad() B " << err << std::endl;
                        glEnableVertexAttribArray(0);
                        glEnableVertexAttribArray(1);
                        glBindVertexArray(quadVAO);
                        glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
                        glBindVertexArray(0);
                     };

   auto render_to_shadow_map = [] () {

                                 if (shadow_strength == 0.0) return;
                                  graphics_info_t di;
                                  glViewport(0, 0, shadow_texture_width, shadow_texture_height);
                                  glBindFramebuffer(GL_FRAMEBUFFER, shadow_depthMap_framebuffer);

                                  glClear(GL_DEPTH_BUFFER_BIT);
                                  unsigned int light_index = 0;
                                  // make these static at some stage
                                  di.draw_Models_for_shadow_map(light_index); // for coot non-Model objects
                                  di.draw_molecules_for_shadow_map(light_index);

                               };

   graphics_info_t di;
   GtkWidget *gl_area = graphics_info_t::glareas[0];

   int ssao_blur_size = graphics_info_t::ssao_blur_size;

   if (di.displayed_image_type == SHOW_AO_SCENE) {

      // std::cout << "DEBUG:: render_scene_with_depth_blur() ------------------------------- " << std::endl;

      // bind SSAO gbuffer framebuffer
      // draw_models_for_ssao()
      // bind SSAO framebuffer
      // renderQuad()
      // bind SSAO blur framebuffer
      // renderQuad()
      // bind depthMap framebuffer // call this shadow_depthMap
      // draw_models_for_shadow_map()
      // bind effects framebuffer
      // draw_models_with_shadows()
      // bind blur_y framebuffer
      // ... send ssao texture/flag/strength
      // glDrawArrays()
      // bind blur_x framebuffer
      // render_scene_with_y_blur()
      // bind combine_textures_using_depth_framebuffer
      // render_scene_with_x_blur();
      // di.attach_buffers()
      // render_scene_with_texture_combination_for_depth_blur();

      {

         // GtkAllocation allocation;
         // gtk_widget_get_allocation(di.gl_area, &allocation);
         int w = width;
         int h = height;

         // render
         // ------

         // 1. geometry pass: render scene's geometry/color data into gbuffer
         // -----------------------------------------------------------------
         // glBindFramebuffer(GL_FRAMEBUFFER, di.gBufferFBO);

         // std::cout << "render_scene_sans_depth_blur() here 1 " << std::endl;

         di.framebuffer_for_ssao_gbuffer.bind();
         glClearColor(0.0, 0.0, 0.0, 1.0);
         glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
         bool do_orthographic_projection = ! perspective_projection_flag;
         glm::mat4 projection_matrix = di.get_projection_matrix(do_orthographic_projection, w, h);
         glm::mat4 view_matrix       = di.get_view_matrix();
         glm::mat4 model_matrix      = di.get_model_matrix();
         model_matrix = glm::translate(model_matrix, glm::vec3(1.1f, 4.0f, 2.2f));
         model_matrix = glm::scale(model_matrix, glm::vec3(0.5f, 0.5f, 0.5f));
         di.shaderGeometryPass.Use();
         di.shaderGeometryPass.set_mat4_for_uniform("model", model_matrix);
         di.shaderGeometryPass.set_mat4_for_uniform("view", view_matrix);
         di.shaderGeometryPass.set_mat4_for_uniform("projection", projection_matrix);

         if (false) {
            std::cout << "sending model      " << glm::to_string(model_matrix) << std::endl;
            std::cout << "sending view       " << glm::to_string(view_matrix) << std::endl;
            std::cout << "sending projection " << glm::to_string(projection_matrix) << std::endl;
         }

         di.shaderGeometryPass.set_int_for_uniform("invertedNormals", 0);

         model_matrix = di.get_model_matrix();
         model_matrix = glm::translate(model_matrix, glm::vec3(0.0f, 0.5f, 0.0));
         model_matrix = glm::rotate(model_matrix, glm::radians(-90.0f), glm::vec3(1.0, 0.0, 0.0));
         model_matrix = glm::scale(model_matrix, glm::vec3(1.0f));
         di.shaderGeometryPass.set_mat4_for_uniform("model", model_matrix);

         // render_3d_scene_for_ssao(); // 20220130-PE render_3d_scene() was coot before Models/crows
         draw_models_for_ssao();
         draw_molecules_for_ssao(eye);

         glBindFramebuffer(GL_FRAMEBUFFER, 0);


         // 2. generate SSAO texture
         // ------------------------
         glBindFramebuffer(GL_FRAMEBUFFER, di.ssaoFBO);
         // std::cout << "Here in render_scene_sans_depth_blur() bound ssaFBO frame buffer " << di.ssaoFBO << std::endl;

         glClear(GL_COLOR_BUFFER_BIT);
         di.shaderSSAO.Use();
         di.shaderSSAO.set_int_for_uniform("gPosition", 0);
         di.shaderSSAO.set_int_for_uniform("gNormal",   1);
         di.shaderSSAO.set_int_for_uniform("texNoise",  2);

         // Send kernel + rotation
         di.shaderSSAO.set_int_for_uniform("n_ssao_kernel_samples", di.n_ssao_kernel_samples);
         for (unsigned int i = 0; i < di.n_ssao_kernel_samples; ++i) {
            // std::cout << "kernel sample " << i << " " << glm::to_string(di.ssaoKernel[i]) << std::endl;
            di.shaderSSAO.set_vec3_for_uniform("samples[" + std::to_string(i) + "]", di.ssaoKernel[i]);
         }
         di.shaderSSAO.set_mat4_for_uniform("projection", projection_matrix);

         if (false)
            std::cout << "sending ssao radius " << SSAO_radius << " bias " << SSAO_bias
                      << " n-kernel-samples "  << n_ssao_kernel_samples << std::endl;
         shaderSSAO.set_float_for_uniform("radius", SSAO_radius);
         shaderSSAO.set_float_for_uniform("bias",   SSAO_bias);
         glActiveTexture(GL_TEXTURE0);
         glBindTexture(GL_TEXTURE_2D, di.framebuffer_for_ssao_gbuffer.gPosition);
         glActiveTexture(GL_TEXTURE1);
         glBindTexture(GL_TEXTURE_2D, di.framebuffer_for_ssao_gbuffer.gNormal);
         glActiveTexture(GL_TEXTURE2);
         glBindTexture(GL_TEXTURE_2D, di.noiseTexture);
         if (false)
            std::cout << "debug:: in render_scene() gPosition is " << di.framebuffer_for_ssao_gbuffer.gPosition
                      << " and gNormal is " << di.framebuffer_for_ssao_gbuffer.gNormal << std::endl;
         renderQuad();
         glBindFramebuffer(GL_FRAMEBUFFER, 0);

         // std::cout << "debug gPosition gNormal noiseTexture " << di.gPosition << " "
         //           << di.gNormal << " " << di.noiseTexture << std::endl;


         // 3. blur SSAO texture to remove noise
         // ------------------------------------
         glBindFramebuffer(GL_FRAMEBUFFER, di.ssaoBlurFBO);
         // std::cout << "Here in render_scene_sans_depth_blur() bound ssaFBOblur frame buffer " << di.ssaoBlurFBO
         // << std::endl;
         glClear(GL_COLOR_BUFFER_BIT);
         di.shaderSSAOBlur.Use();
         di.shaderSSAOBlur.set_int_for_uniform("ssaoInput", 0);
         di.shaderSSAOBlur.set_int_for_uniform("blur_size", ssao_blur_size);
         glActiveTexture(GL_TEXTURE0);
         glBindTexture(GL_TEXTURE_2D, di.ssaoColorBuffer);
         renderQuad();
         glBindFramebuffer(GL_FRAMEBUFFER, 0);
         // std::cout << "render_scene_sans_depth_blur() here 4 " << std::endl;

         bool show_shadow_map = false;

         if (show_shadow_map) {

            render_to_shadow_map();

            // now we have a texture in a framebuffer. I want to see what that looks like
            di.attach_buffers(); // switch back to the defaut framebuffer
            glViewport(0, 0, graphics_x_size, graphics_y_size);
            float gamma_pe = 1.2f;
            std::cout << "debug:: clearing with gamma colour" << std::endl;
            glClearColor(pow(0.07f, gamma_pe), pow(0.13f, gamma_pe), pow(0.17f, gamma_pe), 1.0f);
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

            glActiveTexture(GL_TEXTURE0 + 0);
            glBindTexture(GL_TEXTURE_2D, shadow_depthMap_texture);
            di.tmesh_for_shadow_map.draw(&di.shader_for_shadow_map_image_texture_mesh);

         } else {

            render_to_shadow_map();

            {
               glViewport(0,0, width, height);
               di.framebuffer_for_effects.bind();
               // di.attach_buffers(); // hack
               glClearColor(background_colour.r, background_colour.g, background_colour.b, 1.0);
               glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
               glEnable(GL_DEPTH_TEST);
               glDepthFunc(GL_LESS);
               glDisable(GL_BLEND);
               bool draw_shadows = true;

               // std::cout << "show_just_shadows " << show_just_shadows << std::endl;

               // maybe doing the shadows separately is the right approach, in fact.
               // because currently *every* draw call of each of the tmeshes
               // samples the shadow map (for every pixel). And the shadow map
               // samples a 3x3 patch for each pixel.
               // Maybe it's better to draw (like SHOW_BASIC_SCENE) to a framebuffer
               // colour texture without shadows and then use a new shader to combine
               // that image with the shadow map.
               //
               // I will write a shader for the meshes (shader_for_meshes_with_shadows_p)
               // that doesn't draw shadows.

               di.draw_models_with_shadows(eye,
                                           shader_for_tmeshes_with_shadows_p,
                                           shader_for_meshes_with_shadows_p,
                                           width, height,
                                           draw_shadows,
                                           shadow_strength,
                                           show_just_shadows);
               render_3d_scene_with_shadows(eye); // no longer does a glClear()
            }


            // if I enable di.attach_buffers() above and comment out below, then it renders

            // di.attach_buffers(); // for non-blurring
            di.blur_y_framebuffer.bind();

            GLint local_fbo;
            glGetIntegerv(GL_FRAMEBUFFER_BINDING, &local_fbo);

            {
               glClearColor(background_colour.r, background_colour.g, background_colour.b, 1.0);
               glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
               glBindVertexArray(di.screen_AO_quad_vertex_array_id); // maybe make this a HUDTexture?
               glm::vec4 bg_col(background_colour, 1.0f);

               glActiveTexture(GL_TEXTURE0 + 0);
               glBindTexture(GL_TEXTURE_2D, di.framebuffer_for_effects.get_texture_colour());
               glActiveTexture(GL_TEXTURE0 + 1);
               glBindTexture(GL_TEXTURE_2D, di.framebuffer_for_effects.get_texture_depth());
               glActiveTexture(GL_TEXTURE0 + 2);
               glBindTexture(GL_TEXTURE_2D, di.ssaoColorBufferBlur);
               di.shader_for_effects.Use();
               // std::cout << "using shader_for_effects " << shader_for_effects.name << std::endl;
               shader_for_effects.set_int_for_uniform("screenTexture", 0);
               shader_for_effects.set_int_for_uniform("screenDepth",   1);
               shader_for_effects.set_int_for_uniform("ssao",          2); // sampler2D
               shader_for_effects.set_bool_for_uniform("use_ssao", di.use_ssao);
               shader_for_effects.set_float_for_uniform("ssao_strength", di.ssao_strength);
               shader_for_effects.set_bool_for_uniform("show_ssao", di.show_just_ssao);
               shader_for_effects.set_vec4_for_uniform("background_colour", bg_col);
               shader_for_effects.set_int_for_uniform("effects_output_type", effects_shader_output_type);
               shader_for_effects.set_float_for_uniform("brightness", effects_brightness);
               shader_for_effects.set_float_for_uniform("gamma",      effects_gamma);
               GLenum err = glGetError();
               if (err)
                  std::cout << "GL ERROR:: render_scene_sans_depth_blur() D err " << err << std::endl;

               glEnableVertexAttribArray(0);
               glEnableVertexAttribArray(1);
               glDrawArrays(GL_TRIANGLES, 0, 6);
               err = glGetError(); if (err) std::cout << "GL ERROR:: render_scene_sans_depth_blur() E err "
                                                      << err << std::endl;

               // std::cout << "DEBUG:: render_scene_sans_depth_blur() -- done -- " << std::endl;
            }

            di.draw_particles(eye);

            draw_bad_nbc_atom_pair_markers(eye, PASS_TYPE_STANDARD);
            draw_molecules_atom_labels(eye);
            draw_rotation_centre_crosshairs(eye, GTK_GL_AREA(gl_area), PASS_TYPE_STANDARD);

            blur_x_framebuffer.bind();
            render_scene_with_y_blur();

            combine_textures_using_depth_framebuffer.bind();
            render_scene_with_x_blur();

            gtk_gl_area_attach_buffers(GTK_GL_AREA(gl_area));

            render_scene_with_texture_combination_for_depth_blur();

            // this does not belong here - is it rendered twice?
            if (di.show_fps_flag)
               di.draw_hud_fps();

        }
      }
   }
}



// needs args bool to_screendump_framebuffer_flag, const std::string &output_file_name
gboolean
graphics_info_t::render_scene() {

   bool do_mono = true;
   if (display_mode == coot::SIDE_BY_SIDE_STEREO)          do_mono = false;
   if (display_mode == coot::SIDE_BY_SIDE_STEREO_WALL_EYE) do_mono = false;

   if (do_mono) {
      return render_scene_for_eye_internal(stereo_eye_t::MONO);
   } else {
      bool ss = render_scene_for_eye_internal(stereo_eye_t::LEFT_EYE);
      bool st = render_scene_for_eye_internal(stereo_eye_t::RIGHT_EYE);
      return ss;
   }
}


// needs args bool to_screendump_framebuffer_flag, const std::string &output_file_name
gboolean
graphics_info_t::render_scene_for_eye_internal(stereo_eye_t eye) {

   // make this a member data
   // bool draw_background_image = true;

   auto render_scene_basic = [] (stereo_eye_t eye) {

      // std::cout << "--- render_scene_basic() ---------------------------------------- " << std::endl;
      // auto tp_0 = std::chrono::high_resolution_clock::now();

      GLenum err = glGetError();
      if (err)
         std::cout << "GL ERROR:: render_scene_basic() --- start --- " << err << std::endl;

      GtkAllocation allocation;
      auto gl_area = graphics_info_t::glareas[0];
      gtk_widget_get_allocation(gl_area, &allocation);
      int width = allocation.width;
      int height = allocation.height;
      // auto tp_1 = std::chrono::high_resolution_clock::now();
      // auto d10 = std::chrono::duration_cast<std::chrono::microseconds>(tp_1 - tp_0).count();
      // std::cout << "Timings for window allocation " << d10 << " microseconds" << std::endl;  // 0
      int sf = 1;
#ifdef __APPLE__
      sf = 2;
#endif
      // we always want this viewport to be the size of the widget (in the case of APPLE, there
      // is the double resolution issue to handle)

      if (scale_up_graphics != 1) {
         width  *= scale_up_graphics;
         height *= scale_up_graphics;
      }
      if (scale_down_graphics != 1) {
         width /= scale_down_graphics;
         height /= scale_down_graphics;
      }
      // std::cout << "render_scene_basic() " << width << " " << height << std::endl;

      err = glGetError();
      if (err)
         std::cout << "GL ERROR:: render_scene_basic() A " << err << std::endl;

      attach_buffers(); // just GTK things
      err = glGetError();
      if (err)
         std::cout << "GL ERROR:: render_scene_basic() B " << err << std::endl;

      int left = 0;
      int bot = 0;
      int right = width * sf;
      int top = height *sf;
      if (eye == stereo_eye_t::LEFT_EYE) {
         right /= 2;
      }
      if (eye == stereo_eye_t::RIGHT_EYE) {
         left = right/2;
      }
      glViewport(left, bot, right, top);
      err = glGetError();
      if (err)
         std::cout << "GL ERROR:: render_scene_basic() C " << err << std::endl;

      glClearColor(background_colour.r, background_colour.g, background_colour.b, 1.0);
      err = glGetError();
      if (err)
         std::cout << "GL ERROR:: render_scene_basic() D " << err << std::endl;

      if (eye != stereo_eye_t::RIGHT_EYE) // drawn second
         glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
      err = glGetError();
      if (err)
         std::cout << "GL ERROR:: render_scene_basic() E " << err << std::endl;

      glDisable(GL_BLEND);
      err = glGetError();
      if (err)
         std::cout << "GL ERROR:: render_scene_basic() F " << err << std::endl;

      glEnable(GL_DEPTH_TEST);
      err = glGetError();
      if (err)
         std::cout << "GL ERROR:: render_scene_basic() G " << err << std::endl;

      glDepthFunc(GL_LESS);

      if (draw_background_image_flag) {
         // std::cout << "basic render: draw background image " << std::endl;
         texture_for_background_image.Bind(0);
         tmesh_for_background_image.draw(&shader_for_background_image, HUDTextureMesh::TOP_LEFT);
      }

      graphics_info_t g;// needed? Yes.
      err = glGetError();
      if (err)
         std::cout << "GL ERROR:: render_scene_basic() H " << err << std::endl;

      g.draw_models(eye, &shader_for_tmeshes, &shader_for_meshes, nullptr, nullptr, width, height);
      err = glGetError();
      if (err)
         std::cout << "GL ERROR:: render_scene_basic() I " << err << std::endl;

      draw_rotation_centre_crosshairs(eye, GTK_GL_AREA(gl_area), PASS_TYPE_STANDARD);
      render_3d_scene(GTK_GL_AREA(gl_area), eye);
      // HUD things? Not here?
      if (show_fps_flag) {
         err = glGetError();
         if (err)
            std::cout << "GL ERROR:: render_scene_basic() J " << err << std::endl;
         draw_hud_fps();
      }
   };

   // 20250422-PE - if we switch contexts when rendering the Rama plot, we needd
   // to switch back. It seems that here is the place to do it (I might be wrong).
   //
   GtkWidget *gl_area = glareas[0];
   gtk_gl_area_make_current(GTK_GL_AREA(gl_area));
   if (gtk_gl_area_get_error(GTK_GL_AREA(gl_area))) {
      auto mess = gtk_gl_area_get_error(GTK_GL_AREA(gl_area))->message;
      std::cout << "ERROR:: graphics_info_t::render_scene(): error making current GL context: "
                << mess << std::endl;
   };
   // it is disasterous to carry on here if we get the above error

   // crow variable conversion
   Shader *shader_for_tmeshes_p = &shader_for_texture_meshes;
   Shader *shader_for_meshes_p  = &shader_for_meshes;
   Shader *shader_for_tmeshes_with_shadows_p = &shader_for_tmeshes_with_shadows;
   Shader *shader_for_meshes_with_shadows_p  = &shader_for_meshes_with_shadows;

   gboolean status = gboolean(true);

   bool show_basic_scene_state = (displayed_image_type == SHOW_BASIC_SCENE);

   // std::cout << "render_scene(): show_basic_scene_state " << show_basic_scene_state << std::endl;

   if (show_basic_scene_state) {

      GLenum err = glGetError();
      if (err)
         std::cout << "GL ERROR:: render_scene() basic path " << err << std::endl;

      render_scene_basic(eye);
   } else {
      if (shader_do_depth_of_field_blur_flag || shader_do_outline_flag) {
         render_scene_with_depth_blur(eye,
                                      shader_for_tmeshes_p, // or outline
                                      shader_for_meshes_p,
                                      shader_for_tmeshes_with_shadows_p,
                                      shader_for_meshes_with_shadows_p,
                                      graphics_x_size, graphics_y_size);
      } else {
         render_scene_sans_depth_blur(eye,
                                      shader_for_tmeshes_p,
                                      shader_for_meshes_p,
                                      shader_for_tmeshes_with_shadows_p,
                                      shader_for_meshes_with_shadows_p,
                                      graphics_x_size, graphics_y_size);
      }
   }

   return status;
}

