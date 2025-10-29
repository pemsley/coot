/*
 * src/sound.cc
 *
 * Copyright 2018 by Medical Research Council
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

#include <string>
#include <iostream>
#include "sound.hh"

// This will need to be more clever in future
#ifdef WITH_SOUND

#include <vorbis/vorbisfile.h>
#include <AL/al.h>
#include <AL/alc.h>
#endif

#include "utils/coot-utils.hh"
#include "graphics-info.h" // to check for graphics_info_t::use_sounds

#include <thread>
#include <atomic>

std::atomic<unsigned int> n_sound_files_playing(0);

void
play_sound_file(const std::string &file_name) {

#ifdef WITH_SOUND

   std::cout << "play_sound_file_macos() " << file_name << std::endl;

   auto _ = [] (ALenum err) {
     std::string s = std::to_string(err);
     if (err == AL_NO_ERROR)          s = "AL_NO_ERROR";
     if (err == AL_INVALID_NAME)      s = "AL_INVALID_NAME";
     if (err == AL_INVALID_ENUM)      s = "AL_INVALID_ENUM";
     if (err == AL_INVALID_VALUE)     s = "AL_INVALID_VALUE";
     if (err == AL_INVALID_OPERATION) s = "AL_INVALID_OPERATION";
     return s;
   };

   auto play_sound_file_inner = [_] (const std::string &file_name) {

      std::cout << "DEBUG:: play_sound_file_inner: " << file_name << std::endl;

      ALCdevice *m_pDevice = alcOpenDevice(NULL);
      alcCreateContext(m_pDevice, NULL);

      std::cout << "debug:: m_pDevice is " << m_pDevice << std::endl;

      ALenum err = alGetError();
      if (err) std::cout << "AL ERROR:: play_sound_file_inner() A0 " << _(err) << std::endl;
      err = alGetError();
      if (err) std::cout << "AL ERROR:: play_sound_file_inner() A1 " << _(err) << std::endl;
      err = alGetError();
      if (err) std::cout << "AL ERROR:: play_sound_file_inner() A2 " << _(err) << std::endl;
      for (unsigned int i=0; i<10;i++) {
         err = alGetError();
         if (err) std::cout << "AL ERROR:: play_sound_file_inner() Ae " << i << " " << _(err) << std::endl;
      }
      ALuint source;
      alGenSources(1, &source);
      err = alGetError();
      if (err) std::cout << "AL ERROR:: play_sound_file_inner() B1 " << _(err) << std::endl;

      ALuint buffer;
      alGenBuffers(1, &buffer);
      err = alGetError();
      if (err) std::cout << "AL ERROR:: play_sound_file_inner() B2 " << _(err) << std::endl;
      FILE* file = fopen(file_name.c_str(), "rb");
      OggVorbis_File ovf;
      ov_open(file, &ovf, NULL, 0);
      vorbis_info *vi = ov_info(&ovf, -1);

      ALsizei size = vi->channels * vi->rate * 2;
      ALshort* data = new ALshort[size];
      int bitstream = 0;
      ov_read(&ovf, (char*)data, size, 0, 2, 1, &bitstream);
      err = alGetError();
      if (err) std::cout << "AL ERROR:: play_sound_file_inner() C " << _(err) << std::endl;
      alBufferData(buffer, AL_FORMAT_STEREO16, data, size, vi->rate);
      err = alGetError();
      if (err) std::cout << "AL ERROR:: play_sound_file_inner() D " << _(err) << std::endl;

      // Play the sound
      std::cout << "%%%%%%%%%%% play the sound " << file_name << std::endl;
      alSourcei(source, AL_BUFFER, buffer);
      err = alGetError();
      if (err) std::cout << "AL ERROR:: play_sound_file_inner() E " << _(err) << std::endl;
      alSourcePlay(source);
      err = alGetError();
      if (err) std::cout << "AL ERROR:: play_sound_file_inner() F " << _(err) << std::endl;

      ov_clear(&ovf);
      fclose(file);
   };

   std::string fn = file_name;
   if (coot::file_exists(fn)) {
      // don't touch the path then
   } else {
      // try to find it in the installation
      std::string dir = coot::package_data_dir();
      std::string sounds_dir = coot::util::append_dir_dir(dir, "sounds");
      fn = coot::util::append_dir_file(sounds_dir, fn);
   }

   if (!coot::file_exists(fn))
      return;

   std::thread t(play_sound_file_inner, fn);
   t.detach();
#endif
}

// void
// play_sound_file_alsa(const std::string &file_name) {

// #ifdef WITH_SOUND

//    auto play_sound_file_inner = [] (const std::string &file_name) {

//       bool debug = false;

//       FILE *f = fopen(file_name.c_str(), "r");
//       if (!f) {
//          std::cout << "DEBUG:: test_sound(): File " << file_name << " could not be found." << std::endl;
//       } else {
//          OggVorbis_File ovf;
//          const char *initial = 0;
//          long ibytes = 1;
//          int open_status = ov_open(f, &ovf, initial, ibytes);
//          if (open_status < 0) {
//             std::cout << "Failed to VorbisOpen " << file_name << std::endl;
//          } else {
//             if (true) {
//                n_sound_files_playing++;
//                std::string ss = "n_sound_files_playing: " + std::to_string(n_sound_files_playing) + "\n";
//                if (debug)
//                   std::cout << ss << std::endl;
//                if(!ov_seekable(&ovf)) {
//                   std::cout << "Failed to Vorbis Seek " << file_name << std::endl;
//                } else {

//                   snd_pcm_sframes_t frames;
//                   snd_pcm_t *handle;
//                   const char *device = "default";

//                   int err = snd_pcm_open(&handle, device, SND_PCM_STREAM_PLAYBACK, 0);
//                   if (err < 0) {
//                      printf("Playback open error: %s\n", snd_strerror(err));
//                      return;
//                   }
//                   err = snd_pcm_set_params(handle, SND_PCM_FORMAT_S16_LE, SND_PCM_ACCESS_RW_INTERLEAVED, 2, 44100, 1, 500000);  /* 0.5sec */
//                   if (err < 0) {
//                      printf("Playback open error: %s\n", snd_strerror(err));
//                      return;
//                   }

//                   long tt = ov_time_total(&ovf, -1);
//                   if (debug)
//                      std::cout << "OggVorbis total time for " << file_name << " " << tt << std::endl;

//                   for(int i=0; i<ov_streams(&ovf); i++){
//                      vorbis_info *vi=ov_info(&ovf,i);
//                      if (false) {
//                         printf("    logical bitstream section %d information:\n", i+1);
//                         printf("        %ldHz %d channels bitrate %ldkbps serial number=%ld\n",
//                                vi->rate,vi->channels,ov_bitrate(&ovf,i)/1000,
//                                ov_serialnumber(&ovf,i));
//                         printf("        compressed length: %ld bytes ",(long)(ov_raw_total(&ovf,i)));
//                         printf(" play time: %lds\n",(long)ov_time_total(&ovf,i));
//                      }

//                      int eof = 0;
//                      char pcmout[4096];
//                      int current_section;
//                      while (!eof) {
//                         // std::cout << "################ really playing sound" << std::endl;
//                         long ret= ov_read(&ovf, pcmout, sizeof(pcmout), 0, 2, 1, &current_section);
//                         if (ret == 0) {
//                            /* EOF */
//                            eof=1;
//                         } else if (ret < 0) {
//                            std::cout << "test_sound(): error in the stream" << std::endl;
//                         } else {

//                            frames = snd_pcm_writei(handle, pcmout, ret/4);
//                            snd_pcm_wait(handle, 20000); // max wait time in ms.
//                            if (frames < 0)
//                               frames = snd_pcm_recover(handle, frames, 0);
//                            if (frames < 0) {
//                               printf("snd_pcm_writei failed: %s\n", snd_strerror(err));
//                               break;
//                            }
//                            // test for dropped frames here, comparing ret and frames.

//                            // fwrite(pcmout,1,ret,stdout);
//                         }
//                      }
//                   }
//                }
//             }
//             std::string ss = "reducing n_sound_files_playing\n";
//             if (false)
//                std::cout << ss << std::endl;
//             n_sound_files_playing--;
//          }
//          fclose(f);
//          // close oggvorbisfile?
//          // ov_clear(&ovf);
//       }
//    };


//    std::string fn = file_name;
//    if (coot::file_exists(fn)) {
//       // don't touch the path then
//    } else {
//       // try to find it in the installation
//       std::string dir = coot::package_data_dir();
//       std::string sounds_dir = coot::util::append_dir_dir(dir, "sounds");
//       fn = coot::util::append_dir_file(sounds_dir, fn);
//    }
//    std::thread t(play_sound_file_inner, fn);
//    t.detach();

// #endif

// }

void play_sound(const std::string &type) {

   if (graphics_info_t::use_sounds) {
      if (type == "SUCCESS") play_sound_file("538554_3725923-lq-Sjonas88-success.ogg");
      if (type == "CLICK")   play_sound_file("538548_3725923-lq-Sjonas-Select-3.ogg");
      if (type == "TINK")    play_sound_file("538549_3725923-lq-Sjonas-Select-2.ogg");
      if (type == "STARS")   play_sound_file("538553_3725923-lq-Sjonas88-Stars.ogg");
      if (type == "OOPS")    play_sound_file("538550_3725923-lq-Sjonas88-Deep-tone.ogg");
      if (type == "diego-gone-pop") play_sound_file("pop-dodrio-554022_1433422-lq.ogg");
      if (type == "diego-arrives")  play_sound_file("cdonahueucsd-337132_5955158-lq.ogg");
   }

}


int test_sound(int argc, char **argv) {

   int status = 0;
   std::string fn("test.ogg");
   std::cout << "################ playing sound " << fn << std::endl;
   play_sound_file(fn);
   return status;
}


