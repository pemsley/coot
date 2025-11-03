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

#ifdef WITH_SOUND
#include <mutex>
class OpenALState {
   public:
   ALCdevice *device;
   ALCcontext *context;

   OpenALState(ALCdevice * m_device, ALCcontext* m_context) : device(m_device), context(m_context) {
      // empty
   }

   ~OpenALState() {
      if(context) {
         alcMakeContextCurrent(NULL);
         alcDestroyContext(context);
      }
      if(device) {
         alcCloseDevice(device);
      }
   }
};

inline std::mutex openal_init_mutex;
inline OpenALState* openal_state_ptr = nullptr;

#endif

void
play_sound_file(const std::string &file_name) {

#ifdef WITH_SOUND

   std::cout << "play_sound_file() " << file_name << std::endl;

   auto print_alerror = [] (ALenum err) {
     std::string s = std::to_string(err);
     if (err == AL_NO_ERROR)          s = "AL_NO_ERROR";
     if (err == AL_INVALID_NAME)      s = "AL_INVALID_NAME";
     if (err == AL_INVALID_ENUM)      s = "AL_INVALID_ENUM";
     if (err == AL_INVALID_VALUE)     s = "AL_INVALID_VALUE";
     if (err == AL_INVALID_OPERATION) s = "AL_INVALID_OPERATION";
     if (err == ALC_INVALID_CONTEXT) s = "ALC_INVALID_CONTEXT";
     if (err == ALC_INVALID_DEVICE) s = "ALC_INVALID_DEVICE";
     return s;
   };

   auto check_alerror = [print_alerror] (const std::string &where) {
      ALenum err = alGetError();
      if (err != AL_NO_ERROR) {
         std::cout << "AL ERROR:: " << where << " : " << print_alerror(err) << std::endl;
         return true;
      }
      return false;
   };

   auto play_sound_file_inner = [print_alerror, check_alerror] (const std::string &file_name) {

      std::cout << "DEBUG:: play_sound_file_inner: " << file_name << std::endl;

      ALCdevice *m_pDevice;
      ALCcontext *m_pContext;

      openal_init_mutex.lock();
      if(openal_state_ptr != nullptr) {
         m_pDevice = openal_state_ptr->device;
         m_pContext = openal_state_ptr->context;
      } else { // Initialize OpenAL
         m_pDevice = alcOpenDevice(NULL);
         if(!m_pDevice) {
            std::cout << "ERROR:: play_sound_file_inner() giving up after device opening failure." << std::endl;
            return;
         }
         m_pContext = alcCreateContext(m_pDevice, NULL);
         // this error check seems to be unreliable
         // check_alerror("Context creation");
         if(!m_pContext) {
            // We don't close the device here, as it may be in use by other threads
            std::cout << "ERROR:: play_sound_file_inner() giving up after context creation failure." << std::endl;
            return;
         }

         alcMakeContextCurrent(m_pContext);
         check_alerror("make context current");
         
         OpenALState* initial_state = new OpenALState(m_pDevice, m_pContext);
         openal_state_ptr = initial_state;

         const auto* default_device = alcGetString(NULL, ALC_DEFAULT_DEVICE_SPECIFIER);
         if(default_device) {
            std::cout << "DEBUG:: Default sound device: " << default_device;
         } else {
            std::cout << "DEBUG:: No default sound device found. ";
         }
         const auto* devices = alcGetString(NULL, ALC_DEVICE_SPECIFIER);
         std::cout << "; Available sound devices: ";
         while(*devices) {
            std::cout << " " << devices;
            devices += strlen(devices) + 1;
         }
         std::cout << std::endl;
      }
      openal_init_mutex.unlock();


      std::cout << "debug:: m_pDevice is " << m_pDevice << std::endl;
      
      if(!m_pDevice || !m_pContext) {
         std::cout << "ERROR:: play_sound_file_inner() device or context is null. Giving up" << std::endl;
         return;
      }
      
      ALuint source;
      alGenSources(1, &source);
      check_alerror("source generation");
      
      ALuint buffer;
      alGenBuffers(1, &buffer);
      check_alerror("buffer generation");

      auto openal_cleanup = [&] () {
         alDeleteBuffers(1, &buffer);
         alDeleteSources(1, &source);
      };

      FILE* file = fopen(file_name.c_str(), "rb");

      if(!file) {
         std::cout << "ERROR:: play_sound_file_inner() could not open sound file. Giving up. File: " << file_name << std::endl;
         openal_cleanup();
         return;
      }

      OggVorbis_File ovf;
      if (ov_open(file, &ovf, NULL, 0) < 0) {
         std::cout << "ERROR:: play_sound_file_inner() could not open OggVorbis file. Giving up. File: " << file_name << std::endl;
         fclose(file);
         openal_cleanup();
         return;
      }

      vorbis_info *vi = ov_info(&ovf, -1);

      ALenum format = vi->channels == 1 ? AL_FORMAT_MONO16 : AL_FORMAT_STEREO16;
      ALsizei size_in_samples = ov_pcm_total(&ovf, -1) * vi->channels;
      // Both AL_FORMAT_MONO16 and AL_FORMAT_STEREO16 use two bytes per sample.
      ALsizei size_in_bytes = size_in_samples * 2;
      ALshort* data = new ALshort[size_in_samples];
      int bitstream = 0;
      unsigned int offset = 0;
      while(true) {
         long ov_read_res = ov_read(&ovf, (char*)(data) + offset, size_in_bytes - offset, 0, 2, 1, &bitstream);
         if (ov_read_res > 0) {
            offset += ov_read_res;
         } else if (ov_read_res == 0) {
            // EOF
            break;
         } else {
            std::cout << "ERROR:: play_sound_file_inner() error while reading OggVorbis file. Giving up. File: " << file_name << std::endl;
            ov_clear(&ovf);
            openal_cleanup();
            delete[] data;
            return;
         }
      }

      alBufferData(buffer, format, data, size_in_bytes, vi->rate);
      check_alerror("Write buffer data");

      alSourcei(source, AL_BUFFER, buffer);
      check_alerror("Attach buffer to source");
      std::cout << "DEBUG:: Playing the sound " << file_name << std::endl;
      alSourcePlay(source);
      check_alerror("source play");

      ALint source_state;
      do {
         // std::this_thread::yield();
         std::this_thread::sleep_for(std::chrono::milliseconds(100));
         alGetSourcei(source, AL_SOURCE_STATE, &source_state);
      } while (source_state == AL_PLAYING);

      std::cout << "DEBUG:: Playback finished " << file_name << std::endl;

      ov_clear(&ovf);
      // this appear to be already freed above
      // fclose(file);
      openal_cleanup();
      delete[] data;
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


