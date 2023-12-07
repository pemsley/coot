
#include <string>
#include <iostream>
#include <thread>

// This will need to be more clever in future
#ifdef WITH_SOUND

#include <vorbis/vorbisfile.h>
#ifdef __APPLE__
#include <AL/al.h>
#include <AL/alc.h>
#else
#include <alsa/asoundlib.h>
#endif

#endif

#include "utils/coot-utils.hh"
#include "graphics-info.h" // to check for graphics_info_t::use_sounds

#include <atomic>

std::atomic<unsigned int> n_sound_files_playing(0);

void
play_sound_file_macos(const std::string &file_name) {
#ifdef WITH_SOUND
#ifdef __APPLE__
   auto play_sound_file_inner = [] (const std::string &file_name) {

      ALuint source;
      alGenSources(1, &source);

      ALuint buffer;
      alGenBuffers(1, &buffer);
      FILE* file = fopen(file_name.c_str(), "rb");
      OggVorbis_File ovf;
      ov_open(file, &ovf, NULL, 0);
      vorbis_info *vi = ov_info(&ovf, -1);

      ALsizei size = vi->channels * vi->rate * 2;
      ALshort* data = new ALshort[size];
      int bitstream = 0;
      ov_read(&ovf, (char*)data, size, 0, 2, 1, &bitstream);
      alBufferData(buffer, AL_FORMAT_STEREO16, data, size, vi->rate);

      // Play the sound
      alSourcei(source, AL_BUFFER, buffer);
      alSourcePlay(source);

      ov_clear(&ovf);
      fclose(file);
   };

   if (!coot::file_exists(file_name))
      return;

   std::thread t(play_sound_file_inner, file_name);
   t.detach();
#endif
#endif

}

void
play_sound_file(const std::string &file_name) {

#ifdef WITH_SOUND

#ifdef __APPLE__
#else

   auto play_sound_file_inner = [] (const std::string &file_name) {

      FILE *f = fopen(file_name.c_str(), "r");
      if (!f) {
         std::cout << "DEBUG:: test_sound(): File " << file_name << " could not be found." << std::endl;
      } else {
         OggVorbis_File ovf;
         const char *initial = 0;
         long ibytes = 1;
         int open_status = ov_open(f, &ovf, initial, ibytes);
         if (open_status < 0) {
            std::cout << "Failed to VorbisOpen " << file_name << std::endl;
         } else {
            if (true) {
               n_sound_files_playing++;
               std::string ss = "n_sound_files_playing: " + std::to_string(n_sound_files_playing) + "\n";
               std::cout << ss << std::endl;
               if(!ov_seekable(&ovf)) {
                  std::cout << "Failed to Vorbis Seek " << file_name << std::endl;
               } else {

                  snd_pcm_sframes_t frames;
                  snd_pcm_t *handle;
                  const char *device = "default";

                  int err = snd_pcm_open(&handle, device, SND_PCM_STREAM_PLAYBACK, 0);
                  if (err < 0) {
                     printf("Playback open error: %s\n", snd_strerror(err));
                     return;
                  }
                  err = snd_pcm_set_params(handle, SND_PCM_FORMAT_S16_LE, SND_PCM_ACCESS_RW_INTERLEAVED, 2, 44100, 1, 500000);  /* 0.5sec */
                  if (err < 0) {
                     printf("Playback open error: %s\n", snd_strerror(err));
                     return;
                  }

                  long tt = ov_time_total(&ovf, -1);
                  std::cout << "OggVorbis total time for " << file_name << " " << tt << std::endl;

                  for(int i=0; i<ov_streams(&ovf); i++){
                     vorbis_info *vi=ov_info(&ovf,i);
                     printf("    logical bitstream section %d information:\n", i+1);
                     printf("        %ldHz %d channels bitrate %ldkbps serial number=%ld\n",
                            vi->rate,vi->channels,ov_bitrate(&ovf,i)/1000,
                            ov_serialnumber(&ovf,i));
                     printf("        compressed length: %ld bytes ",(long)(ov_raw_total(&ovf,i)));
                     printf(" play time: %lds\n",(long)ov_time_total(&ovf,i));

                     int eof = 0;
                     char pcmout[4096];
                     int current_section;
                     while (!eof) {
                        // std::cout << "################ really playing sound" << std::endl;
                        long ret= ov_read(&ovf, pcmout, sizeof(pcmout), 0, 2, 1, &current_section);
                        if (ret == 0) {
                           /* EOF */
                           eof=1;
                        } else if (ret < 0) {
                           std::cout << "test_sound(): error in the stream" << std::endl;
                        } else {

                           frames = snd_pcm_writei(handle, pcmout, ret/4);
                           snd_pcm_wait(handle, 20000); // max wait time in ms.
                           if (frames < 0)
                              frames = snd_pcm_recover(handle, frames, 0);
                           if (frames < 0) {
                              printf("snd_pcm_writei failed: %s\n", snd_strerror(err));
                              break;
                           }
                           // test for dropped frames here, comparing ret and frames.

                           // fwrite(pcmout,1,ret,stdout);
                        }
                     }
                  }
               }
            }
            std::string ss = "reducing n_sound_files_playing\n";
            std::cout << ss << std::endl;
            n_sound_files_playing--;
         }
         fclose(f);
         // close oggvorbisfile?
         // ov_clear(&ovf);
      }
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

#endif
#ifdef __APPLE__
   play_sound_file_macos(file_name);
#else
   std::thread t(play_sound_file_inner, fn);
   t.detach();
#endif

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


