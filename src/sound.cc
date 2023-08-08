
#include <string>
#include <iostream>
#include <thread>

#include <vorbis/vorbisfile.h>
#include <alsa/asoundlib.h>

void
play_sound_file(const std::string &file_name) {

#ifdef WITH_SOUND

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
         fclose(f);
         // close oggvorbisfile?
      }
   };

   std::thread t(play_sound_file_inner, file_name);
   t.detach();
#endif

}

int test_sound(int argc, char **argv) {

   int status = 0;
   std::string fn("test.ogg");
   std::cout << "################ playing sound " << fn << std::endl;
   play_sound_file(fn);
   return status;
}


