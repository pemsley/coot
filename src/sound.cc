
#ifdef WITH_SOUND

#include <string>
#include <iostream>

#include <vorbis/vorbisfile.h>

int test_sound(int argc, char **argv) {

   int status = 0;
   std::string fn("test.ogg");
   FILE *f = fopen(fn.c_str(), "r");
   if (!f) {
      std::cout << "DEBUG:: test_sound(): File " << fn << " could not be found." << std::endl;
   } else {
      OggVorbis_File ovf;
      const char *initial = 0;
      long ibytes = 1;
      int open_status = ov_open(f, &ovf, initial, ibytes);
      if (open_status < 0) {
	 std::cout << "Failed to VorbisOpen " << fn << std::endl;
      } else {
	 if(!ov_seekable(&ovf)) {
	    std::cout << "Failed to Vorbis Seek " << fn << std::endl;
	 } else {
	    long tt = ov_time_total(&ovf, -1);
	    std::cout << "OggVorbis total time for " << fn << " " << tt << std::endl;

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

	       while(!eof){
		  long ret= ov_read(&ovf, pcmout, sizeof(pcmout), 0, 2, 1, &current_section);
		  if (ret == 0) {
		     /* EOF */
		     eof=1;
		  } else if (ret < 0) {
		     std::cout << "test_sound(): error in the stream" << std::endl;
		  } else {
		     // fwrite(pcmout,1,ret,stdout);
		  }
	       }
	       
	    }
	 }
      }

      fclose(f);

      // close oggvorbisfile?
   }

   return(0);

}

#endif // WITH_SOUND
