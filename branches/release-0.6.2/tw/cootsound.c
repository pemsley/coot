/* src/cootsound.c
 * 
 * Copyright 2005 by Paul Emsley, The University of York
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA
 */

/* Taken from Ogg Vorbis Decoding Example Code */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "vorbisfile.h"
#include <ao/ao.h>

char pcmout[4096];

/* Inside main(), we declare our primary OggVorbis_File structure. We
   also declare a few other helpful variables to track out progress
   within the file. Also, we make our final concession to Windows
   users by setting the stdin and stdout to binary mode. */

int main(int argc, char **argv){
  OggVorbis_File vf;
  int eof=0;
  int current_section;

  /* Deal with libao functions */
  ao_initialize();
  int driver_id = ao_default_driver_id();
  ao_sample_format sample_format;
  ao_option *options = NULL;


/* ov_open() must be called to initialize the OggVorbis_File structure
   with default values. ov_open() also checks to ensure that we're
   reading Vorbis format and not something else. */

  if(ov_open(stdin, &vf, NULL, 0) < 0) {
      fprintf(stderr,"Input does not appear to be an Ogg bitstream.\n");
      exit(1);
  }

/* We're going to pull the channel and bitrate info from the file
   using ov_info() and show them to the user. We also want to pull out
   and show the user a comment attached to the file using
   ov_comment(). */

  {
    char **ptr=ov_comment(&vf,-1)->user_comments;
    vorbis_info *vi=ov_info(&vf,-1);
    while(*ptr){
      fprintf(stderr,"%s\n",*ptr);
      ++ptr;
    }
    fprintf(stderr,"\nBitstream is %d channel, %ldHz\n",vi->channels,vi->rate);
    fprintf(stderr,"\nDecoded length: %ld samples\n",
            (long)ov_pcm_total(&vf,-1));
    fprintf(stderr,"Encoded by: %s\n\n",ov_comment(&vf,-1)->vendor);
    sample_format.rate = vi->rate;
    sample_format.channels = vi->channels;
    sample_format.bits = 8 * 2;
    sample_format.byte_format = AO_FMT_LITTLE;
  }

  /* ad is where we send audio data */
  ao_device *ad = ao_open_live(driver_id, &sample_format, options);
  
/* Here's the read loop: */

  while(!eof){
    long ret=ov_read(&vf,pcmout,sizeof(pcmout),0,2,1,&current_section);
    if (ret == 0) {
      /* EOF */
      eof=1;
    } else if (ret < 0) {
      /* error in the stream.  Not a problem, just reporting it in
	 case we (the app) cares.  In this case, we don't. */
    } else {
      /* we don't bother dealing with sample rate changes, etc, but
	 you'll have to*/
       /* fwrite(pcmout,1,ret,stdout); */
       ao_play(ad, pcmout, ret);
    }
  }

/* The code is reading blocks of data using ov_read(). Based on the
   value returned, we know if we're at the end of the file or have
   invalid data. If we have valid data, we write it to the pcm
   output. */

/* Now that we've finished playing, we can pack up and go home. It's
   important to call ov_clear() when we're finished. */

  ov_clear(&vf);
  ao_close(ad);  /* close the audio devide */
    
  fprintf(stderr,"Done.\n");
  return(0);
}

