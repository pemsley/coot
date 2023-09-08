
#include <AudioToolbox/AudioToolbox.h>

   int main (int argc, const char * argv[]) {
       SystemSoundID mySSID;
       CFURLRef myURLRef = CFURLCreateWithFileSystemPath (
           kCFAllocatorDefault,
           CFSTR ("../../ComedyHorns.aif"),
           kCFURLPOSIXPathStyle,
           FALSE
       );

       AudioServicesCreateSystemSoundID (myURLRef, &mySSID);
       AudioServicesPlaySystemSound (mySSID);

       CFRunLoopRun ();

       return 0;
   }

