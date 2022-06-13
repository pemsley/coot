
#include "utils/coot-utils.hh"
#include "startup-utils.hh"
extern "C" int git_revision_count(); // 20220612-PE doesn this need to be extern now?

std::string
make_main_window_title() {

   std::string version_string = VERSION;
   std::string main_title = "Coot " + version_string;

#ifdef MAKE_ENHANCED_LIGAND_TOOLS
   // main_title += " EL";
#endif

#ifdef COOT_MAIN_TITLE_EXTRA
   main_title += COOT_MAIN_TITLE_EXTRA;
#else

   // if this is a pre-release, stick in the revision number too
   if (version_string.find("-pre") != std::string::npos) {
      main_title += " (revision count ";
      main_title += coot::util::int_to_string(git_revision_count());
      main_title += ")";
   }
#endif

#ifdef WINDOWS_MINGW
   main_title = "Win" + main_title;
#endif

   return main_title;
}
