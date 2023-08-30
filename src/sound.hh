
#ifndef SOUND_HH
#define SOUND_HH

#include <string>

int test_sound(int argc, char **argv);

void play_sound_file(const std::string &file_name);

void play_sound(const std::string &type);

#endif // SOUND_HH

