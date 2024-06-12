#ifndef BASE64_ENCODE_DECODE_HH
#define BASE64_ENCODE_DECODE_HH

#include <ctype.h>
#include <string>

namespace moorhen_base64 {

   static inline bool is_base64(unsigned char c) {
      return (isalnum(c) || (c == '+') || (c == '/'));
   }

   std::string base64_encode(unsigned char const *bytes_to_encode,
                             unsigned int in_len);

   std::string base64_decode(std::string const &encoded_string);
}


#endif // BASE64_ENCODE_DECODE_HH
