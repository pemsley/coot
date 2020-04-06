
class keyboard_key_t {
public:
   int gdk_key;
   bool ctrl_is_pressed;
   keyboard_key_t(int g) {
      gdk_key = g;
      ctrl_is_pressed = false;
   }
   keyboard_key_t(int g, bool state) {
      gdk_key = g;
      ctrl_is_pressed = state;
   }
   bool operator<(const keyboard_key_t &other) const {
      if (other.gdk_key < gdk_key) {
         return true;
      } else {
         if (other.gdk_key == gdk_key) {
            if (other.ctrl_is_pressed) {
               if (ctrl_is_pressed) {
                  return false;
               } else {
                  return true;
               }
            } else {
               if (ctrl_is_pressed) {
                  return false; // hmm.
               } else {
                  return false;
               }
            }
         } else {
            return false;
         }
      }
   }
};

class key_bindings_t {

public:
  enum binding_type { NONE, SCHEME, PYTHON, BUILT_IN };
  binding_type type;
  std::string scripting_function_text;
  std::string description;
  void (*func)();
  key_bindings_t() { type = NONE; }
  // external
  key_bindings_t(binding_type bt, const std::string &fn, const std::string &description_in) :
    type(bt), scripting_function_text(fn), description(description_in) {}
  // internal
  key_bindings_t(void (*func_in) (), const std::string &description_in) {
    description = description_in;
    type = BUILT_IN;
    func = func_in;
  }

};
