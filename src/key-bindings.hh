

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
