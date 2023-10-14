#include <emscripten/bind.h>
#include <iostream>
using namespace emscripten;

class LhasaDemo {
    int data;

    public:

    std::string demo() {
        return "Hello from demo!";
    }

    LhasaDemo() {

    }

};

float lerp(float a, float b, float t) {
    return (1 - t) * a + t * b;
}

EMSCRIPTEN_BINDINGS(my_module) {
  function("lerp", &lerp);
  class_<LhasaDemo>("LhasaDemo")
    .constructor<>()
    .function("demo", &LhasaDemo::demo);
}

// int main() {

// }