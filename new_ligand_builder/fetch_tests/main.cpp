#include <iostream>
#include <string>
#include "../../src/cc-interface-network.hh"
#include "../../src/coot-setup-python.hh"


void repl() {
    std::cout<<"Coot pyrepl:\n";
    while(true) {
        std::string input;
        std::getline(std::cin,input);
        PyRun_SimpleString(input.c_str());
    }
}

int main(int argc, char** argv) {
    setup_python_basic(argc, argv);
    setup_python_coot_module();
    PyRun_SimpleString("import coot_utils");
    repl();
    return 0;

    std::string input;
    std::cout<<"Drug to fetch: \n";
    std::cin>>input;
    auto res = get_drug_via_wikipedia_and_drugbank_py(input);
    std::cout<<"\nResult: \n"<<res<<'\n';
}