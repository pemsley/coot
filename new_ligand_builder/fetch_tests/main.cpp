#include <iostream>
#include "../../src/cc-interface-network.hh"
#include "../../src/coot-setup-python.hh"

int main(int argc, char** argv) {
    setup_python_basic(argc, argv);
    setup_python_coot_module();
    safe_python_command("import coot_utils");
    std::string input;
    std::cout<<"Drug to fetch: \n";
    std::cin>>input;
    auto res = get_drug_via_wikipedia_and_drugbank_py(input);
    std::cout<<"\nResult: \n"<<res<<'\n';
}