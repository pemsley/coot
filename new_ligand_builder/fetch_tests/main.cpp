#include <iostream>
#include <string>
#include "../../src/cc-interface-network.hh"
#include "../../src/coot-setup-python.hh"


void repl() {
    std::cout<<"Coot pyrepl:\n";
    PyObject *am = PyImport_AddModule("__main__");
    PyObject* d = PyModule_GetDict(am);
    while(true) {
        std::string input;
        std::cout<<"\n> ";
        std::flush(std::cout);
        std::getline(std::cin,input);
        PyObject* source_code = Py_CompileString(input.c_str(), "adhoc", Py_eval_input);
        if(!source_code) {
            std::cout<<"Null source_code.\n";
            PyErr_Print();
            continue;
        }
        PyObject* func = PyFunction_New(source_code, d);
        PyObject* result = PyObject_CallObject(func, PyTuple_New(0));
        if(!result) {
            std::cout<<"Null result.\n";
            PyErr_Print();
            Py_XDECREF(func);
            Py_XDECREF(source_code);
            continue;
        }
        const char* s = PyUnicode_AsUTF8(result);
        if(!s) {
            std::cout<<"Result is most likely not a string.\n";
        } else {
            std::cout<<"Result: "<<s<<'\n';
        }
        Py_XDECREF(result);
        Py_XDECREF(func);
        Py_XDECREF(source_code);
    }
}

int main(int argc, char** argv) {
    setup_python_basic(argc, argv);
    setup_python_coot_module();
    PyRun_SimpleString("import coot_utils");

    std::string input;
    std::cout<<"Drug to fetch: \n";
    std::getline(std::cin,input);
    auto res = get_drug_via_wikipedia_and_drugbank_py(input);
    std::cout<<"\nResult: \n"<<res<<'\n';
    repl();
}