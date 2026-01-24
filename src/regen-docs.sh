
echo doxygen...
doxygen coot-dox.cfg > doxygen.log 2>&1
echo conveting from xml to feature docstrings
python3 xml-to-swig-docstrings.py --output       coot_c_interface_docs.i doxygen/xml/c-interface_8h.xml
python3 xml-to-swig-docstrings.py --output      coot_cc_interface_docs.i doxygen/xml/cc-interface_8hh.xml
python3 xml-to-swig-docstrings.py --output     coot_rsr_functions_docs.i doxygen/xml/rsr-functions_8hh.xml
python3 xml-to-swig-docstrings.py --output     coot_read_molecule_docs.i doxygen/xml/read-molecule_8hh.xml
python3 xml-to-swig-docstrings.py --output coot_network_functions_docs.i doxygen/xml/network-get_8hh.xml
python3 xml-to-swig-docstrings.py --output coot_cc_interface_molecular_representation.i  doxygen/xml/cc-interface-molecular-representation_8hh.xml
python3 xml-to-swig-docstrings.py --output coot_cc_interface_user_defined_atom_colours.i doxygen/xml/cc-interface-user-defined-atom-colours_8hh.xml
