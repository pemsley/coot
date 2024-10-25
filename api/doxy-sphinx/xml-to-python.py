# This script is an attempt to use the XML output of doxygen to
# create a python file that contains stubs of functions
# which can then be used to generate python documentation
# with Sphinx.
#
# It does some of that work, but to get it working properly
# will be about two weeks of work - or so.
# And spending that time to save maybe 3 hours of boring
# editing is not worth it.
#
# It is better to just create a python stub script now
# and modify it by tracking changes made to the nanobinds file.

import xml.etree.ElementTree as ET
mytree = ET.parse('xml/classmolecules__container__t.xml')
myroot = mytree.getroot()

def convert_type(tt: str) -> str:
    if tt == "const std::string &": tt = "str"
    if tt == "std::string": tt = "str"
    if tt == "unsigned int": tt = "int"
    if tt == "double": tt = "float"
    if tt == 'std::vector< ': tt = "list"
    if tt == 'const std::vector< ': tt = "list"
    if tt == "const std::vector< std::string > &": tt = "list"
    if tt == "const std::vector< std::string > &": tt = "list"
    if tt == 'const std::map< unsigned int, std::array< float, 3 > > &': tt = "dict"
    if tt == 'const std::vector< std::pair< std::string, unsigned int > > &': tt = "list"
    if tt == 'const std::vector< std::pair< std::string, unsigned int > > &': tt = "list"
    if tt == 'const std::vector< std::pair< bool, mmdb::Residue * > > &, links: const std::vector< mmdb::Link > &': tt = "list"
    if tt == 'std::vector<': tt == list
    return tt

def make_paren_string(function: dict) -> str:
    if function['kind'] == "function":
        r = ""
        try:
            for idx,param in enumerate(function['params']):
                if idx > 0:
                    r += ", "
                r += param['declname']
                r += ": "
                t = convert_type(param['type'])
                r += t
            r = "(self, " + r + ")"
            return r
        except KeyError as e:
            return "(self)"
    else:
        return ""

def make_return_type(function: dict) -> str:
    return_type = " -> float"
    return_type = ""
    return return_type

def make_python_script(functions: list) -> None:

    f = open("chapi-functions.py", "w")
    f.write("class molecules_container_t:\n")
    for function in functions:
        print("\n--- make_python_script(): Handling function: ", function, ":")

        parens = ""
        def_ = ""
        if function['kind'] == "function":
            parens = make_paren_string(function)
            return_type = make_return_type(function)
            print("debug args: function ", function, "made args:", parens)
            def_ = "def "
            s = f"    {def_}{function['name']}{parens}{return_type}:\n"
            f.write(s)

            done_brief = False
            done_detailed = False
            try:
                d = function["briefdescription"]
                print(f"debug:: briefdescription:{d}:")
                f.write('        """ ')
                f.write(d)
                f.write('"""')
                f.write('\n')
                done_brief = True
            except KeyError as e:
                # print("No briefdescription")
                pass
            except TypeError as e:
                pass
            # maybe use a for loop for this and the above ["briefdescription", "detaileddescription"]
            try:
                d = function["detaileddescription"]
                print(f"debug:: detaileddescription:{d}:")
                f.write('        """ ')
                f.write(d)
                f.write('"""')
                f.write('\n')
                done_detailed = True
            except KeyError as e:
                # print("No detaileddescription")
                pass
            except TypeError as e:
                pass
            if not done_brief:
                if not done_detailed:
                    # we need some doc string for the sphinx to pick up the function
                    f.write('        """Sphinx-Doc-Placeholder"""\n')
            f.write("        pass\n")

        f.write("\n")
    f.close()

functions = []

for x in myroot.iter('sectiondef'):

    print("#### x:", x)
    # print(dir(x))
    h = x.find('header')
    print(h)
    if h is not None:
        ht = h.text
        print("####### header text: ", ht)
    for child in x:
        if child.tag == "memberdef":
            try:
                name = "--unset--"
                a_function = {}
                kind = child.attrib['kind']
                print("   child kind:", kind)
                a_function['kind'] = kind
                if kind == "function":
                    print("\n -----  Handling function")
                    # if function.name ==  "molecules_container_t": continue
                    # if function.name == "~molecules_container_t": continue
                for ii,ch in enumerate(child):
                    print("      ch.tag ", ii, ch.tag, ch.text, ":")
                    if ch.tag == "definition":
                        print('============ definition', ch.tag)
                        if ch.text == "molecules_container_t::~molecules_container_t":
                            print('breaking out')
                            break
                            # next memberdef
                        if ch.text == "molecules_container_t::molecules_container_t":
                            print('breaking out')
                            break
                    if ch.tag == "param":
                        print("   found a param!")
                        t = ch.find("type")
                        print("   debug:: param type:", t)
                        tt = t.text
                        print("    tt: '" + tt + "'")
                        dn = ch.find("declname")
                        dn = dn.text
                        print("    dn", dn)
                        if tt:
                            try:
                                a_function['params'].append({"type": tt, "declname": dn})
                            except KeyError as e:
                                a_function['params'] = [{"type": tt, "declname": dn}]
                    if ch.tag == "name":
                        name = ch.text
                        a_function["name"] = name
                    if ch.tag == "type":
                        a_function["type"] = ch.text
                    if ch.tag == "briefdescription":
                        for c in ch:
                          if c.tag == "para":
                              brief_descr = c.text
                              a_function["briefdescription"] = brief_descr
                    if ch.tag == "detaileddescription":
                        for idx,c in enumerate(ch):
                          print("   detaileddescription: item", idx, "is:", c)
                          if c.tag == "para":
                              # a para can have a parameter list and no text
                              descr = c.text
                              if descr:
                                  print("Here with descr", descr, " for name ", name)
                                  a_function["detaileddescription"] = descr
                if a_function:
                    functions.append(a_function)

            except AttributeError as e:
                print(e)

make_python_script(functions)
