# This script is an attempt to use the XML output of doxygen to
# create a python file that contains stubs of functions
# which can then be used to generate python documenation
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
            r = "(" + r + ")"
            return r
        except KeyError as e:
            return "()"
    else:
        return ""


def make_python_script(functions: list) -> None:

    f = open("chapi-functions.py", "w")
    f.write("def molecules_container_t():\n\n")
    for function in functions:
        print("Handling function: ", function, ":")
        try:
            d = function["briefdescription"]
            print(f"debug:: briefdescription:{d}:")
            f.write("    # ")
            f.write(d)
            f.write('\n')
        except KeyError as e:
            # print("No briefdescription")
            pass
        except TypeError as e:
            pass
        # maybe use a for loop for this and the above ["briefdescription", "detaileddescription"]
        try:
            d = function["detaileddescription"]
            print(f"debug:: detaileddescription:{d}:")
            f.write("    # ")
            f.write(d)
            f.write('\n')
        except KeyError as e:
            # print("No detaileddescription")
            pass
        except TypeError as e:
            pass

        parens = ""
        def_ = ""
        if function['kind'] == "function":
            parens = make_paren_string(function)
            def_ = "def "

        s = f"    {def_}{function['name']}{parens}\n"
        f.write(s)
        if function['kind'] == "function":
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
                    print("   Handling function")
                for ch in child:
                    print("      ch.tag ", ch.tag, ch.text, ":")
                    if ch.tag == "param":
                        print("   found a param!")
                        t = ch.find("type")
                        print("   debug:: param type:", t)
                        tt = t.text
                        print("    tt", tt)
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
                        for c in ch:
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
