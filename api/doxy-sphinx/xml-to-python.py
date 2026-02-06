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
mytree = ET.parse('xml/namespacegemmi.xml')
myroot = mytree.getroot()

def convert_type(tt: str) -> str:
    if tt == "const ": tt = "str" # needed for a coot::colour_t
    if tt == "const std::string &": tt = "str"
    if tt == "std::string": tt = "str"
    if tt == "unsigned int": tt = "int"
    if tt == "double": tt = "float"
    if tt == 'std::vector< ': tt = "list"
    if tt == 'const std::vector< ': tt = "list"
    # if tt == 'std::pair< int, unsigned int >': tt = "tuple"
    if tt == 'const std::vector< float > &': tt = "list"
    if tt == 'std::vector< float > &': tt = "list"
    if tt == "const std::vector< std::string > &": tt = "list"
    if tt == "const std::vector< std::string > &": tt = "list"
    if tt == 'const std::map< unsigned int, std::array< float, 3 > > &': tt = "dict"
    if tt == 'const std::map< unsigned int, std::array< float, 4 > > &': tt = "dict"
    if tt == 'const std::vector< std::pair< std::string, unsigned int > > &': tt = "list"
    if tt == 'const std::vector< std::pair< std::string, unsigned int > > &': tt = "list"
    if tt == 'const std::vector< std::pair< bool, mmdb::Residue * > > &, links: const std::vector< mmdb::Link > &': tt = "list"
    if tt == 'std::vector<': tt = "list"
    if tt == 'const coot::residue_spec_t &': tt = 'str'
    if tt == 'const coot::colour_t &': tt = 'list'   #use the proper type later
    if tt == 'const coot::atom_spec_t &': tt = 'str'
    if tt == 'coot::residue_spec_t &': tt = 'list'
    if tt == 'std::vector< coot::api::moved_atom_t > &': tt = 'list'
    if tt == 'const std::vector< coot::api::moved_residue_t > &': tt = 'list'

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
    # print("--- make_return_type() dict is ", function)
    return_type = ""
    rt = ""
    t = str(function['type'])
    if t == 'int':   rt = "int"
    if t == 'bool':  rt = "bool"
    if t == 'void':  rt = "None"
    if t == 'float': rt = "float"
    if t == 'std::string':  rt = "str"
    if t == 'unsigned int': rt = "int"
    if t == 'std::pair< int, unsigned int >': rt = "tuple"
    if rt:
        return_type = " -> " + rt
    return rt, return_type


def make_python_script(functions: list, output_file_name: str) -> None:

    f = open(output_file_name, "w")
    f.write("class molecules_container_t:\n")
    for function in functions:
        func_name = ""
        try:
            func_name = function['name']
        except AttributeError as e:
            pass
        print("\n--- make_python_script(): Handling function: ", func_name)

        parens = ""
        def_ = ""
        if function['kind'] == "function":
            parens = make_paren_string(function)
            return_type, return_type_with_arrow = make_return_type(function)
            # print("debug args: function ", function, "made args:", parens)
            def_ = "def "
            s = f"    {def_}{function['name']}{parens}{return_type_with_arrow}:\n"
            f.write(s)

            done_brief = False
            done_detailed = False
            f.write('        """ ')
            try:
                d = function["briefdescription"]
                # print(f"debug:: briefdescription:{d}:")
                f.write(d)
                # f.write('\n')
                done_brief = True
            except KeyError as e:
                # print("No briefdescription")
                pass
            except TypeError as e:
                pass
            # maybe use a for loop for this and the above ["briefdescription", "detaileddescription"]
            try:
                d = function["detaileddescription"]
                # print(f"debug:: detaileddescription:{d}:")
                # f.write('        """ ')
                f.write(d)
                # f.write('\n')
                done_detailed = True
            except KeyError as e:
                # print("No detaileddescription")
                pass
            except TypeError as e:
                pass
            if not done_brief:
                if not done_detailed:
                    # we need some doc string for the sphinx to pick up the function
                    f.write('        Sphinx-Doc-Placeholder')
            f.write('"""')
            f.write('\n')
            if not return_type:               f.write("        pass\n")
            if return_type == "int":          f.write("        return 0\n")
            if return_type == "float":        f.write("        return 0.0\n")
            if return_type == "str":          f.write("        return 'a-string'\n")
            if return_type == "bool":         f.write("        return True\n")

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
            keep_going = True
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
                        if ch.text == "molecules_container_t::~molecules_container_t":
                            keep_going = False
                            break
                            # next memberdef
                        if ch.text == "molecules_container_t::molecules_container_t":
                            keep_going = False
                            break
                    if ch.tag == "param":
                        t = ch.find("type")
                        print("   param type:", t)
                        tt = t.text
                        print("    tt: '" + str(tt) + "'")
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
                    parts = ''
                    if ch.tag == "detaileddescription":
                        for idx,c in enumerate(ch):
                            print("     detaileddescription: item", idx, "is:", c, c.text)
                            if c.tag == "para":
                                # a para can have a parameter list and no text
                                #   descr = c.text
                                #   if descr:
                                #       print("      detailed descr", descr)
                                #       a_function["detaileddescription"] = descr
                                if c.text:
                                    if parts:
                                        parts += "\n\n        "
                                        parts += c.text
                                    else:
                                        parts = c.text
                                # print(dir(c))
                                for jj, chunk in enumerate(c):
                                    print('      chunk', jj, ":", chunk, "len:", len(chunk))
                                    if chunk.tag == "parameterlist":
                                        print("        handle parameterlist here")

                                        for pk,pchunk in enumerate(chunk):
                                               print("pk:", pk, 'pchunk', pchunk)
                                               for ck,cchunk in enumerate(pchunk):
                                                   print("ck:", ck, 'cchunk', cchunk)

                                                   if cchunk.tag == "parameternamelist":
                                                      parts += "\n"
                                                      for kk,kchunk in enumerate(cchunk):
                                                           if kchunk.tag == "parametername":
                                                              print("kchunk", kk, "text:", kchunk.text)
                                                              if kchunk.text:
                                                                 parts += "\n       " + " :param " + kchunk.text + ": "

                                                   if cchunk.tag == "parameterdescription":
                                                       print('len(cchunk)', len(cchunk))
                                                       for kk,kchunk in enumerate(cchunk):
                                                           if kchunk.tag == "para":
                                                               if kchunk.text:
                                                                   print("kchunk", kk, "text:", kchunk.text)
                                                                   parts += " " + kchunk.text
                                                                   for kk,pchunk in enumerate(kchunk):
                                                                      if pchunk.tag == "computeroutput":
                                                                         if pchunk.text:
                                                                             if pchunk.tail:
                                                                                 print("        computeroutput:", pchunk.text)
                                                                                 parts += "`"
                                                                                 parts += pchunk.text
                                                                                 parts += "`"
                                                                                 parts += pchunk.tail

                                    if chunk.tag == "computeroutput":
                                        print("        computeroutput:", chunk.text)
                                        if chunk.text:
                                            if chunk.tail:
                                                parts += "`"
                                                parts += chunk.text
                                                parts += "`"
                                                parts += chunk.tail
                                    if chunk.tag == "simplesect":
                                        for kk,kchunk in enumerate(chunk):
                                            # print("kk:", kk, kchunk)
                                            if kchunk.tag == "para":
                                                print("kchunk", kk, "text:", kchunk.text)
                                                #parts += "Return" + " " + kchunk.text
                                                if kchunk.text:
                                                    parts += "\n\n       " + " :return: " + kchunk.text
                                                    for kk,pchunk in enumerate(kchunk):
                                                        if pchunk.tag == "computeroutput":
                                                           if pchunk.text:
                                                               if pchunk.tail:
                                                                   print("        computeroutput:", pchunk.text)
                                                                   parts += "`"
                                                                   parts += pchunk.text
                                                                   parts += "`"
                                                                   parts += pchunk.tail


                    if parts:
                        a_function["detaileddescription"] = parts
                if a_function:
                    if keep_going:
                        functions.append(a_function)

            except AttributeError as e:
                print(e)

fn = "gemmi_funcs.py"
make_python_script(functions, fn)
