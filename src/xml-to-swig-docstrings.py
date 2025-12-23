#!/usr/bin/env python3
"""
Improved version - extracts sections properly without duplication
"""

import xml.etree.ElementTree as ET
import sys
import argparse

def convert_type(tt: str) -> str:
    """Convert C++ types to Python type hints for documentation"""
    type_map = {
        "const ": "str",
        "const std::string &": "str",
        "std::string": "str",
        "unsigned int": "int",
        "double": "float",
        'PyObject *': 'object',
        'char const *': 'str',
        'short': 'int',
        'bool': 'bool',
        'int': 'int',
        'float': 'float',
        'void': 'None',
    }
    
    # Try exact match first
    if tt in type_map:
        return type_map[tt]
    
    # Pattern matching for complex types
    if 'std::vector<' in tt or 'const std::vector<' in tt:
        return "list"
    if 'std::map<' in tt:
        return "dict"
    if 'coot::residue_spec_t' in tt:
        return 'residue_spec'
    if 'coot::colour_t' in tt:
        return 'colour'
    if 'coot::atom_spec_t' in tt:
        return 'atom_spec'
    
    return tt

def extract_all_text(elem):
    """
    Recursively extract ALL text from an element, including nested elements.
    Returns plain text with newlines preserved where appropriate.
    """
    parts = []
    
    if elem.text:
        parts.append(elem.text.strip())
    
    for child in elem:
        # Recursively get text from children
        child_text = extract_all_text(child)
        if child_text:
            parts.append(child_text)
        
        # Get tail text (after closing tag)
        if child.tail:
            tail = child.tail.strip()
            if tail:
                parts.append(tail)
    
    return " ".join(parts)

def format_as_list(itemizedlist_elem, indent=0):
    """Format an itemizedlist element as a bulleted list"""
    items = []
    indent_str = "  " * indent
    
    for listitem in itemizedlist_elem.findall("listitem"):
        # Get all text from this list item
        item_text = extract_all_text(listitem)
        if item_text:
            # Check for nested lists
            nested_lists = listitem.findall(".//itemizedlist")
            if nested_lists:
                # Add the main item
                items.append(f"{indent_str}- {item_text.split('  -')[0].strip()}")
                # Add nested items with more indentation
                for nested in nested_lists:
                    nested_items = format_as_list(nested, indent + 1)
                    items.extend(nested_items)
            else:
                items.append(f"{indent_str}- {item_text}")
    
    return items

def parse_detailed_description(detaileddesc_elem):
    """
    Parse the detaileddescription element and return structured content.
    Returns dict with 'intro' and 'sections' keys.
    """
    result = {
        'intro': '',
        'sections': []
    }
    
    # Get introductory text (direct para children before any sections)
    intro_parts = []
    for child in detaileddesc_elem:
        if child.tag == "para":
            # Extract text, handling computeroutput
            text = extract_all_text(child)
            if text:
                intro_parts.append(text)
        elif child.tag in ["sect1", "sect2"]:
            # Stop when we hit sections
            break
    
    if intro_parts:
        result['intro'] = "\n\n".join(intro_parts)
    
    # Now process all sect1/sect2 sections
    for sect in detaileddesc_elem.findall(".//sect2"):
        section = {'title': '', 'content': []}
        
        # Get title
        title_elem = sect.find("title")
        if title_elem is not None and title_elem.text:
            section['title'] = title_elem.text.strip()
        
        # Process para elements in this section
        for para in sect.findall("para"):
            # Check for special structures
            paramlist = para.find("parameterlist")
            if paramlist is not None:
                # Skip - we'll handle parameters separately
                continue
            
            simplesect = para.find("simplesect")
            if simplesect is not None and simplesect.attrib.get('kind') == 'return':
                # Get return description
                ret_text = extract_all_text(simplesect)
                if ret_text:
                    section['content'].append(ret_text)
            
            # Check for itemizedlist
            itemlist = para.find("itemizedlist")
            if itemlist is not None:
                list_items = format_as_list(itemlist)
                if list_items:
                    section['content'].extend(list_items)
            else:
                # Regular paragraph
                text = extract_all_text(para)
                if text:
                    section['content'].append(text)
        
        if section['title'] or section['content']:
            result['sections'].append(section)
    
    return result

def make_swig_docstring(function: dict) -> str:
    """Generate docstring content for SWIG %feature directive"""
    
    doc_parts = []
    
    # Brief description
    if 'briefdescription' in function and function['briefdescription']:
        doc_parts.append(function['briefdescription'].strip())
    
    # Introductory text from detailed description
    if 'intro' in function and function['intro']:
        if doc_parts:
            doc_parts.append("")
        doc_parts.append(function['intro'])
    
    # Process structured sections
    if 'sections' in function and function['sections']:
        for section in function['sections']:
            title = section.get('title', '')
            content = section.get('content', [])
            
            if not (title or content):
                continue
            
            doc_parts.append("")
            if title:
                doc_parts.append(title)
                # Add underline for major sections
                if title in ['Parameters', 'Returns', 'Return value', 'Notes', 
                             'Examples', 'Reference counting']:
                    doc_parts.append("-" * len(title))
            
            # Add content
            if content:
                for line in content:
                    doc_parts.append(line)
    
    # Add parameter info if we have it
    if 'params' in function and function['params']:
        # Check if we already have a Parameters section
        has_param_section = any(
            s.get('title') == 'Parameters' 
            for s in function.get('sections', [])
        )
        
        if not has_param_section:
            doc_parts.append("")
            doc_parts.append("Parameters")
            doc_parts.append("----------")
            for param in function['params']:
                param_type = convert_type(param['type'])
                param_name = param['declname']
                doc_parts.append(f"{param_name} : {param_type}")
                if 'description' in param and param['description']:
                    desc_lines = param['description'].strip().split('\n')
                    for line in desc_lines:
                        doc_parts.append(f"    {line}")
    
    return "\n".join(doc_parts)

def escape_docstring(docstring: str) -> str:
    """Escape special characters for SWIG docstring"""
    docstring = docstring.replace('\\', '\\\\')
    docstring = docstring.replace('"', '\\"')
    return docstring

def make_swig_interface(functions: list, output_file: str = "coot_docstrings.i") -> None:
    """Generate SWIG interface file with docstring features"""
    
    with open(output_file, "w") as f:
        f.write("/* Auto-generated SWIG docstrings from Doxygen XML */\n")
        f.write("/* Generated by xml-to-swig-docstrings.py */\n\n")
        
        function_count = 0
        
        for function in functions:
            if function.get('kind') != "function":
                continue
            
            func_name = function.get('name', '')
            if not func_name or func_name.startswith('~'):
                continue
            
            docstring = make_swig_docstring(function)
            
            if docstring.strip():
                escaped_docstring = escape_docstring(docstring)
                f.write(f'%feature("docstring") {func_name} "\n')
                f.write(escaped_docstring)
                f.write('\n";\n\n')
                function_count += 1
        
        print(f"✓ Generated SWIG docstrings for {function_count} functions")
        print(f"✓ Output written to: {output_file}")

def parse_doxygen_xml(xml_file: str) -> list:
    """Parse Doxygen XML and extract function information"""
    
    try:
        mytree = ET.parse(xml_file)
    except FileNotFoundError:
        print(f"Error: Cannot find XML file: {xml_file}")
        sys.exit(1)
    except ET.ParseError as e:
        print(f"Error parsing XML: {e}")
        sys.exit(1)
    
    myroot = mytree.getroot()
    functions = []
    
    for sectiondef in myroot.iter('sectiondef'):
        for memberdef in sectiondef:
            if memberdef.tag != "memberdef":
                continue
            
            if memberdef.attrib.get('kind') != 'function':
                continue
            
            try:
                a_function = {'kind': 'function'}
                
                # Get function name
                name_elem = memberdef.find("name")
                if name_elem is None or not name_elem.text:
                    continue
                a_function['name'] = name_elem.text
                
                # Skip destructors
                definition_elem = memberdef.find("definition")
                if definition_elem is not None and definition_elem.text and "~" in definition_elem.text:
                    continue
                
                # Get return type
                type_elem = memberdef.find("type")
                if type_elem is not None:
                    a_function['type'] = type_elem.text if type_elem.text else ""
                
                # Get parameters
                params = []
                for param in memberdef.findall("param"):
                    t = param.find("type")
                    dn = param.find("declname")
                    params.append({
                        "type": t.text if t is not None and t.text else "",
                        "declname": dn.text if dn is not None and dn.text else ""
                    })
                if params:
                    a_function['params'] = params
                
                # Get brief description
                briefdesc = memberdef.find("briefdescription")
                if briefdesc is not None:
                    for para in briefdesc.findall("para"):
                        text = extract_all_text(para)
                        if text:
                            a_function['briefdescription'] = text
                
                # Get detailed description with structured sections
                detaileddesc = memberdef.find("detaileddescription")
                if detaileddesc is not None:
                    parsed = parse_detailed_description(detaileddesc)
                    if parsed['intro']:
                        a_function['intro'] = parsed['intro']
                    if parsed['sections']:
                        a_function['sections'] = parsed['sections']
                    
                    # Also extract parameter descriptions to add to params
                    for para in detaileddesc.findall(".//para"):
                        paramlist = para.find("parameterlist")
                        if paramlist is not None:
                            for paramitem in paramlist.findall("parameteritem"):
                                param_name_elem = paramitem.find(".//parametername")
                                param_desc_elem = paramitem.find(".//parameterdescription")
                                
                                if param_name_elem is not None and param_desc_elem is not None:
                                    param_name = param_name_elem.text
                                    param_desc = extract_all_text(param_desc_elem)
                                    
                                    # Add description to matching parameter
                                    if 'params' in a_function:
                                        for p in a_function['params']:
                                            if p['declname'] == param_name:
                                                p['description'] = param_desc
                
                functions.append(a_function)
            
            except Exception as e:
                func_name = a_function.get('name', 'unknown')
                print(f"Warning: Error processing function {func_name}: {e}")
                import traceback
                traceback.print_exc()
                continue
    
    return functions

def main():
    parser = argparse.ArgumentParser(description='Convert Doxygen XML to SWIG docstrings')
    parser.add_argument('xml_file', help='Doxygen XML file to process')
    parser.add_argument('-o', '--output', default='coot_docstrings.i', help='Output .i file')
    
    args = parser.parse_args()
    
    print(f"Parsing Doxygen XML: {args.xml_file}")
    functions = parse_doxygen_xml(args.xml_file)
    
    print(f"Found {len(functions)} functions")
    
    if functions:
        make_swig_interface(functions, args.output)
    else:
        print("Warning: No functions found in XML file")
        return 1
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
