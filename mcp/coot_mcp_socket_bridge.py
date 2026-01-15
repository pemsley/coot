import socket
import json
import struct
import os
import itertools
from pathlib import Path
from mcp.server.fastmcp import FastMCP

# CONFIGURATION
COOT_HOST = "localhost"
COOT_PORT = 9090

# Initialize MCP Server
mcp    = FastMCP("Coot Socket Bridge")
server = FastMCP("Knowing about stuff")

_request_id = itertools.count(0)

def send_coot_rpc(method: str, params: dict) -> dict:
    """
    Connects to Coot via socket, sends a framed JSON-RPC request,
    and awaits a framed response.
    """
    payload = {
        "jsonrpc": "2.0",
        "method": method,
        "params": params or {},
        "id": next(_request_id)
    }

    json_str = json.dumps(payload)
    encoded_msg = json_str.encode('utf-8')

    try:
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
            sock.connect((COOT_HOST, COOT_PORT))

            # --- FRAMING LOGIC (Adjust this if Coot uses \n delimiters) ---

            # OPTION A: 4-Byte Length Prefix (Standard "Framed" approach)
            # We pack the length of the message as a 4-byte big-endian integer
            header = struct.pack('>I', len(encoded_msg))
            sock.sendall(header + encoded_msg)

            # OPTION B: Newline Delimited (Uncomment if Coot uses this)
            # sock.sendall(encoded_msg + b'\n')

            # --- RECEIVING RESPONSE ---

            # 1. Read the length prefix first (4 bytes)
            # If using Option B, you would instead read_until(b'\n')
            header_data = _recv_exact(sock, 4)
            if not header_data:
                return {"error": "No response header from Coot"}

            msg_length = struct.unpack('>I', header_data)[0]

            # 2. Read the exact number of bytes for the body
            response_data = _recv_exact(sock, msg_length)

            if response_data:
                r = response_data.decode('utf-8')
                return json.loads(r)
            else:
                return {}

    except ConnectionRefusedError:
        return {"error": f"Could not connect to Coot at {COOT_HOST}:{COOT_PORT}"}
    except Exception as e:
        return {"error": str(e)}

def _recv_exact(sock, n):
    """Helper to ensure we read exactly n bytes from the stream"""
    data = b''
    while len(data) < n:
        packet = sock.recv(n - len(data))
        if not packet:
            return None
        data += packet
    return data

# --- MCP TOOLS ---

@server.resource("coot://skill")
async def read_resource():
    dirname = Path(os.path.dirname(__file__))
    skill_path = dirname / "docs" / "skills" / "best-practices" / "SKILL.md"
    return skill_path.read_text()

@mcp.tool()
def run_python(code: str) -> str:
    """
    Execute single liness of Python code inside Coot.
    Returns the result of the expression (which might be None)
    """
    # Note: We pass 'code' directly. Ensure Coot's python.exe expects
    # a single string argument or a dict like {"code": ...}
    response = send_coot_rpc("python.exec", {"code": code})

    if "error" in response:
        return f"RPC Error: {response['error']}"

    # Check for internal execution errors from Coot
    if "result" in response:
        return str(response["result"])

    return "No result returned."

@mcp.tool()
def run_python_multiline(code: str) -> str:
    """
    Execute a multi-line Python code block inside Coot.
    Returns None
    Typically, this tool would be used for multi-line code that does not need a return value 
    (but can use print() to display output via stdout) or to define a function that wraps a 
    block of code that should return a value. That function is then called by run_python() 
    so that a value is returned.
    """
    # Note: We pass 'code' directly. Ensure Coot's python.exe expects
    # a single string argument or a dict like {"code": ...}
    response = send_coot_rpc("python.exec_multiline", {"code": code})

    if "error" in response:
        return f"RPC Error: {response['error']}"

    # Check for internal execution errors from Coot
    if "result" in response:
        return str(response["result"])

    return "No result returned."

def get_start_text():
    return '''The Coot API has over 2000 functions.

Use coot_ping() to verify Coot is responsive after:
  - A break in the conversation (no Coot commands for several minutes)
  - Coot has been recompiled or restarted
  - You receive an error suggesting Coot might not be responding
Do NOT call coot_ping() before every command - it's only needed to check if Coot has crashed or been restarted, not for routine operation.

Use search_coot_functions(pattern) to find specific ones (where pattern is a (potentially) multi-words pattern (using space-separated fields). 

Use list_available_tools_in_block(block_index) where block_index varies from 0 to 4 (inclusive) to get each of the api documentation blocks.

run_python() only returns values if the code is a single line.

If you need a return value from a block of code then define a wrapper function in one call to run_python_multiline() (that will return None) and
then run that function in the next call using run_python().

You can try to print values, because you have access to the standard output (using the "stdout" key).

You must call coot.set_refinement_immediate_replacement() before running refinment functions (once is enough) - that should make the refinement synchronous.

For more saftey, call `coot.accept_moving_atoms_py()` after refinement - it should do nothing - but it might be ueful if there is a refinement synchronous bug.

If a model-building tool moves the atoms in a way that you later deem "worse than before" you can call coot.apply_undo() to restore the previous model.

Never try to code that writes to disk - instead, write code that returns a string.

Note that coot.active_atom_spec_py() is a useful function to determine the "selected" residue (i.e. the "active" residue that will be acted on by the tools in the interface).

The coot module is already imported, the coot_utils module will need to be imported first if you want to use a function in that module.

Do not use matplotlib for graphs, instead use pygal.
'''

@mcp.tool()
def list_available_tools() -> str:
    """Ask Coot what tools it supports"""
    # response = send_coot_rpc("mcp.list_tools")
    # return str(response.get("result", []))
    return get_start_text()

@mcp.tool()
def list_available_tools_in_block(block_index: int) -> str:
    """Ask Coot what tools it supports. There are about 5 tool blocks, each returning the descriptions of about 400 functions"""
    response = send_coot_rpc("mcp.list_tools_block", {"block_index": block_index})
    return str(response.get("result", []))

@mcp.tool()
def coot_info() -> str:
    """Returns: usage string'"""
    return get_start_text()

@mcp.tool()
def search_coot_functions(pattern: str) -> str:
    """Search for Coot functions by pattern

    Use this to find specific functions. Common searches:
    - 'correlation' - for map-to-model fit functions
    - 'ramachandran' - for backbone validation
    - 'rotamer' - for side-chain validation

    Returns matching functions with docs (max 40 results)"""
    response = send_coot_rpc("mcp.search", {"pattern": pattern})
    return str(response.get("result", []))

@mcp.tool()
def list_coot_categories() -> str:
    """Returns: ['load', 'read', 'display', 'refinement', 'validation', 'ligand', 'util']"""
    return  str(['load', 'read', 'display', 'refinement', 'validation', 'ligand', 'util'])

@mcp.tool()
def get_functions_in_category(category: str) -> str:
    """Returns functions in that category (maybe 50-200 per category)"""
    return search_coot_functions(category)

@mcp.tool()
def coot_ping() -> str:
    """
    Quick health check - verifies Coot is responsive by having it compute 2+2.
    Returns '4' if Coot is alive and responding correctly.
    """
    response = send_coot_rpc("python.exec", {"code": "2 + 2"})

    if "error" in response:
        return f"Coot not responding: {response['error']}"

    if "result" in response:
        result = response["result"]
        # Check if we got the expected value
        if str(result.get("value")) == "4":
            return "4"
        else:
            return f"Unexpected result from Coot: {result}"

    return "No result from Coot"


if __name__ == "__main__":

    if True:
        mcp.run()

    else:
        # Test connection
        print("Testing mcp.list_tools...")
        result = send_coot_rpc("mcp.list_tools")
        print(json.dumps(result, indent=2))

        print("\nTesting python.exec...")
        result = send_coot_rpc("python.exec", {"code": "1 + 1"})
        print(json.dumps(result, indent=2))


