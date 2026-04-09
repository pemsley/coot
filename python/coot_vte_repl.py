#!/usr/bin/env python3
#
# coot_vte_repl.py
#
# A lightweight Python REPL helper for the Coot VTE terminal.
# This script runs as a subprocess inside the VTE terminal widget.
# It communicates with the main Coot process over a Unix socket
# to execute Python code on the main thread.
#
# Usage: python3 coot_vte_repl.py <socket_fd>
#

import sys
import os
import json
import code
import signal
import readline

def main():

   if len(sys.argv) < 2:
      print("Usage: coot_vte_repl.py <socket_fd>", file=sys.stderr)
      sys.exit(1)

   socket_fd = int(sys.argv[1])

   # Open the socket as file objects for convenient I/O
   # We duplicate the fd so we can have separate read/write streams
   sock_r = os.fdopen(os.dup(socket_fd), 'r', buffering=1)  # line-buffered
   sock_w = os.fdopen(socket_fd, 'w', buffering=1)           # line-buffered

   def send_and_receive(code_str, multiline=False):
      """Send code to the main Coot process and return the result."""
      request = json.dumps({"code": code_str, "multiline": multiline})
      try:
         sock_w.write(request + "\n")
         sock_w.flush()
      except (BrokenPipeError, OSError):
         print("\n[VTE REPL: Connection to Coot lost]", file=sys.stderr)
         sys.exit(1)

      try:
         response_line = sock_r.readline()
         if not response_line:
            print("\n[VTE REPL: Connection to Coot closed]", file=sys.stderr)
            sys.exit(1)
         return json.loads(response_line)
      except (json.JSONDecodeError, OSError) as e:
         print(f"\n[VTE REPL: Error reading response: {e}]", file=sys.stderr)
         return {"stdout": "", "result": "", "error": str(e)}

   def request_completions(text):
      """Ask Coot for tab-completions of text."""
      request = json.dumps({"complete": text})
      try:
         sock_w.write(request + "\n")
         sock_w.flush()
         response_line = sock_r.readline()
         if not response_line:
            return [], None
         response = json.loads(response_line)
         return response.get("completions", []), response.get("doc")
      except (BrokenPipeError, OSError, json.JSONDecodeError):
         return [], None

   # Set up readline tab-completion
   completion_cache = []
   completion_doc = None

   def completer(text, state):
      nonlocal completion_cache, completion_doc
      if state == 0:
         # New completion request — query Coot
         completion_cache, completion_doc = request_completions(text)
         # If single match with doc, show it below the current line
         if len(completion_cache) == 1 and completion_doc:
            # Save cursor, move to new line, print doc, restore
            buf = readline.get_line_buffer()
            print(f"\n\033[36m{completion_doc}\033[0m")  # cyan
            # Redisplay prompt + buffer so readline isn't confused
            prompt = ps2 if buffer_lines else ps1
            sys.stdout.write(prompt + buf)
            sys.stdout.flush()
      if state < len(completion_cache):
         return completion_cache[state]
      return None

   def display_matches(substitution, matches, longest_match_length):
      """Custom display hook for readline to show multiple matches."""
      print()
      for match in matches:
         print(f"  {match}")
      # Re-display the prompt and current input
      prompt = ps2 if buffer_lines else ps1
      sys.stdout.write(prompt + readline.get_line_buffer())
      sys.stdout.flush()

   readline.set_completer(completer)
   readline.set_completer_delims(' \t\n`~!@#$%^&*()-=+[{]}\\|;:\'",<>/?')
   # macOS uses libedit which needs different binding syntax
   if 'libedit' in readline.__doc__:
      readline.parse_and_bind("bind ^I rl_complete")
   else:
      readline.parse_and_bind("tab: complete")
   readline.set_completion_display_matches_hook(display_matches)

   # Handle Ctrl-C gracefully
   signal.signal(signal.SIGINT, signal.SIG_IGN)

   print("Coot Python Console: type Python commands to execute in Coot.")
   print()

   ps1 = ">>> "
   ps2 = "... "
   buffer_lines = []

   while True:
      try:
         prompt = ps2 if buffer_lines else ps1
         line = input(prompt)
      except EOFError:
         print()
         break
      except KeyboardInterrupt:
         print()
         buffer_lines.clear()
         continue

      if buffer_lines:
         # We are in a multi-line block
         if line.strip() == "":
            # Empty line ends the block
            full_code = "\n".join(buffer_lines)
            buffer_lines.clear()
            response = send_and_receive(full_code, multiline=True)
         else:
            buffer_lines.append(line)
            continue
      else:
         # Check if this line starts a multi-line block
         # (ends with ':', or is an incomplete statement)
         stripped = line.strip()
         if stripped == "":
            continue

         try:
            compiled = code.compile_command(line)
            if compiled is None:
               # Incomplete statement — start multi-line mode
               buffer_lines.append(line)
               continue
         except SyntaxError:
            # Let Coot handle the syntax error for proper reporting
            pass

         # Check for block-starting keywords with trailing colon
         if stripped.endswith(":") and any(
            stripped.startswith(kw) for kw in
            ["if ", "if(", "for ", "for(", "while ", "while(",
             "def ", "class ", "try:", "except", "else:", "elif ",
             "finally:", "with ", "async "]):
            buffer_lines.append(line)
            continue

         # Single-line expression or statement
         # Try as expression first (non-multiline), fall back to statement
         response = send_and_receive(line, multiline=False)

         # If expression eval failed, try as statement
         if response.get("error", "") and not response.get("stdout", ""):
            response2 = send_and_receive(line, multiline=True)
            # Use the statement result if it succeeded, or if both failed
            # prefer the statement error (usually more relevant)
            if not response2.get("error", ""):
               response = response2
            elif "SyntaxError" not in response.get("error", ""):
               response = response2

      # Display the response
      stdout_text = response.get("stdout", "")
      result_text = response.get("result", "")
      error_text = response.get("error", "")

      if stdout_text:
         print(stdout_text, end="")
         if not stdout_text.endswith("\n"):
            print()
      if result_text:
         print(result_text)
      if error_text:
         print(f"\033[38;2;204;102;102m{error_text}\033[0m")  # Firebrick red for errors


if __name__ == "__main__":
   main()
