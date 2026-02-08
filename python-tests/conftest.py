# conftest.py - pytest configuration for Coot tests
# Run with: ccp4-python -m pytest . --coot=/path/to/coot
from __future__ import print_function
import pytest
import sys
import os
import glob
import subprocess
import tempfile
import json
import ast

# Test*.py files need special handling - see pytest_ignore_collect and pytest_collect_file hooks below
_test_star_files = glob.glob(os.path.join(os.path.dirname(__file__), "Test*.py"))

# Add paths for coot modules
test_dir = os.path.dirname(os.path.abspath(__file__))
python_dir = os.path.join(os.path.dirname(test_dir), "python")


# =============================================================================
# Dynamic collection of unittest tests from Test*.py files
# =============================================================================

def extract_unittest_tests(filepath):
    """
    Parse a Test*.py file and extract test class and method names without importing.
    Returns list of (class_name, method_name) tuples.
    """
    tests = []
    try:
        with open(filepath, 'r') as f:
            source = f.read()
        tree = ast.parse(source)
        
        for node in ast.walk(tree):
            if isinstance(node, ast.ClassDef):
                # Check if it's a unittest.TestCase subclass
                class_name = node.name
                for item in node.body:
                    if isinstance(item, ast.FunctionDef) and item.name.startswith('test'):
                        tests.append((class_name, item.name))
    except Exception as e:
        print("Warning: Could not parse %s: %s" % (filepath, e))
    return tests


class CootUnittestItem(pytest.Item):
    """A pytest Item that runs a single unittest test method inside coot."""
    
    def __init__(self, name, parent, test_file, class_name, method_name):
        super().__init__(name, parent)
        self.test_file = test_file
        self.class_name = class_name
        self.method_name = method_name
    
    def runtest(self):
        """Run the test via coot subprocess."""
        coot_path = self.config.getoption("--coot")
        coot_timeout = self.config.getoption("--coot-timeout")
        
        if coot_path:
            coot_path = os.path.expandvars(coot_path)
        if not coot_path or not os.path.exists(coot_path):
            import shutil
            coot_path = shutil.which("coot")
        if not coot_path or not os.path.exists(coot_path):
            pytest.skip("coot executable not found. Use --coot=/path/to/coot")
        
        # Create script to run single unittest test
        result_file = tempfile.mktemp(suffix='.json')
        module_name = os.path.splitext(os.path.basename(self.test_file))[0]
        
        script_content = '''
# Auto-generated unittest runner
from __future__ import print_function
import sys
import os
import json
import unittest

# Add test directory to path
test_dir = "{test_dir}"
if test_dir not in sys.path:
    sys.path.insert(0, test_dir)

python_dir = os.path.join(os.path.dirname(test_dir), "python")
if os.path.isdir(python_dir) and python_dir not in sys.path:
    sys.path.insert(0, python_dir)

result_file = "{result_file}"
result = {{"success": False, "message": ""}}

try:
    import coot
    import coot_utils
    import coot_testing_utils
except ImportError as e:
    result["message"] = "Failed to import coot modules: " + str(e)
    with open(result_file, 'w') as f:
        json.dump(result, f)
    coot.coot_real_exit(1)

try:
    # Import the test module
    from {module_name} import {class_name}
    
    # Create test suite with single test
    suite = unittest.TestSuite()
    suite.addTest({class_name}("{method_name}"))
    
    # Run with capturing output
    import io
    stream = io.StringIO()
    runner = unittest.TextTestRunner(stream=stream, verbosity=2)
    test_result = runner.run(suite)
    
    output = stream.getvalue()
    
    if test_result.wasSuccessful():
        result["success"] = True
        result["message"] = output
    else:
        result["success"] = False
        # Collect failure/error messages
        messages = []
        for test, traceback in test_result.failures + test_result.errors:
            messages.append(traceback)
        result["message"] = "\\n".join(messages) if messages else output

except Exception as e:
    import traceback
    result["message"] = str(e) + "\\n" + traceback.format_exc()

with open(result_file, 'w') as f:
    json.dump(result, f)

coot.coot_real_exit(0 if result["success"] else 1)
'''.format(
            test_dir=test_dir,
            result_file=result_file,
            module_name=module_name,
            class_name=self.class_name,
            method_name=self.method_name
        )
        
        # Write script to temp file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.py', delete=False) as f:
            f.write(script_content)
            script_path = f.name
        
        try:
            # Run coot with the script
            proc_result = subprocess.run(
                [coot_path, "--no-graphics", "--script", script_path],
                capture_output=True,
                text=True,
                timeout=coot_timeout
            )
            
            # Read result from file
            if os.path.exists(result_file):
                with open(result_file, 'r') as f:
                    data = json.load(f)
                    if not data.get("success", False):
                        raise CootTestFailure(data.get("message", "Test failed"))
            else:
                raise CootTestFailure(
                    "Result file not created - coot may have crashed.\n"
                    "stdout: %s\nstderr: %s" % (proc_result.stdout[-1000:], proc_result.stderr[-1000:])
                )
                
        except subprocess.TimeoutExpired:
            raise CootTestFailure("Test timed out after %d seconds" % coot_timeout)
        finally:
            for f in [script_path, result_file]:
                try:
                    if os.path.exists(f):
                        os.unlink(f)
                except:
                    pass
    
    def repr_failure(self, excinfo):
        """Custom failure representation."""
        if isinstance(excinfo.value, CootTestFailure):
            return str(excinfo.value)
        return super().repr_failure(excinfo)
    
    def reportinfo(self):
        return self.fspath, None, "%s::%s::%s" % (
            os.path.basename(self.test_file), self.class_name, self.method_name
        )


class CootTestFailure(Exception):
    """Custom exception for coot test failures."""
    pass


class CootUnittestFile(pytest.File):
    """A pytest collector for Test*.py unittest files."""
    
    def collect(self):
        """Collect test items from the file."""
        tests = extract_unittest_tests(str(self.fspath))
        for class_name, method_name in tests:
            name = "%s::%s" % (class_name, method_name)
            yield CootUnittestItem.from_parent(
                self,
                name=name,
                test_file=str(self.fspath),
                class_name=class_name,
                method_name=method_name
            )


def pytest_ignore_collect(collection_path, config):
    """Don't ignore Test*.py - we handle them specially in pytest_collect_file."""
    # Return None to let normal collection logic apply (but pytest_collect_file intercepts first)
    return None


def pytest_collect_file(parent, file_path):
    """Hook to collect Test*.py files as coot unittest files.
    
    This hook runs early and returns our custom collector, which prevents
    pytest from trying to import Test*.py as normal Python modules.
    """
    if file_path.suffix == ".py" and file_path.name.startswith("Test"):
        return CootUnittestFile.from_parent(parent, path=file_path)
    return None


def pytest_addoption(parser):
    """Add --coot command line option"""
    parser.addoption(
        "--coot",
        action="store",
        default=None,
        help="Path to coot executable (e.g., /usr/bin/coot or coot)"
    )
    parser.addoption(
        "--coot-timeout",
        action="store",
        default=60,
        type=int,
        help="Timeout in seconds for each coot test (default: 60)"
    )


def pytest_configure(config):
    """Configure pytest markers"""
    config.addinivalue_line("markers", "slow: marks tests as slow")
    config.addinivalue_line("markers", "requires_coot: marks tests that require coot runtime")
    config.addinivalue_line("markers", "requires_ccp4: marks tests that require CCP4")


@pytest.fixture(scope="session")
def coot_path(request):
    """Get the coot executable path"""
    path = request.config.getoption("--coot")
    if path is not None:
        # Expand environment variables like $CCP4
        path = os.path.expandvars(path)
    if path is None or not os.path.exists(path):
        # Try to find coot in PATH
        import shutil
        path = shutil.which("coot")
    if path is None or not os.path.exists(path):
        pytest.skip("coot executable not found. Use --coot=/path/to/coot")
    return path


@pytest.fixture(scope="session")
def coot_timeout(request):
    """Get the timeout for coot tests"""
    return request.config.getoption("--coot-timeout")


@pytest.fixture(scope="session")
def coot_runner(coot_path, coot_timeout):
    """Fixture that provides a function to run code in coot"""
    def run_in_coot(test_code, test_name="test"):
        """
        Run Python code inside coot and return (success, output, error).
        
        The test_code should set a variable `_test_result` to True/False
        and optionally `_test_message` for failure details.
        """
        # Create a temporary result file
        result_file = tempfile.mktemp(suffix='.json')
        
        # Create a temporary script
        script_content = '''
# Auto-generated test script
from __future__ import print_function
import sys
import json

_test_result = False
_test_message = ""
_result_file = "{result_file}"

try:
    import coot
    import coot_utils
except ImportError as e:
    _test_message = "Failed to import coot: " + str(e)
    with open(_result_file, 'w') as f:
        json.dump({{"success": False, "message": _test_message}}, f)
    coot.coot_real_exit(1)

try:
{indented_code}
except Exception as e:
    import traceback
    _test_result = False
    _test_message = str(e) + "\\n" + traceback.format_exc()

# Write result to file
result = {{"success": _test_result, "message": _test_message}}
with open(_result_file, 'w') as f:
    json.dump(result, f)

# Exit coot
coot.coot_real_exit(0 if _test_result else 1)
'''.format(
            result_file=result_file,
            indented_code="\n".join("    " + line for line in test_code.split("\n"))
        )

        # Write to temp file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.py', delete=False) as f:
            f.write(script_content)
            script_path = f.name

        try:
            # Run coot with the script
            result = subprocess.run(
                [coot_path, "--no-graphics", "--script", script_path],
                capture_output=True,
                text=True,
                timeout=coot_timeout
            )
            
            output = result.stdout + result.stderr
            
            # Read result from file
            success = False
            message = ""
            
            if os.path.exists(result_file):
                try:
                    with open(result_file, 'r') as f:
                        data = json.load(f)
                        success = data.get("success", False)
                        message = data.get("message", "")
                except (json.JSONDecodeError, IOError) as e:
                    message = "Failed to read result file: %s" % e
            else:
                message = "Result file not created - coot may have crashed"
                # Check if there's useful info in output
                if "Error" in output or "error" in output:
                    message += "\nOutput: " + output[-500:]
            
            return success, output, message
            
        except subprocess.TimeoutExpired:
            return False, "", "Test timed out after %d seconds" % coot_timeout
        except Exception as e:
            return False, "", str(e)
        finally:
            # Clean up temp files
            for f in [script_path, result_file]:
                try:
                    if os.path.exists(f):
                        os.unlink(f)
                except:
                    pass
    
    return run_in_coot


@pytest.fixture(scope="session")
def unittest_data_dir():
    """Provide the test data directory"""
    d = os.environ.get("COOT_TEST_DATA_DIR")
    if not d:
        home = os.path.expanduser("~")
        d = os.path.join(home, "data", "greg-data")
    return d


# For backwards compatibility - these will skip if coot not available as module
@pytest.fixture(scope="session")  
def coot_module(coot_runner):
    """Provide coot module - only works when running via coot binary"""
    try:
        import coot
        return coot
    except ImportError:
        pytest.skip("coot module not available - use --coot flag to run tests")


@pytest.fixture(scope="session")
def coot_utils_module(coot_runner):
    """Provide coot_utils module"""
    try:
        import coot_utils
        return coot_utils
    except ImportError:
        pytest.skip("coot_utils module not available - use --coot flag to run tests")
