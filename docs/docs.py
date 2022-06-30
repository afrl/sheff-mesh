import argparse      # For parsing command line arguments
import os
import shutil

def build():
    """Builds the documentation"""

    # Check that requirements exist
    if shutil.which("doxygen") is None:
        raise Exception("Doxygen was not found in the path and is required to build the documentation." \
                        "Binaries can be downloaded from https://www.doxygen.nl/download.html")

    if shutil.which("dot") is None:
        raise Exception("Dot was not found in the path and is required to build the documentation." \
                        "Binaries can be downloaded from https://graphviz.org/download/")


    if shutil.which("doxyrest") is None:
        raise Exception("Doxyrest was not found in the path and is required to build the documentation." \
                        "Binaries can be downloaded from " \
                        "https://github.com/vovkos/doxyrest/releases/tag/doxyrest-2.1.0""")

    if 'DOXYREST_DIR' not in os.environ:
        raise Exception("DOXYREST_DIR was not defined in the environment variables and is required " \
                        "to build the documentation.  Createa an environment variable that points to " \
                        "the installation directory of Doxyrest.")

    if shutil.which("sphinx-build") is None:
        raise Exception("sphinx-build was not found in the path and is required to build the documentation." \
                        "It can be installed with 'pip install sphinx-build'")


    os.system("doxygen Doxyfile")
    os.system("doxyrest -c doxyrest-config.lua")

    shutil.copyfile('conf.py', './rst_dir/conf.py')
    os.system("sphinx-build -b html rst_dir sphinx_html")

def view():
    if not os.path.exists("./sphinx_html/index.html"):
        raise Exception("The documentation was not build yet for this module.")

    os.chdir('./sphinx_html')
    os.system("index.html")

parser = argparse.ArgumentParser(description='A utility for building documentation for BSAM++.')

command_group = parser.add_argument_group('Commands')
command_group.add_argument('--build',
                           action='store_true',
                           help='Build documentation (requires: Doxygen, Dot, Doxyrest, sphinx-build and sphinx_rtd_theme)')

command_group.add_argument('--view',
                           action='store_true',
                           help='Views documentation in the default webbrowser.')
        
args = parser.parse_args()
if args.build:
    build()
if args.view:
    view()