from cx_Freeze import setup, Executable
import os.path
PYTHON_INSTALL_DIR = os.path.dirname(os.path.dirname(os.__file__))
os.environ['TCL_LIBRARY'] = os.path.join(PYTHON_INSTALL_DIR, 'tcl', 'tcl8.6')
os.environ['TK_LIBRARY'] = os.path.join(PYTHON_INSTALL_DIR, 'tcl', 'tk8.6')
build_exe_options = {"packages": ["errno", "os", "re", "stat", "subprocess","collections", "pprint","shutil", "humanize","pycallgraph"],"includes": ['numpy.core._methods', 'numpy.lib.format','matplotlib.backends.backend_qt5agg']}

setup( name = "Standalone FlareNet PPC", 
       version = "0.1", 
       options = {"build_exe": build_exe_options}, 
       description = "FlareNet PPC V0.1",
	   author = 'Keyhan Babaee',
	   author_email = 'KeyhanBabaee@gmail.com',
       executables = [Executable("FlareNet_Particle_Penetration_Calculator.py")])