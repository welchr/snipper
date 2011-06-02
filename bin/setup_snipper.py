#!/usr/bin/env python

#===============================================================================
# This file is part of the Snipper program. 
# 
# Copyright (C) 2010 Ryan Welch
# 
# Snipper is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# Snipper is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#===============================================================================

import sys
import os.path as path
import os
import re
import time
import platform
from textwrap import fill
from ConfigParser import ConfigParser

# Setup virtual environment. 
def createVirtualEnv():
  script_dir = os.path.dirname(os.path.realpath(sys.argv[0]));
  env_dir = os.path.join(script_dir,"../pyenv");
  venv_script_path = os.path.abspath("./virtualenv.py");  
  
  cmd = "%s %s %s" % (sys.executable,venv_script_path,env_dir);
  code = os.system(cmd);

  if code != 0:
    sys.exit("Error: cannot continue, virtualenv installation failed.");
    
  d = {};
  
  if sys.platform == 'win32':
    binary = os.path.join(env_dir,"Scripts/python.exe");
    easy_install = os.path.join(env_dir,"Scripts/easy_install.exe");
  else:
    binary = os.path.join(env_dir,"bin/python");
    easy_install = os.path.join(env_dir,"bin/easy_install");

  d['binary'] = os.path.abspath(binary);
  d['easy_install'] = os.path.abspath(easy_install);
  return d;

def install_deps():
  print "Installing dependencies for Snipper..";

  # Create virtual environment. 
  env_paths = createVirtualEnv();

  # Install required packages. 
  os.system("%s %s" % (env_paths['easy_install'],"sphinx")); 
  
def test_deps():
  # Try to import sqlite3. 
  try:
    import sqlite3
  except:
    sys.exit("Error: your installed python interpreter does not have "
            "sqlite installed (tried to import sqlite3 and failed.) This "
            "is unusual. Please ask your administrator to try rebuilding "
            "the latest python interpreter, or if on a personal machine, "
            "re-download and install the latest python 2.x interpreter "
            "from http://www.python.org/. ");

def findConf():
  script_dir = path.abspath(path.dirname(sys.argv[0]));
  conf = path.abspath(path.join(script_dir,path.sep.join(['..','conf',"snipper.conf"])));
  
  if not path.isfile(conf):
    sys.exit("Cannot find snipper.conf location, tried: " % conf);

  return conf;

def getEmail():
  print "";
  print fill("NCBI requires that users of their Entrez Utilities encode "
        "their email address in the query string. This allows NCBI to "
        "contact you should there be a problem with the number of "
        "queries, or the type of queries, that you are executing. "
        "Please enter your email address below. ")
  print "";
  
  email = raw_input("Email: ");
  return email;

def updateConfEmail(conf,email):
  p = ConfigParser();
  p.read(conf);
  
  if not p.has_section('user'):
    p.add_section('user');
    
  p.set('user','email',email);

  f = open(conf,'w');
  p.write(f);
  f.close();

def check_python():
  pyv = sys.version_info;
  major = pyv[0];
  minor = pyv[1];
  version_string = ".".join([str(i) for i in sys.version_info[0:3]]);
  
  if major >= 2:
    if major < 3:
      if minor >= 6:
        return True;
      else: 
        sys.exit(fill("Error: your python version must be greater than 2.6 on the "
                      "2.0 branch, you are running: %s. Please download an updated "
                      "2.x python interpreter from http://www.python.org/download/. "
                      % version_string));
    else:
      sys.exit(fill("Error: you are currently running the 3.0 branch of python, "
                    "Snipper requires the 2.x branch of python. Please download "
                    "from http://www.python.org/download/. "));
  else:
    sys.exit(fill("Error: you appear to be running a very early release of "
                  "python - please download the latest 2.x interpreter from "
                  "http://www.python.org/download/. "));

def main():
  conf = findConf();
  
  # Check the python interpreter for the proper version. 
  check_python();

  # Install Snipper dependencies. 
  install_deps();
  test_deps();

  # WRONG! This is the developer's email, NOT the end user's email. 
  # See: http://www.ncbi.nlm.nih.gov/books/NBK25497/
#  # Ask user for their email address to encode in the NCBI queries. 
#  email = getEmail();
#  updateConfEmail(conf,email);

  print "Setup complete. Snipper is ready for launch - use bin/snipper.py.";

if __name__ == "__main__":
  main();
