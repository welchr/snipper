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
import os

script_dir = os.path.dirname(os.path.realpath(sys.argv[0]));

# Find virtualenv python interpreter.
if sys.platform == 'win32':
  pyex = os.path.join(script_dir,"../pyenv/Scripts/python.exe");
else:
  pyex = os.path.join(script_dir,"../pyenv/bin/python");

if not os.path.exists(pyex):
  print >> sys.stderr, "Error: could not find python interpreter at snipper/pyenv. Did you run the setup script bin/setup_snipper.py?";
  sys.exit(1);

# Fix args to have quotes. 
args = [];
for i in sys.argv[1:]:
  if i[0] == "-":
    args.append(i);
    continue;

  args.append('"' + str(i) + '"');
  
# First, try to find snipper at a relative location to where this script is located. 
snipper = os.path.join(script_dir,"../src/main.py");
if os.path.isfile(snipper):
  e = pyex + " -OO " + snipper + " " + " ".join(args);
  os.system(e);

# That apparently failed, so now let's try the environment variable "SNIPPER_PATH". 
snipper_path = os.environ.get("SNIPPER_PATH");
if snipper_path != None:
  snipper_abspath = os.path.join(snipper_path,"src/main.py");
  e = pyex + " -OO " + snipper_abspath + " " + " ".join(args);
  os.system(e);
