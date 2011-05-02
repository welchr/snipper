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

from optparse import OptionParser, SUPPRESS_HELP
from textwrap import *
import sys

# Slightly modified OptionParser to make the help output look reasonable.
class VerboseParser(OptionParser):
  def print_help(self):
    # Print usage.
    print "usage: %s [options]" % sys.argv[0];
    print "";

    # Print options.
    for option in self.option_list:
      if option.help == "SUPPRESSHELP":
        continue;

      if option.type != None:
        print fill(", ".join(option._short_opts + option._long_opts) + " <%s>" % option.type,\
          initial_indent="  ",subsequent_indent="  ");
      else:
        print fill(", ".join(option._short_opts + option._long_opts),\
          initial_indent="  ",subsequent_indent="  ");

      for line in option.help.split("\n"):
        print fill(line,initial_indent="    ",subsequent_indent="    ");
      print "";