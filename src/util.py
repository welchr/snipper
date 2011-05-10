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

from BeautifulSoup import BeautifulStoneSoup as bss
from urllib2 import urlopen
from StringIO import StringIO
import sys
import re
import os
import warnings

def safe_int(x):
  try:
    return int(x);
  except:
    return None;

def mute_stderr():
  sys.stderr = StringIO();
  
def unmute_stderr():
  sys.stderr = sys.__stderr__;
  
def mute_std():
  sys.stderr = StringIO();
  sys.stdout = StringIO();

def unmute_std():
  sys.stderr = sys.__stderr__;
  sys.stdout = sys.__stdout__;

def markup_term(text,terms):
  if text == None or terms == None or len(terms) == 0:
    return text;
  else:
    for term in terms:
      p = re.compile(term,re.I);
      match = p.search(text);
      
      if match:
        matched_word = match.group();
        text = p.sub("\ ``%s``\ " % matched_word,text);
          
  return text;

# Return chrom/position from 1000G SNP. 
def parse1000G(snp):
  temp = snp.replace("chr","").split(":");
  if len(temp) != 2:
    return (None,None);
  
  chr = safe_int(temp[0]);
  pos = safe_int(temp[1]);
  
  return (chr,pos);

# Region specification. 
class ChromRegion:
  def __init__(self,chr=None,start=None,end=None,name=None):
    self.chr = chr;
    self.start = start;
    self.end = end;
    
    if name == None:
      if None not in (chr,start,end):
        self.name = str(self);
    else:
      self.name = name;
  
  def __str__(self):
    return "chr%s:%s-%s" % (str(self.chr),str(self.start),str(self.end));

  def __hash__(self):
    return hash(" ".join(map(str,[self.chr,self.start,self.end,self.name])));

  def __cmp__(self,other):
    if self.chr > other.chr:
      return 1;
    elif self.chr < other.chr:
      return -1;
    else:
      if self.start == other.start:
        return 0;
      elif self.start > other.start:
        return 1;
      else:
        return -1;

#  def __eq__(self,other):
#    if self.chr == other.chr and self.start == other.start and self.end == other.end:
#      return True;
#    else:
#      return False;
#  
#  def __ne__(self,other):
#    if self.chr != other.chr:
#      return True;
#    else:
#      if self.start != other.start:
#        return True;
#      elif self.end != other.end:
#        return True;
#      else:
#        return False;

  @staticmethod
  def from_str(string):
    region = ChromRegion();
    try:
      string = string.replace("chr","");
      chr = chrom2chr(string.split(":")[0]);
      start = string.split(":")[1].split("-")[0];
      end = string.split(":")[1].split("-")[1];
      
      start = start.replace(",","");
      end = end.replace(",","");
      
      start = int(start);
      end = int(end);
      
      region.chr = chr;
      region.start = start;
      region.end = end;
      region.name = str(region);
      
      return region;
    except:
      raise ValueError, "Could not correctly create region from string: %s" % string;
  
  @staticmethod
  def from_snp(snp,chr,pos):
    region = ChromRegion();
    region.name = snp;
    region.chr = safe_int(chr);
    region.start = safe_int(pos);
    region.end = safe_int(pos);
    return region;
  
  def get_chrom(self):
    if self.chr < 23:
      return "chr" + str(self.chr);
    elif self.chr == 23:
      return "chrX";
    elif self.chr == 24:
      return "chrY";

# Function for sorting a list of regions. 
def sort_regions(regions):
  def cmp_genome(x,y):
    if x.chr > y.chr:
      return 1;
    elif x.chr < y.chr:
      return -1;
    else:
      if x.start == y.start:
        return 0;
      elif x.start > y.start:
        return 1;
      else:
        return -1;
  
  return sorted(regions,cmp=cmp_genome);

# Singleton decorator. 
def singleton(cls):
    instance_container = []
    def getinstance():
        if not len(instance_container):
            instance_container.append(cls())
        return instance_container[0]
    return getinstance

@singleton
class ChromConverter:
  def __init__(self):
    self.pattern = re.compile("(chrom|chr|)(\w+)");
  
  def __call__(self,chr):
    if chr == None:
      return None;
    
    chr = str(chr);
    search = self.pattern.search(chr);
    chr_int = None;
    if search != None:
      (chr_string,chr_val) = search.groups();
      if chr_val != None:
        try:
          chr_int = int(chr_val);
        except:
          if chr_val == 'X':
            chr_int = 23;
          elif chr_val == 'Y':
            chr_int = 24;
          elif chr_val == 'mito':
            chr_int = 25;
          elif chr_val == 'XY':
            chr_int = 26;
          else:
            chr_int = None;
    
    return chr_int;

# Call this as if it were a function:
# chrom2chr('chrX');
# chrom2chr('chrom2');
chrom2chr = ChromConverter();

def chr2ucsc(c):
  c = int(c);
  if c < 23 and c > 0:
    return "chr%i" % c;
  elif c == 23:
    return 'chrX';
  elif c == 24:
    return 'chrY';
  elif c == 25:
    return 'chrmito';
  elif c == 26:
    return 'chrXY';
  else:
    return None;

def chrom2ucsc(chrom):
  try:
    chr = chrom2chr(chrom);
    ucsc = chr2ucsc(chr);
    return ucsc;
  except:
    return None;

# Deprecated decorator. 
def deprecated(func):
  """This is a decorator which can be used to mark functions
  as deprecated. It will result in a warning being emmitted
  when the function is used."""
  def newFunc(*args, **kwargs):
      warnings.warn("Call to deprecated function %s." % func.__name__,
                    category=DeprecationWarning)
      return func(*args, **kwargs)
  newFunc.__name__ = func.__name__
  newFunc.__doc__ = func.__doc__
  newFunc.__dict__.update(func.__dict__)
  return newFunc

# Courtesy of the Python mailing list..
def ioctl_GWINSZ(fd):                  #### TABULATION FUNCTIONS
     try:                                ### Discover terminal width
         import fcntl, termios, struct, os
         cr = struct.unpack('hh',
                            fcntl.ioctl(fd, termios.TIOCGWINSZ, '1234'))
     except:
         return None
     return cr

def terminal_size():
     ### decide on *some* terminal size
     # try open fds
     cr = ioctl_GWINSZ(0) or ioctl_GWINSZ(1) or ioctl_GWINSZ(2)
     if not cr:
         # ...then ctty
         try:
             fd = os.open(os.ctermid(), os.O_RDONLY)
             cr = ioctl_GWINSZ(fd)
             os.close(fd)
         except:
             pass
     if not cr:
         # env vars or finally defaults
         try:
             cr = (env['LINES'], env['COLUMNS'])
         except:
             cr = (25, 80)
     # reverse rows, cols
     return int(cr[1]), int(cr[0])

# Displays the message to the user and asks for yes/no confirmation. 
def confirm(message,out=sys.stderr):
  result = raw_input(message);
  result = result.strip();
  result = result.lower();

  if result == "" or result.lower() in ("yes","y"):
    return True;
  else:
    return False;

def soupify(url):
  return bss(urlopen(url));

def periodCheck(string):
  if string:
    if string[-1] != '.':
      string += '.';

  return string;

# Die with a message.                                                                                   
def die(msg):                                                                                           
  print msg;
  sys.exit(1);

# Iterator to flatten nested structures. 
def iter_flatten(iterable):
  it = iter(iterable)
  for e in it:
    if isinstance(e, (list, tuple)):
      for f in iter_flatten(e):
        yield f
    else:
      yield e

def isIterable(object):
  return hasattr(object,'__iter__');

def make_list(object):
  if not isIterable(object):
    return [ object ];
  else:
    return list(object);

def subsets(x,step):
  x = make_list(x);
  i = 0
  for j in xrange(step,len(x)+step,step):
    yield x[i:j]
    i += step

def ifnone(object,alternative):
  if object == None:
    return alternative;
  else:
    return object;

def convertFlank(flank):
  iFlank = None;

  try:
    iFlank = int(flank);
    return iFlank;
  except:
    pass;

  p = re.compile("(.+)(kb|KB|Kb|kB|MB|Mb|mb|mB)")
  match = p.search(flank);
  if match:
    digits = match.groups()[0];
    suffix = match.groups()[1];

    if suffix in ('kb','KB','Kb','kB'):
      iFlank = float(digits)*1000;
    elif suffix in ('MB','Mb','mb','mB'):
      iFlank = float(digits)*1000000;

    iFlank = int(round(iFlank));

  return iFlank;