#!/usr/bin/env python
import os
import sys

sys.path.insert(0,"../src");
from main import PROG_VERSION

tarfile = "snipper_source_%s.tar.gz" % PROG_VERSION;

if os.path.exists(tarfile):
  print "Removing existing tar.."
  os.remove(tarfile);

command = "tar zcfh %s" % tarfile;
paths = """
../bin
../conf
../doc
../data/genome
../example
../src
""".split();

args = [];

print "Options enabled: "
if 0:
  print ".. excluding database files";
  args.append("--exclude *.db");

if 1:
  print ".. adding snipper/ to beginning of each file path";
  args.append("--xform 's|^|snipper/|'");

if 1:
  print ".. excluding compiled python bytecode files";
  args.append("--exclude *.pyc");
  args.append("--exclude *.pyo");

final = "%s %s" % (command," ".join(paths + args));

print "Packing.."
os.system(final);