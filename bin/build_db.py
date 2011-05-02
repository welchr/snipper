#!/usr/bin/env python
import sys
import os.path as path
import os
import gzip
import re
import csv
import sqlite3
import time
import platform
import urllib
import urllib2
from ConfigParser import ConfigParser

# Constants. 
SQLITE_SNP_POS = "snp_pos";
SQLITE_TRANS = "refsnp_trans";
SQLITE_REFFLAT = "refFlat";
RS_MERGE_ARCH_URL = "ftp://ftp.ncbi.nih.gov/snp/database/organism_data/human_9606/RsMergeArch.bcp.gz";

# Add snipper source to import path. 
script_dir = os.path.dirname(os.path.realpath(sys.argv[0]));
src_dir = path.abspath(path.join(script_dir,path.sep.join(['..','src'])));
sys.path.append(src_dir);

# Snipper imports.
from intersection import *
from verboseparser import VerboseParser

# Parse command line args. 
def get_settings():
  parser = VerboseParser();

  build_help = "Genome build to download tables from. This uses UCSC's " +\
    "naming convention, i.e. hg17, hg18, etc. Defaults to the latest human genome build.";
  parser.add_option("--build",dest="build",help=build_help);
  opts,args = parser.parse_args();

  # Check build convention.
  if opts.build != None and opts.build[0:2] != 'hg':
    print >> sys.stderr, "Error: build must look like \'hg19\', \'hg18\', etc."
    sys.exit(1);

  return (opts,args);

# Extract a filename from a full path (removes directories and extension.) 
def filename(path):
  return os.path.splitext(os.path.basename(path))[0];

def remove_noerror(file):
  try:
    os.remove(file);
  except:
    pass;

# Function to see if a given URL actually exists. 
def exists(url):
  try:
    urllib2.urlopen(url);
  except urllib2.HTTPError, e:
    if e.code == 404:
      return False;

  return True;

def dl_hook(count,block_size,total_size):
  percent = count*block_size*100.0/total_size;
  sys.stdout.write("\r%.1f%%" % percent)
  sys.stdout.flush()

class UCSCManager:
  UCSC_MAIN_URL = "http://hgdownload.cse.ucsc.edu/goldenPath";
  UCSC_FTP_URL = "ftp://hgdownload.cse.ucsc.edu/goldenPath/";
  
  @staticmethod
  def getLatestHumanBuild():
    latest_hg = None;
    
    resp = urllib2.urlopen(UCSCManager.UCSC_FTP_URL);
    lines = resp.readlines();
    dirs = [i.rstrip().split()[8] for i in lines];
    p = re.compile("hg(\d+)");
    hg = filter(lambda x: p.search(x) != None,dirs);
    hg_versions = map(lambda x: int(p.search(x).groups()[0]),hg);
    latest_hg = sorted(hg_versions,reverse=True)[0];
    
    return "hg" + str(latest_hg);
  
  @staticmethod
  def getLatestSNPTable(build):
    p = re.compile("snp(\d+?).sql");
    resp = urllib.urlopen(UCSCManager.UCSC_MAIN_URL + "/" + build + "/" + "database");
    tables = set();
    for line in resp:
      m = p.search(line);
      if m != None:
        table = "snp" + str(m.groups()[0]);
        tables.add(table);

    return max(tables);

  @staticmethod
  def downloadLatestSNPTable(dir,build):
    latest_table = UCSCManager.getLatestSNPTable(build);

    url = "/".join([UCSCManager.UCSC_MAIN_URL,build,'database',latest_table + '.txt.gz']);

    file = path.join(dir,latest_table + ".gz");
    #progress = urlgrabber.progress.TextMeter();
    #grabber = urlgrabber.grabber.URLGrabber(progress_obj=progress,timeout=30);
    
    urllib.urlretrieve(url,file);
    #grabber.urlgrab(url,file);
    
    return file;

  @staticmethod
  def downloadLatestRefFlat(dir,build):
    url = "/".join([UCSCManager.UCSC_MAIN_URL,build,'database','refFlat.txt.gz']);
    file = path.join(dir,'refFlat_' + build + '.gz');
    
    #progress = urlgrabber.progress.TextMeter();
    #grabber = urlgrabber.grabber.URLGrabber(progress_obj=progress,timeout=30);
    urllib.urlretrieve(url,file,reporthook=dl_hook);
    #grabber.urlgrab(url,file);
 
    return file;

  @staticmethod
  def download_snp_table(dir,build,table):
    url = "/".join([UCSCManager.UCSC_MAIN_URL,build,'database',table + '.txt.gz']);
    file = path.join(dir,table + ".gz");

    if not exists(url):
      print >> sys.stderr, "Could not find SNP table %s at UCSC - check your table name." % table;
      print >> sys.stderr, "URL attempted was: %s" % url;
      sys.exit(1);

    try:
      urllib.urlretrieve(url,file,reporthook=dl_hook);
      #progress = urlgrabber.progress.TextMeter();
      #grabber = urlgrabber.grabber.URLGrabber(progress_obj=progress,timeout=30);
      #grabber.urlgrab(url,file);
    except IOError:
      print >> sys.stderr, "A network connection to the UCSC data repository could not be made.";
      sys.exit(1);

    return file;

def find_genome_path():
  script_dir = path.abspath(path.dirname(sys.argv[0]));
  genome_dir = path.abspath(path.join(script_dir,path.sep.join(['..','data','genome'])));
    
  if path.isdir(genome_dir):
    return genome_dir;

  else:
    sys.exit("Cannot find genome directory path - please see installation instructions.");

def get_our_snp_version(conf,build):
  p = ConfigParser();
  p.read(conf);
  try:
    snp_table = p.get(build,'snp_version');
  except:
    snp_table = None;
    
  return snp_table;

def update_conf(conf,build,db_file,snp_table_version):
  p = ConfigParser();
  p.read(conf);

  try:
    p.add_section(build);
  except:
    pass

  db_file_path = "data/genome/";
  db_file_path += os.path.split(db_file)[1];

  p.set(build,'db_file',db_file_path)
  p.set(build,'snp_version',snp_table_version);
  
  f = open(conf,'w');
  p.write(f);
  f.close();
  
def find_conf():
  script_dir = path.abspath(path.dirname(sys.argv[0]));
  conf = path.abspath(path.join(script_dir,path.sep.join(['..','conf',"snipper.conf"])));
  
  if not path.isfile(conf):
    sys.exit("Cannot find snipper.conf location, tried: " % conf);

  return conf;

#def fix_paths(path_list):
#  new_list = [];
#  for p in path_list:
#    if p == '':
#      new_list.append(p);
#      continue;
#
#    if not os.path.exists(p):
#      continue;
#
#    new_list.append(p);
#
#  return new_list;

#def find_relative(file):
#  full_path = None;
#  
#  # Find the root, using the script's location. 
#  start_loc = os.path.realpath(sys.argv[0]);
#  script_dir = None;
#  if os.path.isdir(start_loc):
#    script_dir = start_loc;
#  else:
#    script_dir = os.path.dirname(start_loc);
#  root_dir = os.path.join(script_dir,"../");
#  
#  # If the file to find has a path, it means it is a path relative
#  # to the root. We need to attach that path to the root. 
#  (file_path,file_name) = os.path.split(file);
#  if file_path != "":
#    root_dir = os.path.abspath(os.path.join(root_dir,file_path));
#  
#  if file_name == "" or file_name == None:
#    if os.path.exists(root_dir):
#      full_path = root_dir;
#  else:
#    temp_path = os.path.join(root_dir,file_name);
#    if os.path.exists(temp_path):
#      full_path = temp_path;
#  
#  return full_path;

def quote(x):
  if hasattr(x,"__iter__"):
    return ['"%s"' % str(i) for i in x];
  else:
    return '"%s"' % x;

def load_snp_table(snp_table_file,db_file):
  if path.splitext(snp_table_file)[1] == ".gz":
    f = gzip.open(snp_table_file);
  else:
    f = open(snp_table_file);
  
  # Open database file. Create it, if it doesn't exist already. 
  try:
    db = sqlite3.connect(db_file);
  except:
    sys.exit("Error: could not open db_file for writing: %s. Do you have permissions to write here?" % db_file);
  
  # Drop original table if it exists. 
  db.execute("DROP TABLE IF EXISTS %s" % SQLITE_SNP_POS);
  
  # Create table. 
  columns = ['chr','pos','snp'];
  types = ['TEXT','INTEGER','TEXT'];
  spec = ", ".join(map(lambda x: " ".join(x),zip(columns,types)));
  cmd = "CREATE TABLE %s ( %s );" % (SQLITE_SNP_POS,spec);
  db.execute(cmd);

  # Load values into table, saving the set of SNPs for later use 
  # in the setup.
  snp_set = set();
  reader = csv.reader(f,delimiter="\t");
  for line in reader:
    values = ",".join([quote(line[i]) for i in (1,3,4)]);
    db.execute("INSERT INTO %s VALUES (%s)" % (
      SQLITE_SNP_POS,
      values
    ));
    
    snp_set.add(line[4]);
  
  # Create indices on proper columns. 
  db.execute("CREATE INDEX index_chrpos ON %s (%s)" % (SQLITE_SNP_POS,"chr,pos"));
  db.execute("CREATE INDEX index_snp ON %s (%s)" % (SQLITE_SNP_POS,"snp"));
  
  return snp_set;

def load_refflat(refflat_file,db_file):
  if path.splitext(refflat_file)[1] == ".gz":
    f = gzip.open(refflat_file);
  else:
    f = open(refflat_file);
  
  # Open database file. Create it, if it doesn't exist already. 
  try:
    db = sqlite3.connect(db_file);
  except:
    sys.exit("Error: could not open db_file for writing: %s. Do you have permissions to write here?" % db_file);
  
  # Create table. 
  columns = ['geneName','chrom','strand','txStart','txEnd'];
  types = ['TEXT','TEXT','TEXT','INTEGER','INTEGER'];
  spec = ", ".join(map(lambda x: " ".join(x),zip(columns,types)));
  cmd = "CREATE TABLE %s ( %s );" % (SQLITE_REFFLAT,spec);
  db.execute(cmd);

  # Load values into table. 
  reader = csv.reader(f,delimiter="\t");
  for line in reader:
    values = ",".join([quote(line[i]) for i in (0,2,3,4,5)]);
    db.execute("INSERT INTO %s VALUES (%s)" % (
      SQLITE_REFFLAT,
      values
    ));
  
  # Create indices on proper columns. 
  db.execute("CREATE INDEX index_chrom_start_end ON %s (%s)" % (SQLITE_REFFLAT,"chrom,txStart,txEnd"));
  db.execute("CREATE INDEX index_geneName ON %s (%s)" % (SQLITE_REFFLAT,"geneName"));
  
  # Return location of database file. 
  return db_file;

class MergeHistory:
  def __init__(self):
    self.graph = {};

  def add_merge(self,source,target):
    self.graph[source] = target;

  def get_merge_target(self,source):
    return self.graph.get(source);

  def find_current(self,source):
    target = source;
    while 1:
      next_t = self.get_merge_target(target);
      if next_t != None:
        target = next_t;
      else:
        break;

    return target;

  def iter_node(self):
    for node in self.graph:
      yield node;

def parse_rsmerge(file,snp_set,snp_build):
  snp_build = int(snp_build.replace("snp",""));
  if os.path.splitext(file)[1] == ".gz":
    f = gzip.open(file);
  else:
    f = open(file);

  hist = MergeHistory();
  for line in f:
    e = line.rstrip().split("\t");
    build = int(e[2]);

    if build > snp_build:
      continue;
    
    rs_high = "rs" + e[0];
    rs_low = "rs" + e[1];
    
    hist.add_merge(rs_high,rs_low);

  refsnp_trans = {};
  for snp in hist.iter_node():
    if snp not in snp_set:
      latest = hist.find_current(snp);
      refsnp_trans[snp] = latest;
      
  return refsnp_trans;

def download_merge_arch(dir):
  file = os.path.join(dir,"RsMergeArch.bcp.gz");
  urllib.urlretrieve(RS_MERGE_ARCH_URL,file,reporthook=dl_hook);
  return file;

def load_snp_trans(rs_merge_file,snp_set,snp_table_version,db_file):
  snp_trans = parse_rsmerge(rs_merge_file,snp_set,snp_table_version);
  
  # Open database file. Create it, if it doesn't exist already. 
  try:
    db = sqlite3.connect(db_file);
  except:
    sys.exit("Error: could not open db_file for writing: %s. Do you have permissions to write here?" % db_file);
  
  # Drop original table if it exists. 
  db.execute("DROP TABLE IF EXISTS %s" % SQLITE_TRANS);
  
  # Create table. 
  columns = ['rs_orig','rs_current'];
  types = ['TEXT','TEXT'];
  spec = ", ".join(map(lambda x: " ".join(x),zip(columns,types)));
  cmd = "CREATE TABLE %s ( %s );" % (SQLITE_TRANS,spec);
  db.execute(cmd);

  # Load values into table. 
  for orig,cur in snp_trans.iteritems():
    values = ",".join(quote([orig,cur]));
    db.execute("INSERT INTO %s VALUES (%s)" % (
      SQLITE_TRANS,
      values
    ));
  
  # Create indices on proper columns. 
  db.execute("CREATE INDEX index_rs_orig ON %s (%s)" % (SQLITE_TRANS,"rs_orig"));
  db.execute("CREATE INDEX index_rs_current ON %s (%s)" % (SQLITE_TRANS,"rs_current"));

def main():
  (opts,args) = get_settings();
  conf = find_conf();
  genome_path = find_genome_path();
  
  # If no build was specified, find the latest human genome build from UCSC. 
  build = opts.build;
  if build == None:
    build = UCSCManager.getLatestHumanBuild();
    
  # Database file for this build. 
  db_file = path.join(genome_path,build + ".db");
  
  # Our current latest SNP table for this build. 
  #our_latest_snp_table = get_our_snp_version(conf,build);

  # If we already have the latest SNP table for this build, do nothing. 
  # Otherwise, grab the latest table from UCSC and install it. 
  print "Asking UCSC for latest SNP table in build %s.." % build;
  ucsc_latest_snp_table = UCSCManager.getLatestSNPTable(build);
  print "Found: %s" % ucsc_latest_snp_table;

  print "Downloading SNP table %s.." % ucsc_latest_snp_table;
  snp_table_file = UCSCManager.download_snp_table(genome_path,build,ucsc_latest_snp_table);
  print "\nFinished downloading %s.." % ucsc_latest_snp_table;

  print "Downloading refFlat table..";
  refflat = UCSCManager.downloadLatestRefFlat(genome_path,build);
  print "\nFinished downloading %s.." % refflat;

  print "Inserting SNP positions into database (this can take over 2 hours, please be patient!)..";
  snp_set = load_snp_table(snp_table_file,db_file);

  print "Inserting gene positions into database..";
  load_refflat(refflat,db_file);
  
  print "Downloading rs_merge history..";
  rsmergearch_file = download_merge_arch(genome_path);

  print "Creating SNP translation table from old build to current..";
  load_snp_trans(rsmergearch_file,snp_set,ucsc_latest_snp_table,db_file);

  print "Updating configuration file %s.." % conf;
  update_conf(conf,build,db_file,ucsc_latest_snp_table);

  # Cleanup. 
  remove_noerror(snp_table_file);
  remove_noerror(rsmergearch_file);
  remove_noerror(refflat);

if __name__ == "__main__":
  main();

