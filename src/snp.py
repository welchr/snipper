#!/usr/bin/env python
import sqlite3
from constants import *
from util import *

class SNPDB:
  def __init__(self,db_file):
    self.db_file = db_file;
    try:
      self.db_con = sqlite3.connect(db_file);
    except:
      print >> sys.stderr, ("Error: could not connect to SQLite database. "
                            "Check permissions on the file: %s " % db_file);
      raise;
  
  def get_pos(self,snp):
    if "chr" in snp:
      # This is a 1000G style SNP - get chrom/position from it. 
      temp = snp.replace("chr","").split(":");
      if len(temp) != 2:
        print >> sys.stderr, "Error: cannot recognize SNP format for %s, skipping.." % snp;
        return (None,None);
      
      chr = safe_int(temp[0]);
      pos = safe_int(temp[1]);
      
      return (chr,pos);
      
    else:
      query = "select * from %s where snp='%s'" % (SQLITE_SNP_POS,snp);
      cur = self.db_con.execute(query);
      
      chr = None;
      pos = None;
      while 1:
        row = cur.fetchone();
        if row == None:
          break;
        
        chr = row[0];
        if "_" in chr:
          continue;
        
        pos = row[1];
    
    return (chr,pos);
  
  def get_current_name(self,snp):
    query = "select rs_current from %s where rs_orig='%s'" % (SQLITE_REFSNP_TRANS,snp);
    cur = self.db_con.execute(query);
    
    result = cur.fetchone();
    if result == None:
      return snp;
    else:
      return result[0];
    