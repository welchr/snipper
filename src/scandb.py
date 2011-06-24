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

import httplib
import urllib
import sys
import decimal
import prettytable
from constants import *
from util import make_list
from BeautifulSoup import BeautifulSoup

SCANDB_HOST = "www.scandb.org:80";

class eQTL:
  def __init__(self,snp,gene):
    self.snp = snp;
    self.gene = gene;
    self.pval = None;
    self.tissue = None;
    self.population = None;
    self.organism = None;
    self.source = None;
  
  def __str__(self):
    #format_str = "SNP: %-12s Gene: %-10s P-value: %-.2e   Tissue: %-10s Population: %-4s Organism: %-10s Source: %-10s "
    format_str = "SNP: %s, Gene: %s, P-value: %-.2e, Tissue: %s, Population: %s, Organism: %s, Source: %s"
    return format_str % (
      str(self.snp),
      str(self.gene),
      self.pval,
      str(self.tissue),
      str(self.population),
      str(self.organism),
      str(self.source),
    )
  
  # Return a string representing the "group" this eQTL belongs to.
  # For example, when displaying a list of eQTLs, it's likely that the user
  # wants to sort first by groups (tissue/population/organism/database source),
  # and then by p-value. 
  def group_str(self):
    return "%s,%s,%s,%s" % (
      str(self.tissue),
      str(self.population),
      str(self.organism),
      str(self.source)
    );
  
  # Determines uniqueness of eQTL. An eQTL is defined by the following members.
  # Note that p-value is not used - two eQTLs with exactly the same data
  # but different p-values are still considered equal. 
  def __hash__(self):
    return hash(
      str(self.snp) +
      str(self.gene) +
      str(self.tissue) +
      str(self.population) +
      str(self.source) +
      str(self.organism)
    );

class eQTLSet():
  def __init__(self):
    self.qtls = set();
    self.snp_ptrs = {};
    self.gene_ptrs = {};
  
  def add(self,eqtl):
    eqtl = make_list(eqtl);
    for q in eqtl:
      if isinstance(q,eQTL):
        self.qtls.add(q);
        self.snp_ptrs.setdefault(q.snp,set()).add(q);
        self.gene_ptrs.setdefault(q.gene,set()).add(q);
  
  def forGene(self,gene):
    return self.gene_ptrs.get(gene);
  
  def forSNP(self,snp):
    return self.snp_ptrs.get(snp);
    
  def numGenes(self):
    return len(self.gene_ptrs);
  
  def numSNP(self):
    return len(self.snp_ptrs);
  
  def getGenes(self):
    genes = set();
    for q in self.qtls:
      genes.add(q.gene);
    
    return genes;
  
  def write_rest_table(self,out=sys.stdout):
    fields = ['eQTL','Gene','Tissue','Population','Organism','P-value','Source'];
    t = prettytable.PrettyTable(fields);
    
    for f in fields:
      t.set_field_align(f,'l');
    
    for snp,eqtls in self.snp_ptrs.iteritems():
      for q in eqtls:
        t.add_row([
          q.snp,
          ":ref:`%s`" % q.gene,
          q.tissue,
          q.population,
          q.organism,
          str(q.pval),
          q.source
        ]);
    
    print >> out, t.get_string(hrules=prettytable.ALL,rest=True);
  
  # Return a sorted list of eQTLs. This sorts by group and then by p-value, so
  # clusters of eQTLs of the same "group" come out together, and are sorted
  # by p-value within those groups. 
  def sort(self):
    temp = list(self.qtls);
    temp.sort(key = lambda x: float(x.pval));
    temp.sort(key = lambda x: x.group_str());
    return temp;
  
  def write(self):
    for qtl in self.qtls:
      print qtl;

  def __len__(self):
    return len(self.qtls);
    
  def __iter__(self):
    return self.qtls.__iter__();

# Finds genes associated with given SNPs from ScanDB. 
def scandb_snps(snps,pval):
  snps = make_list(snps);
  params = urllib.urlencode({
    'list' : ",".join(snps),
    'snpinfo' : 1,
    'expr' : 1,
    'pval' : pval,
    'output' : 'tab'
  });
  headers = {
    "Content-type" : "application/x-www-form-urlencoded",
    "Accept": "text/plain"
  };
  
  # Submit query to ScanDB. 
  try:
    conn = httplib.HTTPConnection(SCANDB_HOST,timeout=CON_TIMEOUT);
    conn.request("POST","/newinterface/snpinfo.php",params,headers);
    response = conn.getresponse();
    data = response.read();
    conn.close();
  except:
    print >> sys.stderr, "Error: query to ScanDB failed. The message was:";
    print >> sys.stderr, str(sys.exc_info()[1]);
    print >> sys.stderr, "ScanDB itself may be down, or your internet connection may have failed.";
    conn.close();
    return []; # exit
  
  # If the response status wasn't OK, there's something wrong. Bail out. 
  if response.status != 200:
    print >> sys.stderr, "Error: query to ScanDB failed. The response was:";
    print >> sys.stderr, "%s %s" % (str(response.status),str(response.reason));
    return [];
  # If the query itself failed for some reason, bail out.
  elif data == None:
    return []; # exit
  
  # Checks passed - let's parse our data. 
  parsed = eQTLSet();
  for line in data.split("\n")[1:]:
    if line == "":
      continue;
    
    e = line.split("\t");
    snp = e[0];
    
    for q in e[6].split(":"):
      qs = q.split();
      if qs[0] != "NA":
        eqtl = eQTL(snp,qs[0]);
        eqtl.population = qs[1];
        eqtl.source = "SCAN";
        eqtl.organism = "Homo sapiens"
        eqtl.tissue = "LCL";
        
        # Need to be careful in casting p-value. If it's not a float, it could
        # crash the whole program. 
        try:
          eqtl.pval = float(qs[2]);
        except:
          eqtl.pval = "NA";
        
        parsed.add(eqtl);
  
  return parsed;

def scandb_genes(symbols=set()):
  pass

def main():
  a = eQTL("rs1","X");
  a.pval = "1.4e-04";
  a.source = "ScanDB";
  a.population = "CEU";
  a.tissue = "Liver";
  
  b = eQTL("rs2","Z");
  b.pval = "2.9e-07";
  b.source = "ScanDB";
  b.population = "CEU";
  b.tissue = "Liver";
  
  c = eQTL("rs3","J");
  c.pval = "3.7e-09";
  c.source = "ScanDB";
  c.population = "CEU";
  c.tissue = "Liver";
  
  d = eQTL("rs4","D");
  d.pval = "7.4e-12";
  d.source = "ScanDB";
  d.population = "CEU";
  d.tissue = "Adipose";
  
  e = eQTL("rs5","U");
  e.pval = "4.5e-03";
  e.source = "ScanDB";
  e.population = "CEU";
  e.tissue = "Adipose";
  
  eset = eQTLSet();
  for qtl in [a,b,c,d,e]:
    eset.add(qtl);
  
  for qtl in sorted([a,b,c,d,e],key=lambda x: x.value()):
    print qtl;
    
  # inject into global namespace for inspection
  globals().update(locals());

if __name__ == "__main__":
  main();
  