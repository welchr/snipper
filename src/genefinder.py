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

import os
import sys
import warnings
import shelve
import cPickle
import intersection
import sqlite3
import pdb
from constants import *
from gene import *
from settings import *
from util import *

# Finding genes near SNPs.
class GeneFinder:
  def __init__(self):
    self.genes = dict();
    self.nearest_genes = dict();
  
  def addGene(self,gene,snp_or_region,dist,updown,chr,txStart,txEnd):
    # Add info to our internal tree. 
    node = self.genes.setdefault(gene,{}).setdefault(snp_or_region,{});
    node['dist'] = dist;
    node['updown'] = updown;
    node['txStart'] = txStart;
    node['txEnd'] = txEnd;
    node['chr'] = chr;

    # Get the distance to nearest SNP (seen so far.) 
    nearest = self.nearest_genes.setdefault(gene,None);
    
    # Figure out if this SNP is closer than any seen previously. 
    if nearest == None or dist < nearest:
      self.nearest_genes[gene] = dist;

  # Returns a sorted list of genes by their distance to the nearest SNP. 
  def getGenesByNearest(self):
    return sorted(self.nearest_genes,key=lambda x: self.nearest_genes[x]);

  def loadGeneDB(self):
    # Now, load them with their SNP data. 
    for gene in self.genes:
      for thing in self.genes[gene]:
        node = self.genes[gene][thing];
        g = Gene.valueOf(gene);
        if isinstance(thing,ChromRegion):
          g.addRegion(thing,node['updown']);
        else:
          g.addSNP(
            thing,
            node['dist'],
            node['updown']
          );
        
#        g.chrom = node['chr'];
#        g.txStart = node['txStart'];
#        g.txEnd = node['txEnd'];
        
  def find(self,snps,dist,num):
    raise NotImplementedError;

class SQLiteGeneFinder(GeneFinder):
  def __init__(self,db_file):
    GeneFinder.__init__(self);
    
    self.db_file = db_file;

    # Connect to database file.     
    try:
      self.db_con = sqlite3.connect(db_file);
    except:
      print >> sys.stderr, ("Error: could not connect to SQLite database. "
                            "Check permissions on the file: %s " % db_file);
      raise;
    
    # Insert "if" function for sqlite (I can't believe they don't have this!) 
    def sql_if(cond,true,false):
      if cond:
        return true;
      else:
        return false;
    
    self.db_con.create_function("if",3,sql_if);
  
    def sql_sign(x):
      if x > 0:
        return 1;
      elif x < 0:
        return -1;
      else:
        return 0;
    
    self.db_con.create_function("sign",1,sql_sign);
    
    def sql_least(*args):
      return min(args);
     
    self.db_con.create_function("least",-1,sql_least);
  
  def getPos(self,gene_symb):
    statement = """
    select 
      geneName,
      chrom,
      min(txStart) as txStart,
      max(txEnd) as txEnd 
    from refFlat 
    where 
      geneName='%(gene)s' and 
      chrom not like "%%\_%%" escape '\\' 
    group by geneName
    """ % { 
     'gene' : gene_symb
    };
    
    cursor = self.db_con.execute(statement);
    
    chrom = None;
    txStart = None;
    txEnd = None;
    while 1:
      row = cursor.fetchone();
      if row == None:
        break;
  
      chrom = row[1];
      txStart = row[2];
      txEnd = row[3];
    
    return (chrom,txStart,txEnd);
  
  def findRegions(self,regions):
    for region in regions:
      statement = """
        SELECT
          g.geneName as gene,
          g.chrom as chrom,
          g.strand as strand,
          min(g.txStart) as txStart,
          max(g.txEnd) as txEnd
        FROM
          %(gene_table)s g
        WHERE
          g.chrom = '%(chr)s'
          AND g.txEnd > %(rstart)s AND g.txStart < %(rend)s
        GROUP BY g.geneName
        """ % {
               'chr' : region.get_chrom(),
               'rstart' : region.start,
               'rend' : region.end,
               'gene_table' : SQLITE_REFFLAT
              };

      warnings.simplefilter("ignore");

      # Parse result.
      cursor = self.db_con.execute(statement);

      warnings.resetwarnings();

      # Did we find any results?
      if cursor.rowcount == 0:
        print >> sys.stderr, "Warning: could not find genes near SNP %s.." % snp;

      while 1:
        row = cursor.fetchone();
        if row == None:
          break;

        gene = row[0];
        chr = row[1];
        txStart = row[3];
        txEnd = row[4];
        
        if region.start > txStart and region.end < txEnd:
          direction = "within gene";
        else:
          direction = "overlaps gene endpoint";

        self.addGene(gene,region,0,direction,chr,txStart,txEnd);
  
  def find1000G(self,snps,dist,num):
    for snp in snps:
      (chr,pos) = parse1000G(snp);
      if chr == None or pos == None:
        continue;
      
      chr = chrom2ucsc(chr); # ucsc chrom format
      
      statement = """
        SELECT
          g.geneName as nearest_gene,
          g.strand as strand,
          min(if(sign(g.txEnd - %(pos)i) * sign(g.txStart+1 - %(pos)i) <= 0,0,least(abs(g.txEnd - %(pos)i), abs(g.txStart+1 - %(pos)i)))) as dist_to_gene,
          min(if(%(pos)i BETWEEN txStart+1 and txEnd, 'within',if((%(pos)i < txStart and g.strand='+') or (%(pos)i > txEnd and g.strand='-'),
               'upstream',
               'downstream'))) as direction,
          min(g.txStart) as txStart,
          max(g.txEnd) as txEnd
        FROM
          %(gene_table)s g
        WHERE
          g.chrom = "%(chr)s"
          and g.txEnd  >= %(pos)i - %(radius)s
          and g.txStart < %(pos)i + %(radius)s
        GROUP BY g.geneName
        ORDER BY dist_to_gene
        LIMIT %(limit)s
        """ % {
               'chr' : chr,
               'pos' : pos,
               'radius' : str(dist), 
               'limit' : str(num), 
               'gene_table' : SQLITE_REFFLAT
              };

      warnings.simplefilter("ignore");

      # Parse result.
      cursor = self.db_con.execute(statement);

      warnings.resetwarnings();

      # Did we find any results?
      if cursor.rowcount == 0:
        print >> sys.stderr, "Warning: could not find genes near SNP %s.." % snp;

      while (1):
        row = cursor.fetchone();
        if row == None:
          break;

        gene = row[0];
        dist_to_gene = row[2];
        direction = row[3];
        txStart = row[4];
        txEnd = row[5];

        self.addGene(gene,snp,dist_to_gene,direction,str(chr),txStart,txEnd);
  
  def find(self,snps,dist,num):
    for snp in snps:
      statement = """
        SELECT
          a.snp,
          a.chr,
          a.pos,
          g.geneName as nearest_gene,
          g.strand as strand,
          min(if(sign(g.txEnd - a.pos) * sign(g.txStart+1 - a.pos) <= 0,0,least(abs(g.txEnd - a.pos), abs(g.txStart+1 - a.pos)))) as dist_to_gene,
          min(if(pos BETWEEN txStart+1 and txEnd, 'within',if((pos < txStart and g.strand='+') or (pos > txEnd and g.strand='-'),
               'upstream',
               'downstream'))) as direction,
          min(g.txStart) as txStart,
          max(g.txEnd) as txEnd
        FROM
          %(snp_table)s a,
          %(gene_table)s g
        WHERE
          a.snp = "%(snp)s"
          and g.chrom = a.chr
          and g.txEnd  >= a.pos - %(radius)s
          and g.txStart < a.pos + %(radius)s
        GROUP BY a.snp,g.geneName
        ORDER BY a.snp,dist_to_gene
        LIMIT %(limit)s
        """ % {
               'snp' : snp, 
               'radius' : str(dist), 
               'limit' : str(num), 
               'snp_table' : SQLITE_SNP_POS, 
               'gene_table' : SQLITE_REFFLAT
              };

      warnings.simplefilter("ignore");

      # Parse result.
      cursor = self.db_con.execute(statement);

      warnings.resetwarnings();

      # Did we find any results?
      if cursor.rowcount == 0:
        print >> sys.stderr, "Warning: could not find genes near SNP %s.." % snp;

      while (1):
        row = cursor.fetchone();
        if row == None:
          break;

        # 0 - SNP
        # 1 - chrom
        # 3 - gene
        # 4 - strand
        # 5 - distance
        # 6 - direction
        row_snp = row[0];
        chr = row[1];
        gene = row[3];
        dist_to_gene = row[5];
        direction = row[6];
        txStart = row[7];
        txEnd = row[8];

        self.addGene(gene,row_snp,dist_to_gene,direction,chr,txStart,txEnd);

#class FlatGeneFinder(GeneFinder):
#  def __init__(self,snp_db,gene_file):
#    GeneFinder.__init__(self);
#
#    self.snpdb = None;
#    self.gene_tree = None;
#    
#    # Load SNP position database, and gene interval tree files. 
#    # Both of these must be constructed using the setup.py script found in <snipperpath>/bin. 
#  
#    if not os.path.isfile(snp_db):
#      raise ValueError, "Cannot find SNP database file: %s" % snp_db;
#  
#    if not os.path.isfile(gene_file):
#      raise ValueError, "Cannot find gene database file: %s" % gene_file;  
#  
#    self._load_snpdb(snp_db);
#    self._load_gene_tree(gene_file);
#
#  def _load_snpdb(self,dbfile):
#    self.snpdb = shelve.open(dbfile);
#  
#  def _load_gene_tree(self,gene_file):
#    self.gene_tree = cPickle.load(open(gene_file,"rb"));
#  
#  def __del__(self):
#    if self.snpdb != None:
#      self.snpdb.close();
#
#  def _region_lookup(self,region):
#    chr = region.chr;
#    start = region.start;
#    end = region.end;
#    
#    node = self.gene_tree.get(chr2ucsc(chr));
#    if node == None:
#      return;
#
#    # Genes on positive strand.      
#    overlaps = node['+'].find(start,end);
#    for gene in overlaps:
#      updown = None;
#      dist = None;
#
#      if gene.start < start and gene.end > start:
#        dist = 0;
#        updown = "overlapping-endpoint";
#      elif gene.start < end and gene.end > end:
#        dist = 0;
#        updown = "completely-overlapping";
#      elif gene.start > start and gene.end < end:
#        dist = 0;
#        updown = "completely-contained";
#      elif gene.start > start and gene.end > end:
#        dist = 0;
#        updown = "overlapping-endpoint";
#
#      self.addGene(gene.value,region,dist,updown);
#
#    # Genes on negative strand. 
#    overlaps = node['-'].find(start,end);
#    for gene in overlaps:
#      updown = None;
#      dist = None;
#
#      if gene.start < start and gene.end > start:
#        dist = 0;
#        updown = "overlapping-endpoint";
#      elif gene.start < end and gene.end > end:
#        dist = 0;
#        updown = "completely-overlapping";
#      elif gene.start > start and gene.end < end:
#        dist = 0;
#        updown = "completely-contained";
#      elif gene.start > start and gene.end > end:
#        dist = 0;
#        updown = "overlapping-endpoint";
#
#      self.addGene(gene.value,region,dist,updown);
#
#  def _gene_lookup(self,snp,chr,pos,dist,num):
#    results = [];
#    
#    node = self.gene_tree.get(chr);
#    if node == None:
#      return;
#
#    # Genes on positive strand.      
#    overlaps = node['+'].find(pos-dist,pos+dist);
#    for gene in overlaps:
#      updown = None;
#      dist = None;
#
#      if gene.start < pos and gene.end > pos:
#        dist = 0;
#        updown = 'within';
#      elif pos > gene.end:
#        dist = pos - gene.end;
#        updown = 'downstream';
#      else:
#        dist = gene.start - pos;
#        updown = 'upstream';
#
#      results.append( (dist,updown,gene.value) );
#
#    # Genes on negative strand. 
#    overlaps = node['-'].find(pos-dist,pos+dist);
#    for gene in overlaps:
#      updown = None;
#      dist = None;
#
#      if gene.start < pos and gene.end > pos:
#        dist = 0;
#        updown = 'within';
#      elif pos > gene.end:
#        dist = pos - gene.end;
#        updown = 'upstream';
#      else:
#        dist = gene.start - pos;
#        updown = 'downstream';
#
#      results.append( (dist,updown,gene.value) );
#
#    results = sorted(results);
#    if num > len(results):
#      num = len(results);
#
#    for i in xrange(num):
#      self.addGene(results[i][2],snp,results[i][0],results[i][1]);
# 
#  def find(self,snps,dist,num):
#    for snp in snps:
#      result = self.snpdb.get(snp);
#      if result != None:
#        (chr,pos) = result;
#        self._gene_lookup(snp,chr,pos,dist,num);
#      else:
#        print >> sys.stderr, "Warning: could not find position for SNP %s in database.." % snp; 
#  
#  def findRegions(self,regions):
#    for region in regions:
#      self._region_lookup(region);
      
