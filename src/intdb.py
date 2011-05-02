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

import shelve
import sys
import re
from prettytable import PrettyTable
from constants import *
from util import terminal_size, make_list

# Class to represent a MiMI interaction. 
class GeneInteraction():
  def __init__(self,gene1,gene2):
    self.gene1 = gene1;
    self.gene2 = gene2;
    self.gene_sym1 = None;
    self.gene_sym2 = None;
    self.int_info = {};
    self.components = set();
    self.functions = set();
    self.provenance = set();
    self.processes = set();
    self.types = set();
    self.pmids = set();
    self.taxon_id = None;
  
  def __eq__(self,other):
    if self.gene1 == other.gene1 and self.gene2 == other.gene2:
      return True;
    else:
      return False;
  
  def get_gene_pair(self):
    return (self.gene1,self.gene2);
  
  # Set taxonomy ID. 
  def set_tax(self,id):
    self.taxon_id = id;

  # Set symbols for each gene. 
  def set_symbol(self,id,symbol):
    if id == self.gene1:
      self.gene_sym1 = symbol;
    elif id == self.gene2:
      self.gene_sym2 = symbol;
    else:
      raise ValueError, "Error: ID (%s) supplied does not"\
                        " match (%s) or (%s)" % (id,self.gene1,self.gene2);
    
  def add_component(self,component,go_id):
    if component not in (None,''):
      self.components.add((component,go_id));
  
  def add_function(self,function,go_id):
    if function not in (None,''):
      self.functions.add((function,go_id));
  
  def add_interaction_type(self,type):
    if type not in (None,""):
      self.types.add(type);
  
  def add_provenance(self,source):
    self.provenance.add(source);
  
  def add_pubmed(self,pmid):
    self.pmids.add(pmid);
  
  def add_process(self,process,go_id):
    if process not in (None,''):
      self.processes.add((process,go_id));
  
  def write(self,out=sys.stdout):
    print >> out, str(self);
  
  def __str__(self):
    s = ""
    
    gene_sym1 = self.gene_sym1;
    gene_sym2 = self.gene_sym2;
    if gene_sym1 == None:
      gene_sym1 = "Symbol unknown";
    if gene_sym2 == None:
      gene_sym2 = "Symbol unknown";
  
    s += "%s (ID:%s) <--> %s (ID:%s): \n" % (self.gene_sym1,self.gene1,self.gene_sym2,self.gene2);
    
    if len(self.components) > 0:
      s += "-- Components: %s\n" % ", ".join(["%s [GO:%s]" % (i,j) for i,j in self.components]);
      
    if len(self.functions) > 0:
      s += "-- Functions: %s\n" % ", ".join(["%s [GO:%s]" % (i,j) for i,j in self.functions]);
      
    if len(self.processes) > 0:
      s += "-- Processes: %s\n" % ", ".join(["%s [GO:%s]" % (i,j) for i,j in self.processes]);
      
    if len(self.pmids) > 0:
      s += "-- PubMed IDs: %s\n" % ", ".join(self.pmids);
      
    if len(self.provenance) > 0:
      s += "-- Provenance: %s\n" % ", ".join(self.provenance);
    
    if self.taxon_id != None:
      s += "-- Taxonomy ID: %s\n" % self.taxon_id;
    
    return s;
  
  def search(self,terms):
    data = dict();
    terms = make_list(terms);
    
    if self.gene_sym1 and self.gene_sym2:
      key = "Interaction %s <--> %s" % (self.gene_sym1,self.gene_sym2);
    else:
      key = "Interaction %s <--> %s" % (self.gene1,self.gene2);
    
    for term in terms:
      re_string = "\\b" + str(term) + "\\b";
      pattern = re.compile(re_string,re.I);
      
      for proc in self.processes:
        if pattern.search(proc[0]):
          data.setdefault(key + ' Processes',[]).append(term);
          break;
      
      for comp in self.components:
        if pattern.search(comp[0]):
          data.setdefault(key + ' Components',[]).append(term);
          break;
      
      for func in self.functions:
        if pattern.search(func[0]):
          data.setdefault(key + ' Functions',[]).append(term);
          break;
    
    return data;

class InteractionDB:
  _GENE_ROOT = "gene_interactions";
  _PROT_ROOT = "protein_interactions";
  
  def __init__(self,file=None):
    self.root = dict();
    self.root[InteractionDB._GENE_ROOT] = {};
    self.root[InteractionDB._PROT_ROOT] = {};
    self.gene_interactions = [];  
    
  def add_gene_interactions(self,interactions):
    interactions = make_list(interactions);
    for inter in interactions:
      if isinstance(inter,GeneInteraction):
        (gene1,gene2) = inter.get_gene_pair();
        self.root[InteractionDB._GENE_ROOT].setdefault(gene1,{}).setdefault(gene2,inter);
        self.root[InteractionDB._GENE_ROOT].setdefault(gene2,{}).setdefault(gene1,inter);
        self.gene_interactions.append(inter);
  
  def has_gene_interactions(self):
    if len(self.gene_interactions) > 0:
      return True;
    else:
      return False;
  
  def get_interactions_for_gene(self,gene):
    return self.root[InteractionDB._GENE_ROOT].get(gene);
  
  def iter_gene_interactions(self):
    for i in self.gene_interactions:
      yield i;
  
  def direct_gene_interactions(self,genes):
    if not hasattr(genes,'__iter__'):
      raise ValueError;
    
    genes = list(set(genes));
    
    inters = InteractionDB();
    gene_root = self.root[InteractionDB._GENE_ROOT];
    for i in xrange(len(genes)-1):
      for j in xrange(i+1,len(genes)):
        level1 = gene_root.get(genes[i]);
        if level1 != None:
          level2 = level1.get(genes[j]);
          if level2 != None:
            inters.add_gene_interactions(level2);
    
    return inters;
    
  def print_term(self):
    """
    >>> from mimi import mimi_fetch_interactions
    >>> from intdb import InteractionDB
    >>> test = mimi_fetch_interactions('1436');
    >>> db = InteractionDB();
    >>> db.add_gene_interactions(test);
    >>> db.print_term();
    """
    
    fields = [
      'Gene 1',
      'Gene 2',
      'Components',
      'Processes',
      'Functions',
      'Types',
      'Provenance',
      'PMIDs'
    ];
    
    table = PrettyTable(fields);
    for f in fields:
      table.set_field_align(f,'l');
     
    fracs = [0.05,0.05,0.15,0.15,0.15,0.15,0.15,0.15];
    term_width = int(terminal_size()[0]*1.35)
      
    widths = [];
    for i in fracs:
      widths.append(int(term_width * i));
    
    table.set_field_widths(widths);
    
    gene_root = self.root[InteractionDB._GENE_ROOT];
    for uid1 in gene_root:
      uid1_root = gene_root.get(uid1);
      for uid2 in uid1_root:
        inter = uid1_root.get(uid2);
        row = [
          uid1,
          uid2,
          ["%s [GO:%s]" % (i,j) for i,j in inter.components],
          ["%s [GO:%s]" % (i,j) for i,j in inter.processes],
          ["%s [GO:%s]" % (i,j) for i,j in inter.functions],
          list(inter.types),
          ",".join(inter.provenance),
          ",".join(inter.pmids)
        ];
        table.add_wrapped_row(row);
    
    table.printt(fields=['Gene 1','Gene 2','Types','Provenance','PMIDs']);
    table.printt(fields=['Gene 1','Gene 2','Processes','Components','Functions']);
  
  def print_rest(self,out=sys.stdout):
    fields = [
      'Gene 1',
      'Gene 2',
      'Components',
      'Processes',
      'Functions',
      'Types',
      'Provenance',
      'PMIDs'
    ];
    
    table = PrettyTable(fields);
    for f in fields:
      table.set_field_align(f,'l');
    
    gene_root = self.root[InteractionDB._GENE_ROOT];
    pmids = set();
    for uid1 in gene_root:
      uid1_root = gene_root.get(uid1);
      for uid2 in uid1_root:
        inter = uid1_root.get(uid2);
        
        gene1 = None;
        gene2 = None;
        
        if inter.gene_sym1 != None:
          gene1 = ":ref:`%s`" % inter.gene_sym1;
        else:
          gene1 = uid1;
        
        if inter.gene_sym2 != None:
          gene2 = ":ref:`%s`" % inter.gene_sym2;
        else:
          gene2 = uid2;
        
        row = [
          [gene1],
          [gene2],
          ["* %s [GO:%s]" % (i,j) for i,j in inter.components],
          ["* %s [GO:%s]" % (i,j) for i,j in inter.processes],
          ["* %s [GO:%s]" % (i,j) for i,j in inter.functions],
          ["* %s" % i for i in inter.types],
          ["* %s" % i for i in inter.provenance],
          ["* %s_" % i for i in inter.pmids]
        ];
        table.add_multi_row(row);
        pmids.update( set(inter.pmids) );
    
    print >> out, "The following table lists direct interactions between the protein products of genes found near any input region:"
    print >> out, "";
    print >> out, table.get_string(fields=['Gene 1','Gene 2','Types','Provenance','PMIDs'],rest=True);
    print >> out, "";
    print >> out, "For these interactions, the following ontology terms were found:"
    print >> out, "";
    print >> out, table.get_string(fields=['Gene 1','Gene 2','Processes','Components','Functions'],rest=True);
    print >> out, "";
    
    for pmid in pmids:
      print >> out, ".. _%s: %s" % (pmid,PUBMED_WEB_URL + pmid);
    
    print >> out, "";
    
def main():
  import doctest
  doctest.testmod();
  
if __name__ == "__main__":
  main();    
    
    
  