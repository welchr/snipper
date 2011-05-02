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

from collections import namedtuple
from gene import *

class RankH1:
  def __init__(self):
    self.eqtl_score = 1000 / 10;
    self.eqtl_dist_score = 50 / 100000;
    self.num_sterm_score = 750 / 10;
    self.num_sterm_matches = 400 / 10;
    self.num_sterm_pubmed = 400 / 100;
    self.sterm_omim = 100;
    self.num_inters = 250 / 100;
    self.dist_to_snp = 50 / 100000;
    self.intdb = None;
  
  # Database containing direct interactions among genes to be ranked. 
  def setInteractionDB(self,intdb):
    self.intdb = intdb;
  
  # Ranks genes. Returns list of genes in order by rank. 
  def rank_genes(self,genes):
    scores = [];
    for gene in genes:
      gene_obj = Gene.valueOf(gene);
      score = 0;
      
      # Score for # of eQTLs linked to gene. 
      score += self.eqtl_score * len(gene_obj.eqtls);
      
      # Score for number of search terms that matched the gene's info. 
      score += self.num_sterm_matches * len(gene_obj.terms);
      
      # Score for number of places that search terms matched inside gene info. 
      # Also score for if OMIM was matched. 
      terms = set();
      for where,term_set in gene_obj.terms.iteritems():
        terms.update(term_set);
        if where == Gene.OMIM:
          score += self.sterm_omim;
      score += self.num_sterm_score * len(terms);
      
      # Score for number of interactions with other genes. 
      if self.intdb != None:
        inters = self.intdb.get_interactions_for_gene(gene);
        if inters != None:
          score += self.num_inters * len(inters);
      
      scores.append( (score,gene) );
    
    gene_order = [i[1] for i in sorted(scores,reverse=True)];
    return gene_order;
