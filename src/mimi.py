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

import urllib2
import sys
import re
from constants import *
from intdb import *
from xml.etree.cElementTree import ElementTree

# MiMI interaction fetch URL. 
# The %s should be filled in with a Entrez Gene ID. 
MIMI_INT_URL = "http://mimi.ncibi.org/MimiWeb/fetch.jsp?geneid=%s&type=interactions";

# Given a Entrez gene ID, return a list of interactions with that gene. 
def mimi_fetch_interactions(gene_id,taxid=None):
  gene_id = str(gene_id);
  
  url = MIMI_INT_URL % gene_id;
  
  if _SNIPPER_DEBUG:
    print "DEBUG: executing MiMI URL %s" % url;
  
  xml = urllib2.urlopen(url,timeout=CON_TIMEOUT);
  tree = ElementTree();
  tree.parse(xml);
  
  go_pattern = re.compile("(.+) \[GO:(\d+)\]");
  def extract(pattern,string):
    match = pattern.search(string);
    if match:
      return match.groups();
    else:
      return (None,None);
  
  results = [];
  for int_gene in tree.getroot().findall("MiMI/Response/ResultSet/Result/InteractingGene"):
    other_gene = int_gene.find("GeneID").text;
    interaction = GeneInteraction(gene_id,other_gene);
    
    for element in int_gene.getchildren():
      if element.tag == "TaxonomyID":
        interaction.set_tax(element.text);
      elif element.tag == "InteractionAttribute":
        type = element.get('type');
        if type == "Component":
          tup = extract(go_pattern,element.text);
          interaction.add_component(*tup);
        elif type == "Function":
          tup = extract(go_pattern,element.text);
          interaction.add_function(*tup);
        elif type == "Process":
          tup = extract(go_pattern,element.text);
          interaction.add_process(*tup);
        elif type == "Provenance":
          interaction.add_provenance(element.text);
        elif type == "PubMed":
          interaction.add_pubmed(element.text);
        elif type == "InteractionType":
          interaction.add_interaction_type(element.text);
    
    # Taxonomy ID filter. 
    if taxid != None:
      if interaction.taxon_id != taxid:
        continue;
    
    results.append(interaction);
    
  return results;
        

