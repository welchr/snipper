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

import sys
import os
import re
import prettytable
import shelve
from time import sleep
from ncbi import *
from pubmed import *
from textwrap import *
from util import *
from constants import *
from scandb import *
from snp import *
 
class Gene:
  # References to gene objects. 
  _SYMB = dict();
  _UID = dict();
  
  # Constants for referring to locations in gene info. 
  SUMMARY = "Gene summary";
  GENERIF = "GeneRIF";
  PHENOTYPE = "Phenotype";
  PATHWAY = "Pathway";
  GOTERM = "GO Term";
  OMIM = "OMIM Text";
  PUBMED = "Pubmed";

  # Creates a new gene. You should use valueOf() instead
  # for caching purposes. 
  def __init__(self):
    # Member variables. 
    self.uid = None;
    self.desc = None;
    self.symb = None;
    self.syns = [];
    self.type = None;
    self.loc = None;
    self.omim_id = None;
    self.omim_text = None;
    self.summary = None;
    self.GO = [];
    self.phenotypes = [];
    self.pathways = dict();
    self.terms = dict();
    self.soup = None;
    self.omim_soup = None;
    self.snps = dict();
    self.regions = dict();
    self.eqtls = eQTLSet();
    self.pubmed_ids = set();
    self.pubmed_by_term = dict();
    self.total_pubmed = None;
    self.gene_rifs = [];
    self.user_requested = False;
    self.chrom = None;
    self.txStart = None;
    self.txEnd = None;
    self.genome_region = None;

  # Entrez UID. 
  def setUID(self,uid):
    self.uid = str(uid);

  def getUID(self):
    return self.uid;

  # Long description of gene (long name.)
  def setDesc(self,desc):
    self.desc = str(desc);

  # The official gene symbol used by NCBI. 
  def setSymbol(self,symb):
    self.symb = str(symb);
 
  def getSymbol(self):
    return self.symb;

  # The gene's type (protein-coding, etc.) 
  def setType(self,type):
    self.type = str(type);

  # Gene symbol synonyms. 
  def setSyns(self,syns):
    self.syns = [ str(g) for g in syns ];

  # Chromosomal location. 
  def setLoc(self,loc):
    self.loc = str(loc);
    
  # Set the OMIM ID linked to this gene. 
  def setOMIM(self,omim_id):
    self.omim_id = str(omim_id);
 
  def getOMIM(self):
    return self.omim_id;

  # Set the OMIM text linked to this gene. 
  def setOMIMText(self,omim_text):
    self.omim_text = str(omim_text);

  # Summary text for the gene. 
  def setSummary(self,summary_text):
    self.summary = str(summary_text);

  # Set the GO terms for this gene. Must be a list of strings. 
  def setGO(self,go_terms):
    self.GO = [ str(g) for g in go_terms ];

  # Phenotypes associated with this gene. 
  def setPhenotypes(self,phenotypes):
    self.phenotypes = [ str(g) for g in phenotypes ];

  # KEGG pathways.
  # Pathways should be a dictionary, where: 
  # pathways = {
  #   pathway_name : URL to KEGG pathway }
  def setKEGG(self,pathways):
    self.pathways = pathways;

  def addSearchTerm(self,where,term):
    self.terms.setdefault(where,set()).add(term);

  def getSearchTerms(self):
    return self.terms;

  def getSearchTermWords(self):
    words = set();
    for where,term_set in self.terms.iteritems():
      words.update(term_set);
    
    return list(words);
    
  # Searches information stored in this gene for a search term. 
  def search(self,term):
    # Return data. 
    data = dict();

    re_string = "\\b" + str(term) + "\\b";
    pattern = re.compile(re_string,re.I);

    # Search summary text. 
    if self.summary != None:
      if pattern.search(self.summary):
        data[Gene.SUMMARY] = term;

    # Search GeneRIFs. 
    boom = True;
    for pmid in self.gene_rifs:
      if boom:
        for text in self.gene_rifs[pmid]:
          if pattern.search(text):
            data[Gene.GENERIF] = term;
            boom = False;
            break;

    # Search phenotypes. 
    for pheno in self.phenotypes:
      if pattern.search(pheno):
        data[Gene.PHENOTYPE] = term;
        break;

    # Search KEGG. 
    for pathway in self.pathways:
      if pattern.search(pathway):
        data[Gene.PATHWAY] = term;
        break;

    # Search Gene Ontology. 
    for goterm in self.GO:
      if pattern.search(goterm):
        data[Gene.GOTERM] = term;
        break;

    # Search OMIM record. 
    if self.omim_soup:
      match = self.omim_soup.find(text=pattern);
      if match != None:
        data[Gene.OMIM] = term;

    # Install information into gene. 
    for where in data:
      self.addSearchTerm(where,data[where]);

    return data;

  def countTerms(self):
    unique_terms = set();
    for where in self.terms:
      for term in self.terms[where]:
        unique_terms.add(term);

    return len(unique_terms);

  # Store pubmed IDs linked to the gene. 
  def addPubmed(self,id):
    self.pubmed_ids.add(id);

  def addPubmedByTerm(self,id,term):
    self.pubmed_by_term.setdefault(term,set()).add(id);
    self.addSearchTerm("Pubmed",term);

  # Set the list of pubmed IDs linked to this gene. 
  def setPubmed(self,id_list):
    self.pubmed_ids = id_list;

  # How many articles? 
  def countPubmed(self):
    return len(self.pubmed_ids);
  
  def getPubmed(self):
    return self.pubmed_ids;

  # Set the total number of pubmed articles linked to the gene. 
  # This is not the same as the number of articles currently 
  # stored for the gene. 
  def setTotalPubmed(self,count):
    self.total_pubmed = count;

  def getTotalPubmed(self):
    return self.total_pubmed;

  # Caches the full soupified XML response from NCBI for the gene
  # in case someone else wants to mess with it. 
  def setSoup(self,soup):
    self.soup = soup;

  def getSoup(self):
    return self.soup;

  # Storage of OMIM soup as well. 
  def setOMIMSoup(self,soup):
    self.omim_soup = soup;

  def getOMIMSoup(self):
    return self.omim_soup;

  # Adds information about SNPs that are either near or inside this gene. 
  # If the distance is 0, it is considered to be within the gene. 
  # Otherwise, distance is how far up or downstream in megabases. 
  def addSNP(self,rsid,distance,updown):
    self.snps.setdefault(rsid,{})['distance'] = distance;
    self.snps.setdefault(rsid,{})['updown'] = updown;

  def addRegion(self,region,containment):
    self.regions[region] = containment;

  # Returns SNP information. 
  def countAllSNP(self):
    return len(self.snps);

  # Was this gene populated with data? 
  def isPopulated(self):
    status = True;
    if self.soup == None:
      status = False;

    return status;

  # Print this gene. Takes in a file object to print to. 
  # If nothing is passed, it assumes STDOUT. 
  def write(self,out=sys.stdout,wrap_lines=True):
    # For wrapping lines support, we need to know width of the terminal. 
    term_width = terminal_size()[0];

    # If we're writing to a file, set width to 80. 
    if out != sys.stdout: 
      term_width = 80;

    # Prefix for line "headers"
    prefix = "[+] ";

    # If request is to write a non-populated gene (not sure how this can happen!)
    # print something different. 
    if not self.isPopulated():
      print >> sys.stderr, "Warning: no information found for gene " + str(self.symb);
      return;

    # Description of the gene. 
    if self.desc != None:
      print >> out, prefix + "GENE: " + str(self.desc) + " [" + str(self.symb) + "]";
    else: 
      print >> out, prefix + "GENE: [" + str(self.symb) + "]";
    
    # UID etc. 
    print >> out, prefix + "Entrez Gene UID: " + str(self.uid);
    print >> out, prefix + "Location: " + str(self.loc);
    print >> out, prefix + "Type: " + str(self.type);
    
    # Synonyms. 
    if len(self.syns) > 0:
      out.write(prefix + "Synonyms: ");
      for i in self.syns:
        out.write(i + " ");
      out.write("\n");
    else:
      print >> out, prefix + "Synonyms: None";

    # Print search terms. 
    if len(self.terms) > 0:
      print >> out, "";
      print >> out, prefix + "Search terms matched: "
      for where in self.terms:
        joined = ", ".join(self.terms[where]);
        print >> out, "-- Location: %-15s Terms: %s" % (where,joined);

    # SNP information. 
    if len(self.snps) > 0:
      print >> out, "";
      print >> out, prefix + "SNPs provided by the user: ";
      for snp in sorted(self.snps,key=lambda x: int(self.snps[x]['distance'])):
        print >> out, "SNP: %-15s Distance (bp): %-13s Direction: %-10s" %\
          (snp,self.snps[snp]['distance'],self.snps[snp]['updown']);
    else:
      print >> out, "";
      print >> out, prefix + "SNPs within distance provided by user: None";
  
    if len(self.regions) > 0:
      print >> out, "";
      print >> out, prefix + "Regions provided by the user containing this gene: ";
      for region in self.regions:
        print >> out, str(region);
  
    # eQTL information.
    if len(self.eqtls) > 0:
      print >> out, "";
      print >> out, prefix + "eQTL associations:";
      for eqtl in self.eqtls:
        print >> out, str(eqtl);
  
    # Summary. 
    if self.summary != None:
      print >> out, "";
      for line in wrap(prefix + "Summary: " + self.summary,initial_indent="",subsequent_indent="",width=term_width):
        print >> out, line;
    else:
      print >> out, prefix + "Summary: None";

    # Phenotypes. 
    if len(self.phenotypes) > 0:
      print >> out, "";
      pheno_line = prefix + "Phenotypes: " + "\n".join(["\"" + i + "\"" for i in self.phenotypes]);
      for line in wrap(pheno_line,initial_indent="",subsequent_indent="",width=term_width):
        print >> out, line;
    else:
      print >> out, "";
      print >> out, prefix + "Phenotypes: None";

    # KEGG Pathways. 
    print >> out, "";
    if len(self.pathways) > 0:
      print >> out, prefix + "KEGG Pathways: ";
      for pathway in self.pathways:
        print >> out, "%s [ %s ]" % (pathway,self.pathways[pathway]);
    else:
      print >> out, prefix + "KEGG Pathways: None";

    # GO terms. 
    if len(self.GO) > 0:
      print >> out, "";
      go_line = prefix + "GO Terms: " + "\n".join(["\"" + i + "\"" for i in self.GO]);
      for line in wrap(go_line,initial_indent="",subsequent_indent="",width=term_width):
        print >> out, line;
    else:
      print >> out, "";
      print >> out, prefix + "GO Terms: None";

    # OMIM ID. 
    if self.omim_id != None:
      print >> out, "";
      print >> out, prefix + "OMIM: [" + self.omim_id + "] Link: " + OMIM_URL + "?id=" + self.omim_id;
    else:
      print >> out, "";
      print >> out, prefix + "OMIM: None";

    # OMIM text. 
    if self.omim_text != None:
      omim_line = prefix + "OMIM Text: " + self.omim_text;
      for line in wrap(omim_line,initial_indent="",subsequent_indent="",width=term_width):
        print >> out, line;

    # The number of Pubmed articles. 
    #if self.countPubmed() > 0:
      #print >> out, prefix + "Number of Pubmed articles linked to gene %s: %s" % (str(self.symb),str(self.getTotalPubmed()));
  
    # Print GeneRIFs. 
    if len(self.gene_rifs) > 0:
      print >> out, "\n" + prefix + "Gene references into function for %s:" % (str(self.symb));
      for pmid in self.gene_rifs:
        for text in self.gene_rifs[pmid]:
          if wrap_lines:
            for line in wrap("-- %s %s" % (text,"PMID: " + pmid),subsequent_indent="   ",width=term_width):
              print >> out, line;
          else:
            print >> out, "-- PMID: %-9s Text: %s" % (pmid,text);

    # Print Pubmed articles. 
    if self.countPubmed() > 0:
      print >> out, "\n" + prefix + "Top Pubmed articles linked to gene %s, by date: " % str(self.symb);
      #for pmid in self.pubmed_ids:
        #Pubmed.valueOf(pmid).write();
      Pubmed.writeAll(self.pubmed_ids,out);

      for term in self.pubmed_by_term:
        if term == "any":
          print >> out, "\n" + prefix + "Top Pubmed articles linked to gene %s matching any search term: " % str(self.symb);
        else:
          print >> out, "\n" + prefix + "Top Pubmed articles linked to gene %s matching search term \"%s\":" % (str(self.symb),term);
        
        #for pmid in self.pubmed_by_term[term]:
          #Pubmed.valueOf(pmid).write();
        Pubmed.writeAll(self.pubmed_by_term[term],out);

  # Print information for symbols. 
  @staticmethod
  def writeAll(symbols,out=sys.stdout):
    if not isIterable(symbols):
      raise ValueError, "Must provide list or set of symbols.";

    term_width = terminal_size()[0];
    if out != sys.stdout: 
      term_width = 80;

    # Populate information for genes first. 
    Gene.populate(symbols);

    symbols = list(symbols);
    for i in xrange(len(symbols)):
      symb = symbols[i];
      
      for i in xrange(term_width):
        out.write('=');
      out.write("\n\n");
      
      Gene.valueOf(symb).write(out);

      print >> out, "";

  @staticmethod
  def writeReST(symbols,db_file,search_terms=[],out=sys.stdout):
    print >> out, "Gene Information";
    print >> out, "================";
    print >> out, "";
    
    # Construct mapping from SNP/regions --> genes. 
    region2gene = {};
    for symb in symbols: 
      gene = Gene.valueOf(symb);
      for snp,info in gene.snps.iteritems():
        dist = info['distance'];
        region2gene.setdefault(snp,[]).append((symb,dist,'nearby'));
      for region in gene.regions:
        region2gene.setdefault(region.name,[]).append((symb,0,'nearby'));
    
    # Add in eQTL genes. 
    for symb in symbols:
      gene = Gene.valueOf(symb);
      for eqtl in gene.eqtls:
        region2gene.setdefault(eqtl.snp,[]).append((symb,999999999,'eqtl'));
    
    # Order each list of genes by distance to SNP (region distances are 0.) 
    # We also need to keep track of this list for genes printed out later. 
    for region in region2gene:
      gene_list = sorted(region2gene[region],key = lambda x: x[1]);
      region2gene[region] = gene_list;
    
    # Find ordering of SNPs/regions. 
    regions = [];
    snp_db = SNPDB(db_file);
    for region in region2gene:
      if "chr" in region and "-" in region:
        chrom_region = ChromRegion.from_str(region);
        regions.append(chrom_region);
      elif "rs" in region: 
        chrpos = snp_db.get_pos(region);
        chrom = chrom2chr(chrpos[0]);
        pos = chrpos[1];
        chrom_region = ChromRegion.from_snp(region,chrom,pos);
        regions.append(chrom_region);
      elif "chr" in region:
        chrpos = snp_db.get_pos(region);
        chrom = chrom2chr(chrpos[0]);
        pos = chrpos[1];
        chrom_region = ChromRegion.from_snp(region,chrom,pos);
        regions.append(chrom_region);
    
    sorted_regions = sort_regions(regions);
    
    # Find user-requested genes. 
    user_genes = [];
    for symb in symbols:
      gene = Gene.valueOf(symb);
      if gene.user_requested:
        user_genes.append(gene);
    
    # Sort user-requested genes by genome position. 
    user_genes = sorted(user_genes,key=lambda x: x.genome_region);
    
    # Gene summary table. 
    # Print out a description of the table. 
    symbols_ordered = [];
    if len(sorted_regions) > 0:
      print >> out, "Region Table";
      print >> out, "------------";
      print >> out, "";
      print >> out, "The table below lists each user-provided SNP or chromosomal region, sorted by position in the genome.\n";
      print >> out, "For each SNP, nearby genes are listed first, sorted by distance to the SNP. eQTL genes are listed after nearby genes.\n";
      print >> out, "Search terms that match information other than PubMed searches are listed individually. If at least 1 term matches an article in PubMed, the word \"pubmed\" is listed in the \"Search Terms Matched\" column.";
      print >> out, "";
      
      summary_t = prettytable.PrettyTable(['SNP/Region','Chrom','Nearby Gene','eQTL Gene','Search Terms Matched','PubMed Articles']);
      for region in sorted_regions:
        gene_list = region2gene[region.name];
        for item in gene_list:
          gene = Gene.valueOf(item[0]);
          symb = gene.getSymbol();
        
          # Remember order that symbols are listed in the summary table. But don't include duplicates!
          # Duplicates cause "duplicate target name" errors. 
          if symb not in symbols_ordered:
            symbols_ordered.append(symb); 
  
          #symb = symb + "_"; # ReST hyperlink format
          symb = ":ref:`%s`" % symb;        
  
          row = [""]*6;
          row[0] = region.name;
          row[1] = region.chr;
          if item[2] == 'nearby':
            row[2] = symb;
          elif item[2] == 'eqtl':
            row[3] = symb;
          
          # If "any" is in the search term list, replace that text with something
          # more descriptive. 
          search_terms = gene.getSearchTermWords();
  #        if "any" in search_terms and len(search_terms) == 1:
  #          search_terms = ['pubmed'];
  #        elif "any" in search_terms:
  #          ind = search_terms.index("any");
  #          del search_terms[ind];
          try:
            ind = search_terms.index("any");
            search_terms[ind] = "pubmed";
          except:
            pass;
  
          row[4] = ", ".join(search_terms);
          row[5] = str(gene.getTotalPubmed());
          summary_t.add_row(row);
      
      print >> out, summary_t.get_string(hrules=prettytable.ALL,rest=True)
      print >> out, "";
    
    # Print section for user-requested genes. 
    if len(user_genes) > 0:
      print >> out, "User-requested Genes";
      print >> out, "--------------------";
      print >> out, "";
      print >> out, ("The table below lists all genes that were explicitly "
                    "requested by the user using the -g or --gene command line "
                    "arguments. ");
      print >> out, "";
      
      user_t = prettytable.PrettyTable(['Chrom','Gene','Search Terms Matched','PubMed Articles']);
      for gene in user_genes:
        symb = gene.getSymbol();
      
        # Remember order that symbols are listed in the summary table. But don't include duplicates!
        # Duplicates cause "duplicate target name" errors. 
        if symb not in symbols_ordered:
          symbols_ordered.append(symb); 

        symb = ":ref:`%s`" % symb;        

        row = [""]*4;
        row[0] = gene.chrom;
        row[1] = symb;
        
        # If "any" is in the search term list, replace that text with something
        # more descriptive. 
        search_terms = gene.getSearchTermWords();
        try:
          ind = search_terms.index("any");
          search_terms[ind] = "pubmed";
        except:
          pass;

        row[2] = ", ".join(search_terms);
        row[3] = str(gene.getTotalPubmed());
        user_t.add_row(row);
      
      print >> out, user_t.get_string(hrules=prettytable.ALL,rest=True)
      print >> out, "";
    
    # Print section for each gene. 
    print >> out, "Genes";
    print >> out, "-----";
    print >> out, "";
    
    omim_links = [];
    for symb in symbols_ordered:
      gene = Gene.valueOf(symb);
      symb = gene.getSymbol(); # override symbol with primary symbol
      
      print >> out, ".. _%s:" % symb;
      print >> out, "";
      
      print >> out, symb;
      print >> out, '^' * len(symb);
      print >> out, "";
      print >> out, "General information regarding this gene: ";
      print >> out, "";
      
      print >> out, ":Full gene name: " + str(gene.desc);
      print >> out, ":Entrez Gene ID: " + str(gene.uid);
      print >> out, ":Location: " + str(gene.loc);
      print >> out, ":Synonyms: " + ", ".join(gene.syns);
      print >> out, ":Type: " + str(gene.type);
      print >> out, "";
      
      # SNP table near gene. 
      if len(gene.snps) > 0:
        snp_table = prettytable.PrettyTable(['SNP','Distance (bp)','Direction']);
        print >> out, "SNPs given by the user that are near or inside this gene: ";
        print >> out, "";
        
        for snp in sorted(gene.snps,key=lambda x: int(gene.snps[x]['distance'])):
          snp_table.add_row([snp,gene.snps[snp]['distance'],gene.snps[snp]['updown']]);
        
        print >> out, snp_table.get_string(hrules=prettytable.ALL,rest=True);
        print >> out, "";
      
      # Region table. 
      if len(gene.regions) > 0:
        region_table = prettytable.PrettyTable(['Region','Containment']);
        print >> out, "Regions given by the user overlapping this gene:";
        print >> out, "";
        
        for region,contain in gene.regions.iteritems():
          region_table.add_row([str(region),contain]);
        
        print >> out, region_table.get_string(hrules=prettytable.ALL,rest=True);
        print >> out, "";
      
      # eQTLs. 
      if len(gene.eqtls) > 0:
        print >> out, "This gene is a target of the following eQTL associations:";
        print >> out, "";
        gene.eqtls.write_rest_table(out);
        print >> out, "";
      
      # Summary. 
      print >> out, ".. _%s %s:" % (gene.symb,Gene.SUMMARY);
      print >> out, "";
      print >> out, "**Summary**";
      print >> out, "";
      if gene.summary != None and gene.summary != "None":
        print >> out, markup_term(gene.summary,gene.getSearchTerms().get(Gene.SUMMARY));
      else:
        print >> out, "None available.";
      print >> out, "";
      
      # OMIM summary. 
      if gene.omim_id != None:
        print >> out, ".. _OMIM ID %s: %s" % (str(gene.omim_id),OMIM_WEB_URL + str(gene.omim_id));
        print >> out, "";
        
      print >> out, ".. _%s %s:" % (gene.symb,Gene.OMIM);  
      print >> out, "";
      
      print >> out, "**OMIM Summary**";
      print >> out, "";
      if gene.omim_text != None:
        omim_text = markup_term(gene.omim_text,gene.getSearchTerms().get(Gene.OMIM));
        print >> out, gene.omim_text + " [%s]" % "`OMIM ID %s`_" % str(gene.omim_id); 
      else:
        print >> out, "No summary available for this gene.";
      print >> out, "";
      
      # Phenotypes. 
      print >> out, ".. _%s %s:" % (gene.symb,Gene.PHENOTYPE);
      print >> out, "";
      
      print >> out, "**Phenotypes**";
      print >> out, "";
      if len(gene.phenotypes) > 0:
        for pheno in gene.phenotypes:
          pheno = markup_term(pheno,gene.getSearchTerms().get(Gene.PHENOTYPE));
          print >> out, "* %s" % str(pheno);
        print >> out, "";
      else:
        print >> out, "No phenotypes found linked to this gene.";
        print >> out, "";
        
      # GO Terms. 
      print >> out, ".. _%s %s:" % (gene.symb,Gene.GOTERM);
      print >> out, "";
      
      print >> out, "**Gene Ontology**";
      print >> out, "";
      if len(gene.GO) > 0:
        for goterm in gene.GO:
          goterm = markup_term(goterm,gene.getSearchTerms().get(Gene.GOTERM));
          print >> out, "* %s" % str(goterm);
        print >> out, "";
      else:
        print >> out, "No GO terms found linked to this gene.";
        print >> out, "";  
      
      # KEGG pathways. 
      print >> out, ".. _%s %s:" % (gene.symb,Gene.PATHWAY);
      print >> out, "";
      
      print >> out, "**Pathways**";
      print >> out, "";
      if len(gene.pathways) > 0:
        for pathway,url in gene.pathways.iteritems():
          pathway = markup_term(pathway,gene.getSearchTerms().get(Gene.PATHWAY));
          print >> out, "* `%s <%s>`_" % (str(pathway),str(url));
        print >> out, "";
      else:
        print >> out, "No pathways found linked to this gene.";
        print >> out, "";
      
      # GeneRIFs
      print >> out, ".. _%s %s:" % (gene.symb,Gene.GENERIF);
      print >> out, "";
      
      print >> out, "**GeneRIFs**";
      print >> out, "";
      pubmed_links = [];
      if len(gene.gene_rifs) > 0:
        for pmid in gene.gene_rifs:
          pubmed_links.append(".. _PMID %s: %s" % (str(pmid),PUBMED_WEB_URL + str(pmid)));
          for text in gene.gene_rifs[pmid]:
            text = markup_term(text,gene.getSearchTerms().get(Gene.GENERIF));
            print >> out, "* %s [%s]" % (text,"`PMID %s`_" % str(pmid));
        print >> out, "";
        for link in pubmed_links:
          print >> out, link;
        print >> out, "";
      else:
        print >> out, "None available.";
        print >> out, "";
      
      # PubMed Articles. 
      print >> out, ".. _%s %s:" % (gene.symb,Gene.PUBMED);
      print >> out, "";
      
      print >> out, "**PubMed Articles**";
      print >> out, "";
      
      if len(gene.pubmed_ids) > 0:
        print >> out, "*Recent articles:*"
        print >> out, "";
        Pubmed.writeReST(gene.pubmed_ids,search_terms,out);
      else:
        print >> out, "No recent articles (or PubMed lookup was not enabled.)";
        print >> out, "";
      
      for term in gene.pubmed_by_term:
        if term == "any":
          print >> out, "*Top Pubmed articles linked to gene %s matching any search term:* " % str(gene.symb);
          print >> out, "";
        else:
          print >> out, "*Top Pubmed articles linked to gene %s matching search term \"%s\":*" % (str(gene.symb),term);
          print >> out, "";
        
        Pubmed.writeReST(gene.pubmed_by_term[term],search_terms,out);
      
  # Load a set of eQTLs into the gene database.
  @staticmethod
  def updateEQTL(eqtl_set):
    for q in eqtl_set:
      gene = Gene.valueOf(q.gene);
      gene.eqtls.add(q);

  @staticmethod
  def updateUserGenes(genes):
    for g in genes:
      gene = Gene.valueOf(g);
      gene.user_requested = True;

  # Load GeneRIF information for specified genes. 
  # This is only loaded optionally.. it's a lot of junk. 
  @staticmethod
  def loadGeneRIF(symbols=set()):
    if not isIterable(symbols):
      raise Exception, "You must pass a list or set of symbols to loadGeneRIF.";
 
    for symb in symbols:
      gene = Gene.valueOf(symb);
      soup = gene.getSoup();

      if soup != None:
        gene.gene_rifs = extractGeneRIF( soup );

  @staticmethod
  def loadPositions(finder): 
    seen = {};
    for symb,gene in Gene._SYMB.iteritems():
      check = seen.get(gene);
      if check != None:
        continue;
      
      seen[gene] = 1;
      
      (chrom,start,stop) = finder.getPos(gene.symb);
      gene.chrom = chrom;
      gene.txStart = start;
      gene.txEnd = stop;
      
      chr = chrom2chr(chrom);
      gene.genome_region = ChromRegion(chr,start,stop);

  # Load OMIM information for specified genes. 
  # If the genes specified have not been populated previously, they will be after this call. 
  @staticmethod
  def loadOMIM(symbols=set()):
    if not isIterable(symbols):
      raise Exception, "You must pass a list or set of either symbols or uids into loadOMIM.";
    
    Gene.populate(symbols);

    omim_info = dict();
    for symb in symbols:
      gene = Gene.valueOf(symb);
      omim_id = gene.getOMIM();
      if omim_id != None:
        omim_info[omim_id] = symb;

    # Some genes may not have OMIM IDs. 
    # Remove the 'None' key. 
    if omim_info.has_key('None'):
      omim_info.pop('None');

    # Get the soup for all OMIM IDs. 
    soup = fetchOMIMSoup(omim_info.keys());

    # Process the soup and store info. 
    for omim_id in omim_info:
      root = findOMIMRoot(soup,omim_id);
      if root != None:
        gene = Gene.valueOf(omim_info[omim_id]);
        gene.setOMIMSoup(root);
        gene.setOMIMText( extractOMIMText(root) );

  # Load Pubmed information for specified genes. 
  # If the genes have not been populated previously, they will be after this call. 
  @staticmethod 
  def loadPubmed(symbols=set(),terms=set(),pnum=10,per_term=False):
    if not isIterable(symbols):
      raise ValueError, "Contact developer.";
  
    symbols = set(symbols);

    # Make sure they're all populated. 
    bad_symbs = Gene.populate(symbols);
    #symbols = symbols.difference(bad_symbs);

    # We need a UID --> symbol translation. 
    uid_symbol = {};
    for symb in symbols:
      uid = Gene.valueOf(symb).getUID();
      if uid != None:
        uid_symbol[uid] = symb;

    # Next, lookup articles for each UID. 
    pubmed_data = findPubmed( uid_symbol.keys() );

    # Store the number of articles found, and 
    # also store the top "pnum" articles. 
    pop_articles = set();
    for uid in pubmed_data:
      symb = uid_symbol[uid];
      gene = Gene.valueOf(symb);
      gene.setTotalPubmed( len(pubmed_data[uid]) );
      for i in xrange(len(pubmed_data[uid])):
        pmid = pubmed_data[uid][i];
        if i < pnum:
          pop_articles.add(pmid);
          gene.addPubmed(pmid);
        else:
          break;

    # Search Pubmed articles for each gene and set of search terms. 
    if len(terms) > 0:
      if per_term:
        num_genes = len(pubmed_data);
        num_terms = len(terms);
        
        # This could take a really long time.. are they sure they want to do this? 
        msg = "** WARNING: Searching Pubmed for each search term + gene combination with %s genes and %s search terms will require approximately %s minutes - do you wish to proceed [y/n]? " %\
        (str(num_genes),str(num_terms),str( 3*num_genes*num_terms/60.0 ));

        answer = confirm(msg);

        if answer:
          for uid in pubmed_data:
            gene = Gene.valueOf(uid_symbol[uid]);
            # Search using NCBI gene <--> pubmed links. 
            result = searchArticlesByTerm(pubmed_data[uid],terms);
            for term in result:
              for id in list(result[term])[0:pnum:1]:
                gene.addPubmedByTerm(id,term);
                pop_articles.add(id);

            # Search using Naive Pubmed search. 
            for term in terms:
              naive_articleids = naiveSearchArticlesByTerm(gene.symb,term,pnum);
              for id in list(naive_articleids)[0:pnum:1]:
                gene.addPubmedByTerm(id,term);
                pop_articles.add(id);
              
              # Wait 3 seconds. 
              sleep(3);
              
        else:
          Gene.loadPubmed(symbols,terms,pnum,False);
      else:
        for uid in pubmed_data:
          gene = Gene.valueOf(uid_symbol[uid]);

          # Search using NCBI gene <--> pubmed. 
          result = searchArticlesAnyTerm(pubmed_data[uid],terms);
          for id in list(result)[0:pnum:1]:
            gene.addPubmedByTerm(id,"any");
            pop_articles.add(id);

          # Search using naive approach. 
          naive_articleids = naiveSearchArticlesAnyTerm(gene.symb,terms,pnum);
          for id in list(naive_articleids)[0:pnum:1]:
            gene.addPubmedByTerm(id,"any");
            pop_articles.add(id);

          # Wait 3 seconds. 
          sleep(3);

    # Populate pubmed database of articles. 
    Pubmed.populate(pop_articles);

  # Flyweight method. Returns a gene if already seen before, 
  # otherwise returns a new gene instance. 
  @staticmethod
  def valueOf(symb=None,uid=None):
    if symb != None:
      gene = Gene._SYMB.get(symb);
    elif uid != None:
      gene = Gene._UID.get(uid);
    else:
      raise ValueError;

    if gene == None:
      gene = Gene();
      if symb != None:
        gene.setSymbol(symb);
      if uid != None:
        gene.setUID(uid);
        
    return gene;

  # Populates information for many genes at once. 
  # This is much more efficient than doing genes individually.
  # You can pass in a set of UIDs, or symbols. 
  # If you pass in a symbol or uid we've already seen, this method 
  # will simply do nothing for it. 
  # If you pass in a symbol or uid that cannot be found, 
  # it will be returned in a list. 
  @staticmethod
  def populate(symbols=set(),uids=set(),verbose=True):
    # Filter out any gene symbols we already know about. 
    symbol_list = filter(lambda x: not Gene._SYMB.has_key(x),list(symbols));

    # For gene symbols we don't know about, these must be new - need to 
    # get their UIDs so we can download information about them. 
    converted_uids = symb2uid(symbol_list);

    # Combine list of UIDs for gene symbols with UIDs passed into function. 
    all_uids = list(uids) + converted_uids;

    # Make sure none of the UIDs we've gathered so far have already been seen. 
    final_uids = filter(lambda x: not Gene._UID.has_key(x),all_uids);

    # Finally, make sure there's no duplicate UIDs. 
    final_uids = list(set(final_uids));

    # If there are no UIDs remaining at this point, a query is not needed.
    if len(final_uids) == 0:
      return [];

    # Ask NCBI for XML data for all uids.
    try:
      if _SNIPPER_DEBUG:
        print >> sys.stderr, "DEBUG: Gene.populate() is connecting to NCBI..";
        print >> sys.stderr, "DEBUG: UIDs were: %s" % str(final_uids);
    except:
      pass
    
    soup = fetchGeneSoup(final_uids);
    for uid in final_uids:
      root = findUIDRoot(soup,uid);

      # Has this gene record been discontinued? 
      if isDiscontinued(root):
        continue;

      symbol = extractSymb(root);
      syns = extractSyns(root);

      # If the symbol is null, this is a weird gene. 
      # I don't know why NCBI would populate them without a symbol, but they do. 
      # So, we'll just set the symbol to be whatever was passed in. 
      if symbol == 'None' or symbol == None:
        symbol = findCorrectSymbol(symbol_list,root);

      if verbose:
        print >> sys.stderr, "Parsing gene %s.." % symbol;

      # Grab information for the gene. 
      gene = Gene.valueOf(symbol);
      gene.setUID( uid );
      gene.setSymbol( symbol );
      gene.setSyns( syns );
      gene.setType( extractType(root) );
      gene.setLoc( extractLocus(root) );
      gene.setSummary( extractSummary(root) );
      gene.setGO( extractGO(root) );
      gene.setPhenotypes( extractPhenotypes(root) );
      gene.setOMIM( extractOMIM(root) );
      gene.setDesc( extractDesc(root) );
      gene.setKEGG( extractKEGG(root) );
      gene.setSoup( root );

      # Insert into "database."
      Gene._UID[uid] = gene;
      Gene._SYMB[symbol] = gene;
      for syn in syns:
        Gene._SYMB[syn] = gene;

    # Return a list of symbols that couldn't be found. 
    bad_symbs = list();
    for symb in symbol_list:
      if not Gene._SYMB.has_key(symb):
        bad_symbs.append(symb);
        print >> sys.stderr, "Warning: no information found for gene %s.." % str(symb);
        
    return bad_symbs;
