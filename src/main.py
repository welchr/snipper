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
import pdb
import time
import subprocess
from gene import *
from util import *
from constants import *
from settings import *
from genefinder import *
from scandb import *
from mimi import *
from intdb import *
from ranking import *
from snp import *
from textwrap import fill
from ncbi import *

try:
  import sphinx
except:
  pass

PROG_VERSION = "1.2";
PROG_VERSION_STRING = "Version %s   " % PROG_VERSION;
PROG_DATE = "03-16-2011";
NUM_GENES_WARN = 100;

def get_core_path():
  snipper_src = os.path.dirname(os.path.abspath(sys.argv[0]));
  core_path = os.path.abspath(os.path.join(snipper_src,"../conf/sphinx"));
  
  return core_path;

def print_program_header():
  print >> sys.stderr, "+------------------------------------------------+";
  print >> sys.stderr, "|  Snipper   |   %s  |   %s  |" % (PROG_VERSION_STRING,PROG_DATE);
  print >> sys.stderr, "+------------------------------------------------+";
  print >> sys.stderr, "|     A utility for locating genes near SNPs     |";
  print >> sys.stderr, "|                                                |";
  print >> sys.stderr, "|      Author: Ryan Welch (welchr@umich.edu)     |";
  print >> sys.stderr, "+------------------------------------------------+";
  print >> sys.stderr, "";

def print_genes_per_snp(out,gene_symbols):
  snp_gene = {};
  for symb in gene_symbols:
    gene = Gene.valueOf(symb);
    syns = gene.syns;
    for snp in gene.snps.keys():
      hashed = "/".join([symb] + syns);
      snp_gene.setdefault(snp,[]).append(hashed);

  format_string = "%-12s %s";
  
  print >> out, format_string % ('SNP','Gene/Aliases');
  print >> out, format_string % ('---','------------');
  for snp,genelist in snp_gene.iteritems():
    for genes in genelist:
      print >> out, format_string % (snp,genes);

  print >> out, "";

def print_summary(settings,gene_symbols):
  out = settings.output_file;

  # Create the format string. 
  format_string = "%-10s %-9s";

  # Were search terms specified? 
  found_terms = False;
  if len(settings.terms) > 0:
    found_terms = True;
    format_string += " %-9s";

  # Was Pubmed specified? 
  if settings.pubmed:
    format_string += " %-9s";

  # Print header row. 
  header_list = ["Gene","# SNPs"];
  if found_terms:
    header_list.append("# Terms");
  if settings.pubmed:
    header_list.append("Total Pubmed");
  print >> out, format_string % tuple(header_list);
  dash_list = [];
  for i in xrange(len(header_list)):
    element = "";
    for j in xrange(len(header_list[i])):
      element += "-";
    dash_list.append(element);
  print >> out, format_string % tuple(dash_list);

  for symb in gene_symbols:
    gene = Gene.valueOf(symb);
    symb = gene.getSymbol();
    num_snps = gene.countAllSNP();
    num_terms = gene.countTerms();
    num_total_pubmed = gene.getTotalPubmed();

    output_list = [symb,str(num_snps)];
    if found_terms:
      output_list.append(str(num_terms));
    if settings.pubmed:
      output_list.append(str(num_total_pubmed));

    print >> out, format_string % tuple(output_list);

# Originally from: http://code.activestate.com/recipes/576694/
# Removes duplicates from iterables and maintains order, where order existed
# (sets do not have an ordering, so they are accessed in undefined order.) 
# Modifications for arbitrary number of iterables. 
def remove_duplicates(*args):
  mmap = {};
  oset = [];

  for iterable in args:
    for item in iterable:
      if item not in mmap:
        mmap[item] = 1;
        oset.append(item);

  return oset;

def get_current_snps(snp_list,db_file):
  db = SNPDB(db_file);
  new_list = [];
  d = {};
  for snp in snp_list:
    trans = db.get_current_name(snp);
    if trans != snp:
      print >> sys.stderr, "Warning: SNP %s was merged into %s in the current genome build.." % (snp,trans);
      d[snp] = trans;
    
    new_list.append(trans);
  
  return (new_list,d);

def run_snipper(settings):
  output_file = settings.output_file;
  user_genes = settings.genes;
  
  # Build. Make sure the user knows it. 
  print >> sys.stderr, "Using build %s for SNP and gene positions.." % settings.build;
  
  # Fix user's SNPs to map to current genome build. 
  (settings.snpset,settings.snp_trans) = get_current_snps(settings.snpset,settings.db_file);
  
  # Feedback about user's provided genes. 
  if len(user_genes) > 0:
    print >> sys.stderr, "User provided %i gene(s) for inclusion.." % len(user_genes);
    for gene in user_genes:
      print >> sys.stderr, ".. %s" % str(gene);
  
  # Find the genes lying near the SNPs. 
  finder = SQLiteGeneFinder(settings.db_file);
  print >> sys.stderr, "Locating genes near provided SNPs..";
  
  snps_normal = [];
  snps_1000G = [];
  for snp in settings.snpset:
    if "chr" in snp:
      snps_1000G.append(snp);
    else:
      snps_normal.append(snp);
  
  finder.find(snps_normal,settings.distance,settings.num_genes);
  finder.find1000G(snps_1000G,settings.distance,settings.num_genes);
  
  print >> sys.stderr, "Finding genes within provided chromosomal regions..";
  finder.findRegions(settings.regions);
  
  finder_genes = finder.getGenesByNearest();
  print >> sys.stderr, ".. found %i gene(s)!" % len(finder_genes);

  # If the number of genes is rather large, warn the user that this might be a 
  # bit too large of a task. 
  num_genes = len(user_genes) + len(finder_genes);
  if num_genes > NUM_GENES_WARN:
    print >> sys.stderr, "";
    print >> sys.stderr, fill("Warning: the number of genes in total (genes near SNPs + "
                          "genes in regions + genes requested by the user) is very "
                          "large (> %i). Snipper may not be able to handle a set of this size. "
                          "If you continue, please note the runtime will be very high. " % NUM_GENES_WARN);
    print >> sys.stderr, "";

  # If eQTL settings were enabled, find those genes here.
  scandb_results = set();
  if settings.scandb and len(settings.snpset) > 0:
    print >> sys.stderr, "Searching ScanDB for eQTLs..";
    scandb_results = scandb_snps(settings.snpset,settings.scandb_pval);
    
    # Update user about what happened.
    if scandb_results != None and len(scandb_results) > 0:
      print >> sys.stderr, ".. found %i genes associated with %i input SNPs!" % (scandb_results.numGenes(),scandb_results.numSNP());

  # Figure out list of gene symbols to ask NCBI about.
  scandb_genes = [e.gene for e in sorted(scandb_results,key=lambda x: x.pval)];
  
  # Final list of genes that will be used.
  # They are ordered in terms of "importance." 
  gene_symbols = remove_duplicates(user_genes,scandb_genes,finder_genes);
  
  print >> sys.stderr, "Downloading information from NCBI on %i genes.." % len(gene_symbols);
  Gene.populate(gene_symbols);
  
  # Add SNP position and eQTL info to database, after symbols have been resolved to gene objects.
  finder.loadGeneDB();
  Gene.updateEQTL(scandb_results);
  Gene.updateUserGenes(user_genes);
  Gene.loadPositions(finder);
 
  # Should we look for interactions between genes? 
  direct_ints = None;
  if settings.mimi:
    print >> sys.stderr, "Querying MiMI for direct interactions between genes near SNPs..";
    direct_ints = find_mimi_interactions(gene_symbols);
    
    if direct_ints != None:
      if direct_ints.has_gene_interactions():
        print >> sys.stderr, ".. found interactions!";
      else:
        print >> sys.stderr, ".. no interactions found";
  
  # If OMIM is specified, load OMIM data. 
  if settings.omim:
    print >> sys.stderr, "Loading OMIM information..";
    Gene.loadOMIM(gene_symbols);  

  # If Pubmed specified, load Pubmed information. 
  if settings.pubmed:
    print >> sys.stderr, "Loading Pubmed information..";
    Gene.loadPubmed(gene_symbols,settings.terms,settings.pnum,settings.per_term);

  # If GeneRIF..
  if settings.gene_rif:
    print >> sys.stderr, "Loading GeneRIFs..";
    Gene.loadGeneRIF(gene_symbols);

  # Now, we need to check search terms. Each gene has a soup we can search in. 
  if len(settings.terms) > 0:
    print >> sys.stderr, "Scanning for search terms in gene information..";
    for symb in gene_symbols:
      gene = Gene.valueOf(symb);
      for term in settings.terms:
        gene.search(term);

  # Rank genes. 
#  ranker = RankH1();
#  ranker.setInteractionDB(direct_ints);
#  gene_symbols = ranker.rank_genes(gene_symbols);

  # If we're in HTML mode, write that version of the report. 
  if settings.html:
    make_rest(
      get_core_path(),
      settings,
      gene_symbols,
      scandb_results,
      direct_ints
    );
  else:
    # Console output. 
    # Print our summary table. 
    if output_file == sys.stdout: 
      print >> sys.stderr, "";
      print_sep(sys.stderr);
      print "";
      print >> sys.stderr, "Summary information: \n";
    print_summary(settings,gene_symbols);
    
    # Print gene symbols on per SNP basis. 
    if output_file == sys.stdout:
      print >> sys.stderr, "";
      print_sep(sys.stderr);
      print "";
      print >> sys.stderr, "Genes and their aliases near each SNP: ";
      print "";
    else:
      print >> output_file, "";
    print_genes_per_snp(output_file,gene_symbols);
  
    # Print direct interactions between genes. 
    if direct_ints != None and direct_ints.has_gene_interactions():
      if output_file == sys.stdout:
        print >> sys.stderr, "";
        print_sep(sys.stderr);
        print "";
        print >> sys.stderr, "Direct interactions detected between genes near any given SNP:";
        print "";
      else:
        print >> output_file, "";
        
      print_interactions(direct_ints,output_file);
      
      if len(settings.terms) > 0:
        run_search_interactions(direct_ints,settings.terms,output_file);
  
    # Finally, print the information.
    if output_file == sys.stdout: 
      pass
      #print >> sys.stderr, "";
    Gene.writeAll(gene_symbols,output_file);

def make_rest(core_path,settings,gene_symbols,scandb_results,direct_ints):
  # Directory path to write ReST and HTML files. 
  dir_path = settings.getRestDir();
  
  # Create directory for rst files. 
  rst_path = os.path.join(dir_path,"rst");
  html_path = os.path.join(dir_path,"html");
  os.mkdir(rst_path);
  os.mkdir(html_path);
  
  # Write ReST formatted output. 
  index = os.path.join(rst_path,"index.rst");
  with open(index,'w') as f:
    write_rest_index(f);
  
  input = os.path.join(rst_path,"input.rst");
  with open(input,'w') as f:
    write_rest_input(settings,f);

  genes = os.path.join(rst_path,"genes.rst");
  with open(genes,'w') as f:
    Gene.writeReST(gene_symbols,settings.db_file,settings.terms,f);
  
  terms = os.path.join(rst_path,"terms.rst");
  with open(terms,'w') as f:
    write_rest_terms(gene_symbols,f);
  
  eqtls = os.path.join(rst_path,"eqtls.rst");
  with open(eqtls,'w') as f:
    write_rest_eqtls(scandb_results,f);
  
  inters = os.path.join(rst_path,"interactions.rst");
  with open(inters,'w') as f:
    write_rest_ints(direct_ints,f);
    
  # Run sphinx build process. 
  try:
    from sphinx.cmdline import main as sphinx_build
  except:
    raise Exception("Error: sphinx is not installed, did you run setup_snipper.py?");

# Code does not work in Windows
#  args = "%s -b html -c %s %s %s" % (
#    core_path,
#    core_path,
#    rst_path,
#    html_path
#  );
  
  args = [
    core_path,
    "-b",
    "html",
    "-c",
    core_path,
    rst_path,
    html_path
  ];
  
  # sphinx_build(shlex.split(args)); # doesn't work in Windows
  print "Creating HTML report using Sphinx..";
  mute_std();
  
  try:
    if _SNIPPER_DEBUG: 
      unmute_std();
  except:
    pass
  
  sphinx_build(args);
  unmute_std();
  
  # Write README file. 
  write_readme(os.path.join(dir_path,"README.txt"));
  
  print "Wrote HTML report to: %s" % html_path;
  
  report_index = os.path.join(html_path,"index.html");
  print "Please open %s to view the report." % report_index;

def write_rest_input(settings,out=sys.stdout):
  print >> out, "User Input and Settings";
  print >> out, "=======================";
  print >> out, "";
  
  print >> out, "The following sections detail the input to Snipper, "\
                "as well as the settings supplied.";
  print >> out, "";
  
  print >> out, "SNPs";
  print >> out, "----";
  print >> out, "";
  
  if len(settings.snp_trans) > 0:
    print >> out, ("The following SNPs provided by the user were merged into a "
                   "new SNP name in the current genome build: ");
    print >> out, "";
    
    table = PrettyTable(['Original SNP Name','Current SNP Name']);
    for snp,trans in settings.snp_trans.iteritems():
      table.add_row([snp,trans]);
    
    print >> out, table.get_string(hrules=prettytable.ALL,rest=True);
    print >> out, "";
    print >> out, ("For the remainder of the results, we will refer to the SNPs "
                   "above by the current genome build name. ");
    print >> out, "";
  
  snps = settings.snpset;
  snp_db = SNPDB(settings.db_file);
  if len(snps) > 0:
    table = PrettyTable(['SNP','Chromosome','Position']);
    for snp in snps:
      chrpos = snp_db.get_pos(snp);
      if chrpos != None:
        table.add_row([snp,chrpos[0],chrpos[1]]);
      else:
        table.add_row([snp,'Unknown','Unknown']);
    
    print >> out, table.get_string(hrules=prettytable.ALL,rest=True);
    print >> out, "";
  else:
    print >> out, "No SNPs given by the user.";
    print >> out, "";
  
  regions = settings.regions;
  print >> out, "Regions";
  print >> out, "-------";
  print >> out, "";
  
  if len(regions) > 0:
    region_table = PrettyTable(['Chromosome','Start','End']);
    for region in regions:
      region_table.add_row([
        region.chr,
        region.start,
        region.end
      ]);
    
    print >> out, region_table.get_string(hrules=prettytable.ALL,rest=True);
    print >> out, "";
  else:
    print >> out, "No regions given by the user.";
    print >> out, "";
  
  print >> out, "Genes";
  print >> out, "-----";
  print >> out, "";
  
  genes_ordered = sorted([Gene.valueOf(g) for g in settings.genes],key = lambda x: x.genome_region);
  if len(settings.genes) > 0:
    gene_table = PrettyTable(['Gene','Chrom','txStart','txEnd']);
    for g in genes_ordered:
      gene_table.add_row([
        g.symb,
        g.chrom,
        g.txStart,
        g.txEnd
      ]);
    
    print >> out, gene_table.get_string(hrules=prettytable.ALL,rest=True);
    print >> out, "";
  else:
    print >> out, "No genes explicitly requested by the user.";
    print >> out, "";
  
  print >> out, "Parameters";
  print >> out, "----------";
  print >> out, "";
  
  st = PrettyTable(['Parameter','Value']);
  
  if len(snps) > 0:
    st.add_row(["List of SNPs",", ".join(snps)]);
    st.add_row(["Distance to search from each SNP",str(settings.distance)]);
    st.add_row(["Number of genes to return near each SNP",str(settings.num_genes)]);

  if len(settings.genes) > 0:
    st.add_row(['List of genes',", ".join([i for i in settings.genes])]);
  
  if len(settings.regions) > 0:
    st.add_row(['List of regions',", ".join([str(i) for i in settings.regions])]);
    
  if len(settings.terms) > 0:
    st.add_row(['Search terms',", ".join([i for i in settings.terms])]);
  else:
    st.add_row(['Search terms','None provided.']);
    
  st.add_row(["OMIM",str(settings.omim)]);
  
  st.add_row(["PubMed",str(settings.pubmed)]);
  if settings.pubmed:
    st.add_row(["Maximum number of recent PubMed articles to return",str(settings.pnum)]);
    st.add_row(["Search for each (search term)*(gene) combination",str(settings.per_term)]);
  
  st.add_row(["GeneRIF",str(settings.gene_rif)]);
  
  st.add_row(["ScanDB",str(settings.scandb)]);
  if settings.scandb:
    st.add_row(["ScanDB P-value threshold",str(settings.scandb_pval)]);
  
  st.add_row(["MiMI",str(settings.mimi)]);
  
  st.add_row(["Build for SNP and gene positions",str(settings.build)]);
  
  print >> out, st.get_string(hrules=prettytable.ALL,rest=True);
  print >> out, "";

def write_rest_eqtls(eqtls,out=sys.stdout):
  print >> out, "Expression QTLs (eQTLs)";
  print >> out, "=======================";
  print >> out, "";

  print >> out, "Each SNP provided by the user is looked up in the SCAN "\
                "database for eQTL associations. More information on SCAN "\
                "can be found below.";
  print >> out, "";
  
  print >> out, ".. seealso::"
  print >> out, "";
  print >> out, "  SCAN Website";
  print >> out, "    %s" % SCANDB_WEB_URL;
  print >> out, "  SCAN Publication";
  print >> out, "    %s" % SCANDB_PAPER_URL;
  print >> out, "";
  
  if len(eqtls) > 0:
    print >> out, "The following SNPs were found to be eQTLs by querying SCAN:";
    print >> out, "";
    
    eqtls.write_rest_table(out);
    print >> out, "";
  else:
    print >> out, "None of the input SNPs were detected as eQTLs.";
    print >> out, "";
  
def write_rest_index(out=sys.stdout):
  index = """
Snipper Report
==============

Contents:

.. toctree::
   :maxdepth: 1

   input
   eqtls
   genes
   interactions
   terms

Created |today|.

.. seealso::

  Snipper Website
    http://csg.sph.umich.edu/boehnke/snipper/

  Snipper Documentation
    http://csg.sph.umich.edu/boehnke/snipper/snipper_docs.pdf
""";
  
  print >> out, index;

def get_command_line():
  args = []
  for i in sys.argv[1:]:
    if i[0] == "-":
      args.append(i);
      continue;

    args.append('"' + str(i) + '"');
  
  return "snipper.py %s" % " ".join(args);

def write_readme(path):
  f = open(path,'w');
  text = """
***
To view the Snipper report, please open the "index.html" file under the html folder. 
***

html/
-----
This directory contains the Snipper HTML output. Open "index.html" to view the report. 

rst/
----
This directory contains the source files used to generate the HTML output. These files may be useful to users familiar with the ReST markup language.

Command line
------------
Snipper was executed with the following command: %s

""" % get_command_line();

  print >> f, text.strip();
  f.close();

def write_rest_terms(genes,out=sys.stdout):
  term_links = {};
  for gene in genes:
    g = Gene.valueOf(gene);
    for where,terms in g.getSearchTerms().iteritems():
      for term in terms:
        link = "%s %s" % (g.symb,where);
        link_name = "%s - %s" % (g.symb,where);
        term_links.setdefault(term,[]).append(
          ":ref:`%s <%s>`" % (link_name,link)
        );
  
  print >> out, "Search Terms"
  print >> out, "============"
  print >> out, ""
  
  if len(term_links) > 0:
    sorted_terms = sorted(term_links);
    for term in sorted_terms:
      links = term_links[term];
      if term == "any":
        term = "any search term";
      print >> out, "**%s**" % term;
      print >> out, "";
      for link in links:
        print >> out, "* %s" % link;
      print >> out, "";
  else:
    print >> out, "No search terms were given by the user. Use --terms to supply them.";
    print >> out, "";

def find_mimi_interactions(genes):
  if not hasattr(genes,'__iter__'):
    raise ValueError;
  
  uids = set();
  for gene in genes:
    uids.add( str(Gene.valueOf(gene).getUID()) );
  
  db = InteractionDB();
  try:
    for uid in uids:
      ints = mimi_fetch_interactions(uid,taxid=TAXON_ID_HUMAN);
      db.add_gene_interactions(ints);
      time.sleep(1.25);
  except:
    print >> sys.stderr, "Error: could not connect to MiMI. The message was:";
    print >> sys.stderr, str(sys.exc_info()[1]);
    return None;
    
  direct = db.direct_gene_interactions(uids);
  
  for inter in direct.iter_gene_interactions():
    sym1 = Gene.valueOf(uid=inter.gene1).getSymbol();
    sym2 = Gene.valueOf(uid=inter.gene2).getSymbol();
    inter.gene_sym1 = sym1;
    inter.gene_sym2 = sym2;
  
  return direct;

def write_rest_ints(inters,out=sys.stdout):
  print >> out, "Gene-Gene Interactions";
  print >> out, "======================";
  print >> out, "";
  
  print >> out, "The following sections list different sets of interactions "\
                "between genes found near input SNPs/regions. These "\
                "interactions are collected from NCIBI's MiMI database.";
  print >> out, "";
  print >> out, ".. seealso::"
  print >> out, "";
  print >> out, "  Michigan Molecular Interactions (MiMI)";
  print >> out, "    http://mimi.ncibi.org/";
  print >> out, "";
  
  print >> out, "Direct Interactions";
  print >> out, "-------------------";
  print >> out, "";
  
  if inters != None and inters.has_gene_interactions():
    inters.print_rest(out);
  else:
    print >> out, "No direct interactions found.";
    print >> out, "";
  
def print_interactions(inters,out=sys.stdout):
  for i in inters.iter_gene_interactions():
    print >> out, str(i);

def run_search_interactions(inters,terms,out=sys.stdout):
  results = dict();
  
  for i in inters.iter_gene_interactions():
    for term in terms:
      results.update(i.search(term));
  
  if len(results) > 0:
    print >> out, "Search terms that matched within direct gene interactions:";
    for key,value in results.iteritems():
      term_str = ", ".join(value);
      print >> out, "-- %s: %s" % (key,term_str);
  
  print >> out, "";
  return results;

def print_sep(out=sys.stdout):
  term_width = terminal_size()[0];
  if out not in (sys.stdout,sys.stderr): 
    term_width = 80;
      
  for i in xrange(term_width):
    out.write('=');
  out.write("\n");

#def run_gene(settings):
#  print >> sys.stderr, "Looking up information for genes @ NCBI..";
#  symbols = settings.genes;
#  terms = settings.terms;
#  pnum = settings.pnum;
#  per_term = settings.per_term;
#  bad_symbs = Gene.populate(symbols);
#
#  # Remove symbols that couldn't be found at NCBI. 
#  for symb in bad_symbs:
#    print >> sys.stderr, "Warning: no information found for gene " + str(symb);
#    symbols.remove(symb);
#
#  terms = settings.terms;
#  if settings.omim:
#    print >> sys.stderr, "Loading OMIM information..";
#    Gene.loadOMIM( settings.genes );
#
#  if settings.pubmed:
#    print >> sys.stderr, "Loading Pubmed information..";
#    Gene.loadPubmed( settings.genes,terms,pnum,per_term );
#
#  if settings.gene_rif:
#    print >> sys.stderr, "Loading GeneRIFs..";
#    Gene.loadGeneRIF( settings.genes );
#
#  for symb in symbols:
#    gene = Gene.valueOf(symb);
#    for term in terms:
#      gene.search(term);
#
#  Gene.writeAll(symbols,settings.output_file);

# Main entry point. 
def main(): 
  settings = Settings();
  
  print_program_header();

  # Were there arguments? 
  if len(sys.argv) <= 1:
    die("Nothing to do.. try -h for help.");

  # Check to see if services are online and available. 
  sys.stderr.write("Checking NCBI status.. ");
  ncbi_status = check_ncbi_status();
  if not ncbi_status[0]:
    sys.stderr.write("[FAIL]" + os.linesep);
    print >> sys.stderr, "Error: could not reach NCBI."
    print >> sys.stderr, fill("Snipper requires the ability to query NCBI in order to "
                          "function. It appears NCBI cannot be reached currently. This "
                          "could be caused by your firewall settings, or NCBI may be "
                          "offline temporarily. ");
    sys.exit(1);
  else:
    sys.stderr.write("[OK]" + os.linesep);

  run_snipper(settings);

if __name__ == "__main__":
  main();
