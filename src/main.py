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
import traceback
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
from multiprocessing import freeze_support
import __builtin__

# Global debug flag. 
__builtin__._SNIPPER_DEBUG = False;
__builtin__._SNIPPER_DEV = False;

PROG_VERSION = "1.3";
PROG_VERSION_STRING = "Version %s   " % PROG_VERSION;
PROG_DATE = "03-07-2012";
NUM_GENES_WARN = 50;

def get_core_path():
  snipper_src = os.path.dirname(os.path.abspath(sys.argv[0]));
  core_path = os.path.abspath(os.path.join(snipper_src,"../conf/sphinx"));
  
  return core_path;

def print_program_header():
  print >> sys.stderr, "+------------------------------------------------+";
  print >> sys.stderr, "|  Snipper   |   %s  |   %s  |" % (PROG_VERSION_STRING,PROG_DATE);
  print >> sys.stderr, "+------------------------------------------------+";
  print >> sys.stderr, "|     Contact: Ryan Welch (welchr@umich.edu)     |";
  print >> sys.stderr, "|     Web: csg.sph.umich.edu/boehnke/snipper/    |";
  print >> sys.stderr, "|     Github: github.com/welchr/Snipper/         |";
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
  
  if len(snp_gene) > 0:
    print >> out, format_string % ('SNP','Gene/Aliases');
    print >> out, format_string % ('---','------------');
    for snp,genelist in snp_gene.iteritems():
      for genes in genelist:
        print >> out, format_string % (snp,genes);
  
    print >> out, "";

def print_summary(settings,out,gene_symbols):
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
  new_set = set();
  d = {};
  for snp in snp_list:
    trans = db.get_current_name(snp);
    if trans != snp:
      print >> sys.stderr, "Warning: SNP %s was merged into %s in the current genome build.." % (snp,trans);
      d[snp] = trans;
    
    new_set.add(trans);
  
  return (new_set,d);

def print_debug():
  import traceback
  import pdb
  import inspect
  import platform
  import sys
  import os
  import copy

  def print_err(x):
    print >> sys.stderr, x;

  print_err("An error has occurred. Please provide the developer with the following information.");
  print_err("");
  print_err("Python executable: %s" % sys.executable);
  print_err("Platform: %s" % str(platform.uname()));
  print_err("Linux: %s" % str(platform.linux_distribution()));
  print_err("Windows: %s" % str(platform.win32_ver()));
  print_err("Mac: %s" % str(platform.mac_ver()));
  print_err("Library paths: ");
  for path in sys.path:
    print_err(".. %s" % path);
  print_err("Executing in: %s" % str(os.getcwd()));
  print_err("Called by: %s" % str(sys.argv[0]));
  print_err("");
  traceback.print_exc(file=sys.stderr);
  print_err("");
  for t in inspect.trace():
    print_err("Frame: %s" % str(" ".join(map(str,t[1:]))));
    for local,value in dict(t[0].f_locals).iteritems():
      value = str(value).replace(os.linesep,"");
      print_err("  %s : %s ..." % (local,value[0:120]));

def run_snipper(settings):
  outdir = settings.outdir;
  user_genes = settings.genes;
  
  check_ncbi();
  
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
  
  if len(settings.snpset) > 0:
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
  
  if len(settings.regions) > 0:
    print >> sys.stderr, "Finding genes within provided chromosomal regions..";
    finder.findRegions(settings.regions);
  
  finder_genes = finder.getGenesByNearest();
  if len(settings.snpset) > 0 or len(settings.regions) > 0:
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
  
  print >> sys.stderr, "Downloading gene information from NCBI..";
  Gene.populate(gene_symbols,generif=settings.gene_rif);
  
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
    print >> sys.stderr, "Querying omim.org for OMIM information..";
    Gene.loadOMIM(gene_symbols);  

  # If Pubmed specified, load Pubmed information. 
  if settings.pubmed:
    print >> sys.stderr, "Querying NCBI for Pubmed information..";
    Gene.loadPubmed(gene_symbols,settings.terms,settings.pnum,settings.per_term);

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

  # Write to console, or write to disk? 
  if len(sys.argv) <= 1 or not settings.console:
    print "Creating HTML report..";
    make_rest(
      get_core_path(),
      settings,
      gene_symbols,
      scandb_results,
      direct_ints
    );
    
    make_text(
      settings,
      gene_symbols,
      scandb_results,
      direct_ints
    );
  else:
    make_console(
      settings,
      gene_symbols,
      scandb_results,
      direct_ints
    );

def make_rest(core_path,settings,gene_symbols,scandb_results,direct_ints):
  # Directory path to write ReST and HTML files. 
  dir_path = settings.outdir;
  
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
  Gene.writeReST(gene_symbols,settings.db_file,settings.terms,genes);
  
  terms = os.path.join(rst_path,"terms.rst");
  with open(terms,'w') as f:
    write_rest_terms(gene_symbols,f);
  
  eqtls = os.path.join(rst_path,"eqtls.rst");
  with open(eqtls,'w') as f:
    write_rest_eqtls(scandb_results,f);
  
  inters = os.path.join(rst_path,"interactions.rst");
  with open(inters,'w') as f:
    write_rest_ints(direct_ints,f);
  
  # Write README file. 
  write_readme(os.path.join(dir_path,"README.txt"));
  
  # Run sphinx build process. 
  try:
    from sphinx.cmdline import main as sphinx_build
  except:
    raise Exception("Error: sphinx is not installed.");
  
  args = [
    core_path,
    "-Q",
    "-b",
    "html",
    "-c",
    core_path,
    rst_path,
    html_path
  ];
  
  try:
    sphinx_build(args);

    report_index = os.path.join(html_path,"index.html");
    print "HTML report created: %s." % report_index;
  except:
    raise;

def make_text(settings,gene_symbols,scandb_results,direct_ints):
  dir_path = settings.outdir;
  
  # Create directory for rst files. 
  text_path = os.path.join(dir_path,"text");
  os.mkdir(text_path);
  
  out = open(os.path.join(text_path,"snipper_results.txt"),'w');
  
  # Print our summary table. 
  print_summary(settings,out,gene_symbols);
  
  # Print gene symbols on per SNP basis. 
  print >> out, "";
  print_genes_per_snp(out,gene_symbols);

  # Print direct interactions between genes. 
  if direct_ints != None and direct_ints.has_gene_interactions():
    print >> out, "";  
    print_interactions(out,direct_ints);
    
    if len(settings.terms) > 0:
      run_search_interactions(direct_ints,settings.terms,out);

  # Finally, print the information.
  Gene.writeAll(gene_symbols,out);
  
def make_console(settings,gene_symbols,scandb_results,direct_ints):
  # Print our summary table.  
  print >> sys.stderr, "";
  print_sep(sys.stderr);
  print "";
  print >> sys.stderr, "Summary information: \n";
  print_summary(settings,sys.stdout,gene_symbols);
  
  # Print gene symbols on per SNP basis. 
  print >> sys.stderr, "";
  print_sep(sys.stderr);
  print "";
  print >> sys.stderr, "Genes and their aliases near each SNP: ";
  print "";
  print_genes_per_snp(sys.stdout,gene_symbols);

  # Print direct interactions between genes. 
  if direct_ints != None and direct_ints.has_gene_interactions():
    print >> sys.stderr, "";
    print_sep(sys.stderr);
    print "";
    print >> sys.stderr, "Direct interactions detected between genes near any given SNP:";
    print "";
      
    print_interactions(sys.stdout,direct_ints);
    
    if len(settings.terms) > 0:
      run_search_interactions(direct_ints,settings.terms,sys.stdout);

  # Finally, print the information.
  Gene.writeAll(gene_symbols,sys.stdout);

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
  
  # Figure out which genes provided by the user had positions. 
  # If it doesn't have a position, we won't list it here. 
  user_genes = [];
  user_genes_nopos = [];
  for g in settings.genes:
    gene = Gene.valueOf(g);
    if hasattr(gene.genome_region,'chr'):
      user_genes.append(g);
    else:
      user_genes_nopos.append(g);
  
  genes_ordered = sorted([Gene.valueOf(g) for g in user_genes],key = lambda x: x.genome_region);
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
    
    if len(user_genes_nopos) > 0:
      print >> out, "We could not find positions for the following user-requested genes:";
      print >> out, ""
      for g in user_genes_nopos:
        print >> out, "* %s" % g;
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

text/
-----
This directory contains a textual output file (also what would have been printed had snipper been run with --console.) This file could be useful for parsing information, but is not as useful for reading as the html report. 

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
  
def print_interactions(out,inters):
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

def print_sep(out):
  term_width = terminal_size()[0];
  if out not in (sys.stdout,sys.stderr): 
    term_width = 80;
      
  for i in xrange(term_width):
    out.write('=');
  out.write("\n");

def check_ncbi():
  # Check to see if services are online and available. 
  sys.stderr.write("Checking NCBI status.. ");
  ncbi_status = check_ncbi_status();
  if not ncbi_status[0]:
    sys.stderr.write("[FAIL]" + os.linesep);
    msg = ("Error: could not reach NCBI. Snipper requires the ability to query NCBI in order to "
          "function. It appears NCBI cannot be reached currently. This "
          "could be caused by your firewall settings, or NCBI may be "
          "offline temporarily. ");
    raise Exception, msg;
  else:
    sys.stderr.write("[OK]" + os.linesep);

# Main entry point. 
def main(): 
  settings = Settings();

  # If there are no arguments, try to launch the UI. 
  if len(sys.argv) <= 1:
    try:
      import ui
      ui.main();
    except:
      print >> sys.stderr, fill("I tried to start the GUI since no arguments were given, but "
                        "it failed (you may not have a python interpreter compiled "
                        "against Tk.) Check the documentation for more info on this. "
                        "In the meantime, you can run Snipper in command-line mode by "
                        "supplying arguments. For a list, use -h. ");
      sys.exit(1);
  else:
    # Run in console mode. 
    print_program_header();
    try:
      settings.getCmdLine();
    except SystemExit as e:
      sys.exit(str(e));
    except Exception as e:
      print >> sys.stderr, "An error occurred parsing your command line arguments. The error was: ";
      print >> sys.stderr, str(e);

      print >> sys.stderr, "";
      print >> sys.stderr, "The actual traceback was: ";
      traceback.print_tb(sys.exc_traceback);
      print >> sys.stderr, str(sys.exc_type);
      sys.exit(1);

    try:
      run_snipper(settings);
    except KeyboardInterrupt:
      print >> sys.stderr, "";
      print >> sys.stderr, "Program terminated by user.";
    except:
      if _SNIPPER_DEV:
        print >> sys.stderr, "(dev mode) -- exception thrown!";
        traceback.print_exc();
        pdb.post_mortem();
      elif _SNIPPER_DEBUG:
        print_debug();
      else:
        print >> sys.stderr, fill("An error has occurred. Please provide the following "
                                  "information to the developer, or run snipper with --debug to "
                                  "see additional debugging information. ");
        print >> sys.stderr, "";
        
        print >> sys.stderr, "Traceback:";
        traceback.print_tb(sys.exc_traceback);
        print >> sys.stderr, str(sys.exc_type);
        print >> sys.stderr, "";
        
        print >> sys.stderr, sys.exc_value;

if __name__ == "__main__":
  freeze_support();
  main();
