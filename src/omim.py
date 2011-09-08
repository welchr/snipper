#!/usr/bin/env python
import sys
import re
import string
import traceback
from urllib2 import urlopen
from BeautifulSoup import BeautifulSoup, NavigableString, Tag
from collections import namedtuple

OMIM_SEARCH_URL = "http://omim.org/search?index=entry&start=1&limit=10&search={0}&sort=score+desc%2C+prefix_sort+desc&limit=10&field=title&prefix=*&prefix=%2B&date_created=-&date_updated=-";

OMIM_ENTRY_URL = "http://omim.org/entry/{0}";

def omim_find_gene(gene):
  url = OMIM_SEARCH_URL.format(gene);
  try:
    soup = BeautifulSoup(urlopen(url));
  except:
    raise;

  results_table = soup.find(name="table",attrs={'class':"entry-results-table"});
  if results_table == None:
    return None;

  for row in results_table.findAll("a"):
    href = row.get("href");
    if href != None:
      p = re.compile("\/entry\/(\d+)").search(href);
      if p != None:
        return p.groups()[0];

  return None;

def omim_get_entry(omim_id):
  return BeautifulSoup(urlopen(OMIM_ENTRY_URL.format(omim_id)));

def omim_parse_section(soup,section):
  soup = soup.find("a",{'name':section.code}).next.next.next.next;
  lines = [];
  
  for item in soup.contents:
    if isinstance(item,Tag):
      if item.name == 'br':
        lines += ["\n"];
        continue;
    
    if not isinstance(item,NavigableString):
      lines += [item.text];
    else:
      lines += [str(item)];

  section.text = "".join(lines);
  return section;

def omim_parse(soup):
  ident = omim_get_id(soup);
  title = omim_get_title(soup);

  article = OMIMArticle(ident);
  article.title = title;
 
  sections = omim_find_bio_sections(soup);
  for sect in sections:
    omim_parse_section(soup,sect);
    article.add_section(sect);
  
  for ref in omim_parse_references(soup):
    article.add_reference(ref);

  for pheno in omim_parse_phenos(soup):
    article.add_phenotype(pheno);

  return article;

def omim_find_bio_sections(soup):
  sections = [];
  for tag in soup.findAll("td",{'class':'section-title text-font'}):
    if tag.text == tag.text.upper():
      continue;
    
    name = tag.text;
    code = tag.find('a')['name'];

    sections.append( OMIMSection(code,name) );

  return sections;

def omim_get_title(soup):
  tag = soup.find("title");
  if tag != None:
    return tag.text;
  else:
    return None;

def omim_get_id(soup):
  return soup.find("span",{'class':'prefix'}).next.next;

def omim_parse_references(soup):
  refs = [];
  for tag in soup.findAll("td",{'class':"reference-text text-font"}):
    ref = OMIMReference();
    
    text = "";
    for i in tag.contents[0:6]:
      if isinstance(i,NavigableString):
        if i.strip() == "":
          continue;
        else:
          text += i + " ";
      else:
        text += i.text + " ";

    ref.text = text.rstrip();

    for i in tag.findAll("a"):
      if 'pubmed' in i['href'] and 'related' not in i.text:
        ref.pubmed_url = i['href'];
        ref.pmid = i.text;

    refs.append(ref);

  return refs;

def omim_parse_phenos(soup):
  def pheno_filter(t):
    return t.name == "td" \
      and t.has_key('class') \
      and 'gene-map-value' in t['class'] \
      and 'bottom' not in t['class'] \
      and not t.has_key('rowspan') \
      and 'value-border' in t['class']

  phenos = [];

  tags = soup.findAll(pheno_filter);
  for tag in tags:
    omim_id = tag.next.next.next.text;
    name = tag.text;

    pheno = OMIMPhenotype(omim_id);
    pheno.name = name;
    pheno.omim_url = OMIM_ENTRY_URL.format(omim_id);

    phenos.append(pheno); 

  return phenos;

def omim_get_gene_article(gene_symbol,verbose=False):
  try:
    omim_id = omim_find_gene(gene_symbol);
  except:
    print >> sys.stderr, "Warning: couldn't find OMIM information for gene %s.." % gene_symbol;
    if verbose:
      traceback.print_exc();  
    return None;
  
  try:
    soup = omim_get_entry(omim_id);
  except:
    print >> sys.stderr, "Error: tried to download OMIM article %s, but couldn't!" % omim_id;
    if verbose:
      traceback.print_exc();  
    return None;
    
  try:
    article = omim_parse(soup);
  except:
    print >> sys.stderr, "Error parsing OMIM article (%s) for gene %s, will be skipped.." % (omim_id,gene_symbol);
    if verbose:
      traceback.print_exc();  
    return None;
    
  return article;

class OMIMArticle:
  def __init__(self,id):
    self.id = id;
    self.sections = [];
    self.title = str();
    self.phenotypes = [];  
    self.references = [];

  def add_section(self,section):
    self.sections.append(section);

  def add_reference(self,ref):
    self.references.append(ref);

  def add_phenotype(self,pheno):
    self.phenotypes.append(pheno);

  def __str__(self):
    s = "";
    s += self.title;
    s += "\n\n";

    if len(self.phenotypes) > 0:
      s += "PHENOTYPES:\n\n";
      for pheno in self.phenotypes:
        s += str(pheno) + "\n";
      s += "\n";

    for sec in self.sections:
      s += sec.name.upper() + "\n\n";
      s += str(sec);
      s += "\n\n";

    s += "REFERENCES:\n\n";
    for ref in self.references:
      s += str(ref);
      s +=  "\n\n";

    return s;

class OMIMSection:
  def __init__(self,code,name):
    self.name = name;
    self.code = code;
    self.text = str();

  def __str__(self):
    return self.text;

  def __repr__(self):
    return "<OMIMSection [code:%s,name:%s] @ %s>" % (self.code,self.name,id(self));

class OMIMReference:
  def __init__(self):
    self.text = str();
    self.pmid = str();
    self.pubmed_url = str();

  def __str__(self):
    return self.text;

  def __repr__(self):
    return "<OMIMReference [pmid:%s,textbit:%s] @ %s" % (self.pmid,self.text[0:15],id(self));

class OMIMPhenotype:
  def __init__(self,id):
    self.name = str();
    self.omim_id = id;
    self.omim_url = str();

  def __str__(self):
    return "%s - OMIM: %s - %s" % (self.name,self.omim_id,self.omim_url);

def easy_test():
  from itertools import repeat
  import traceback
  import time
  import random

  genes = """
    RB1
  """.split();

  for gene in genes:
    omim_id = omim_find_gene(gene);
    if omim_id == None:
      print >> sys.stderr, "Warning: couldn't find OMIM ID for gene %s.." % (gene);
      continue;    

    soup = omim_get_entry(omim_id);
    article = omim_parse(soup);
    print article;    
  
    print "".join(repeat('=',80))
    print "";
  
    time.sleep(random.randrange(1,6));

def stress_test():
  from itertools import repeat
  import traceback
  import time
  import random

  genes = """
    RB1
    TCF7L2
    PDE8B
    KCNJ11
    RPA1
    DPH1
    TSR1
    SRR
    ARFIP1
    HIC1
    MNT
    TCF4
    LIPC
    PAX6
    CAV1
    ALMS1
    RBMS1
    SHOC2
    CKAP5
    DDB2
    ACP2
    MADD
    FAT3
    PHLDA2
    HHEX
    ACADS
    AKTIP
    FTO
    RBL2
    KCNQ1
    DUSP14
    SLC6A8
    SLC30A8
    HAUS7
    TREX2
    UPF2
    CDC123
    CAMK1D
    NUDT5
    UPF2
    VTI1A
    TRPM5
    PHLDA2
  """.split();

  for gene in genes:
    article = omim_get_gene_article(gene);
    print article;    
  
    print "".join(repeat('=',80))
    print "";
  
    time.sleep(random.randrange(1,6));

if __name__ == "__main__":
  try:
    stress_test();
  except:
    import pdb
    import traceback
    traceback.print_exc();
    pdb.post_mortem();
