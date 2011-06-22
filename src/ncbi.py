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

import __builtin__
import os
import sys
import urllib2
import re
from time import sleep,strptime
from datetime import date,datetime
from constants import *
from BeautifulSoup import BeautifulStoneSoup
from util import *

# NCBI query URLs.
EMAIL = "welchr@umich.edu";
EFETCH_URL = r"http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?tool=snipper&email=%s" % EMAIL;
ESEARCH_URL = r"http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?tool=snipper&email=%s" % EMAIL;
ELINK_URL = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?tool=snipper&email=%s' % EMAIL;
EPOST_URL = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/epost.fcgi?tool=snipper&email=%s' % EMAIL;
ESUMMARY_URL = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?tool=snipper&email=%s' % EMAIL;
OMIM_URL = r"http://www.ncbi.nlm.nih.gov/entrez/dispomim.cgi"; 

# Check to see if NCBI is online. 
# If this fails, there's either a problem with the internet connection, 
# or NCBI has banned this address (I don't have a way of knowing which, yet.) 
def check_ncbi_status():
  test_query = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=gene&id=23603&retmode=xml";
  try:
    response = urllib2.urlopen(test_query);
  except:
    msg = str(sys.exc_value);
    return (False,msg);
  
  return (True,'');

# Given a list of SNPs, find genomic positions. 
def findGenomicPositions(snps):
  snps = make_list(snps);
  snps = map(str,snps);
  url = EFETCH_URL + "&db=snp&retmode=xml&id=" + "&id=".join([i.replace("rs","") for i in snps]);
  soup = soupify(url);
  results = dict().fromkeys(snps);
  
  for root in soup.findAll("rs"):
    rsid_result = root.get('rsid');
    if rsid_result == None:
      continue;
    
    rsid = 'rs' + str(rsid_result);

    assembly = root.find('assembly',{'grouplabel':'reference'});
    pos = long(assembly.find('maploc')['physmapint'])+1;

    results[rsid] = pos;

  return results;

# Given a list of gene UIDs, find pubmed articles linked to them. 
# Returns dictionary where: 
# { gene UID } = [ list of pubmed ids ] 
def findPubmed(gene_uids):
  gene_uids = map(str,gene_uids);
  data = {};

  # Find pubmed links for gene. 
  final_url = ELINK_URL + "&dbfrom=gene&db=pubmed&sort=pub+date&id=" + "&id=".join(gene_uids);
  link_soup = soupify(final_url);
    
  for uid in gene_uids:
    # Find the root for this uid. 
    root = None;
    result = link_soup.find(text=uid);
    if result != None:
      root = result.parent.parent.parent;
    else:
      #print >> sys.stderr, "Warning: no pubmed links found for gene UID %s.." % uid;
      data[uid] = [];
      continue;

    # Find all pubmed IDs. 
    id_result = root.find(text='gene_pubmed');
    id_root = None;
    if id_result != None:
      id_root = id_result.parent.parent;
    else:
      #print >> sys.stderr, "Warning: no pubmed links found for gene UID %s.." % uid;
      data[uid] = [];
      continue;
    for tag in id_root.findAll('id'):
      data.setdefault(uid,[]).append(str(tag.string));

  return data;

# Easier method for requesting a URL. 
def soupify(url):
  try:
    if _SNIPPER_DEBUG:
      print "DEBUG: query executed from NCBI interface, URL was: %s" % url;
      print "";
  except:
    pass

  try:
    soup = BeautifulStoneSoup(urllib2.urlopen(url));
  except:
    msg = "";
    msg += "** Looks like a network error occurred when trying to contact NCBI.";
    msg += "** Check your network settings and make sure you can ping out.";
    msg += "** It is also entirely possible NCBI eutils are down, try again in a minute or two.";
    msg += "** URL attempted: " + url;
    msg += "** Error was: " + str(sys.exc_info()[1]);
    raise Exception, msg;

  # Rate limit queries. 
  sleep(1);
  
  return soup;

# Given a gene symbol and a search term, find pubmed article IDs that match both. 
def naiveSearchArticlesByTerm(gene_symbol,term,pnum):
  # Run the search, find the article IDs.
  term = gene_symbol + "+AND+" + term;
  search_url = ESEARCH_URL + "&db=pubmed&term=" + term;
  search_soup = soupify(search_url);

  # Pull out count. 
  count = int(search_soup.find('count').string);
  retmax = 20;

  # Iterate over result set. 
  article_ids = set();
  cont = True;
  for retstart in xrange(0,count,retmax):
    if cont:
      soup = soupify(search_url + "&retstart=%s&retmax=%s" % (retstart,retmax));
      for id in soup.findAll('id'):
        article_ids.add( id.string );
        if len(article_ids) >= pnum:
          cont = False;
          break;

  # Return results. 
  return article_ids;

# Given a gene symbol and list of search terms, find pubmed article IDs that match the gene
# and any of the search terms. 
def naiveSearchArticlesAnyTerm(gene_symbol,terms,pnum):
  # Run the search, find the article IDs.
  term = gene_symbol + "+AND+" + "(" + "+OR+".join(terms) + ")";
  search_url = ESEARCH_URL + "&db=pubmed&term=" + term;
  search_soup = soupify(search_url);

  # Pull out count. 
  count = int(search_soup.find('count').string);
  retmax = 20;

  # Iterate over result set. 
  article_ids = set();
  cont = True;
  for retstart in xrange(0,count,retmax):  
    if cont:
      soup = soupify(search_url + "&retstart=%s&retmax=%s" % (retstart,retmax));
      for id in soup.findAll('id'):
        article_ids.add( id.string );
        if len(article_ids) >= pnum:
          cont = False;
          break;

  # Return results. 
  return article_ids;

# Given a list of pubmed articles, search them for given terms. 
# Returns a dictionary, where: 
# { search term } = [ list of pubmed ids ] 
def searchArticlesByTerm(pmids,terms):
  if not isIterable(pmids): 
    raise ValueError, "Must provide list of pmids to search.";

  if not isIterable(terms):
    raise ValueError, "Must provide list of search terms.";
  
  pmids = [str(i) for i in make_list(set(pmids))];
  data = {};

  # First, post the pubmed IDs, and get the query/webenv keys.  
  for pmid_subset in subsets(pmids,100):
    post_url = EPOST_URL + "&db=pubmed&id=" + ",".join(pmid_subset);
    post_soup = soupify(post_url);
    post_querykey = post_soup.find('querykey').string;
    post_webenv = post_soup.find('webenv').string;

    # Next, for each search term, find the pmids that matched it
    # from the supplied list. 
    for term in terms:
      search_url = ESEARCH_URL + "&db=pubmed&term=%23" + "%s+AND+%s&WebEnv=%s&usehistory=y&sort=pub+date" %\
        (post_querykey,term,post_webenv);
      search_soup = soupify(search_url);
      for tag in search_soup.findAll('id'):
        data.setdefault(term,[]).append(str(tag.string));

  return data;

# Given a list of articles and search terms, finds the articles
# which match at least 1 of the search terms. 
# Returns: set of pmids
def searchArticlesAnyTerm(pmids,terms):
  if not isIterable(pmids): 
    raise ValueError, "Must provide list of pmids to search.";

  if not isIterable(terms):
    raise ValueError, "Must provide list of search terms.";
  
  pmids = [str(i) for i in make_list(set(pmids))];
  matching_pmids = set();

  # First, post the pubmed IDs, and get the query/webenv keys.  
  for pmid_subset in subsets(pmids,100):
    post_url = EPOST_URL + "&db=pubmed&id=" + ",".join(pmid_subset);
    post_soup = soupify(post_url);
    post_querykey = post_soup.find('querykey').string;
    post_webenv = post_soup.find('webenv').string;

    # Next, for each search term, find the pmids that matched it
    # from the supplied list. 
    term = "+OR+".join(terms);
    search_url = ESEARCH_URL + "&db=pubmed&term=%23" + "%s+AND+(%s)&WebEnv=%s&usehistory=y&sort=pub+date" %\
      (post_querykey,term,post_webenv);
    search_soup = soupify(search_url);
    for tag in search_soup.findAll('id'):
      matching_pmids.add(str(tag.string));

  return matching_pmids;

# Return esummary soup on pubmed articles. 
def fetchPubmedSoup(pmids): 
  pmids = map(str,pmids); # change each pmid into string
  final_url = ESUMMARY_URL + "&db=pubmed&id=" + ",".join(pmids) + "&retmax=" + str(len(pmids));
  soup = soupify(final_url);
  return soup

# From a pubmed soup, find a root for a given pmid. 
def findPMIDRoot(soup,pmid):
  pmid = str(pmid);
  tag = soup.find(lambda tag: tag.name == 'id' and tag.string == pmid);
  if tag != None:
    return tag.parent;
  else:
    return None;

# Extract pubdate from pubmed root. 
def extractPubDate(root):
  result = root.find(attrs={'name' : 'PubDate'});
  t = None

  if result != None:
    pubdate_str = result.string;
    if pubdate_str != None:
      pubdate_elements = pubdate_str.split();
      try:
        if len(pubdate_elements) == 2:
          t = datetime(*strptime(pubdate_str,"%Y %b")[0:6]);
        elif len(pubdate_elements) == 3:
          t = datetime(*strptime(pubdate_str,"%Y %b %d")[0:6]);
      except:
        pass
  return t;

# Extract EPubDate from pubmed root. 
def extractEPubDate(root):
  result = root.find(attrs={'name' : 'EPubDate'});
  t = None;

  if result != None:
    epubdate_str = result.string;
    if epubdate_str != None:
      epubdate_elements = epubdate_str.split();
      if len(epubdate_elements) == 2:
        t = datetime(*strptime(epubdate_str,"%Y %b")[0:6]);
      elif len(epubdate_elements) == 3:
        t = datetime(*strptime(epubdate_str,"%Y %b %d")[0:6]);

  return t;

# Extract "PubmedDate". 
def extractPubmedDate(root):
  result = root.find(attrs={'name' : 'pubmed', 'type' : 'Date'});
  t = None;

  if result != None:
    pubmed_date_str = result.string;
    if pubmed_date_str != None:
      pubmed_date_elems = pubmed_date_str.split();
      if len(pubmed_date_elems) == 2:
        t = datetime(*strptime(pubmed_date_elems[0],"%Y/%m/%d")[0:6]);

  return t;

# Extract date from pubmed root. 
# Returns a date object. 
def extractDate(root):
  pubdate = extractPubDate(root);
  epubdate = extractEPubDate(root);
  pubmed_date = extractPubmedDate(root);

  if pubdate != None:
    return pubdate;
  elif epubdate != None:
    return epubdate;
  elif pubmed_date != None:
    return pubmed_date;
  else:
    raise ValueError, "Aaaargh - NCBI didn't put a valid pubdate or epubdate!";

# From a pubmed root, extract the title. 
def extractTitle(root):
  if root:
    title = root.find(attrs={'name' : 'Title'});
    if title != None:
      return title.string;
    else:
      return None;
  else:
    return None

# From a pubmed root, extract the SO. 
def extractSO(root):
  SO = None;
  result = root.find(attrs={'name' : 'SO'});
  if result != None:
    SO = result.string;

  return SO;

# Extract first author. 
def extractFirstAuthor(root):
  author = "None";
  result = root.find(attrs={'name' : 'Author'});
  if result != None:
    author = result.string;

  return author;

# Extract journal. 
def extractJournal(root):
  journal = None;
  result = root.find(attrs={'name' : 'Source'});
  if result != None:
    journal = result.string;

  return journal;

# Extract the PMID. 
def extractPMID(root):
  id = None;
  result = root.find('id');
  if result != None:
    return result.string;

  return id;

# Queries Entrez and converts gene symbols into Entrez gene UIDs. 
#def symb2uid(gene_list):
#  uids = list();
#
#  # Fix: if gene_list is empty, and we submit the query.. it will attempt to return 
#  # every Homo sapien gene known to man. Ouch. 
#  if len(gene_list) == 0: 
#    return uids;
#
#  gene_term = "+OR+".join([i + "[SYMB]" for i in gene_list]) + "+AND+Homo+sapiens[orgn]";
#  
#  retmax = 50;
#  init_url = ESEARCH_URL + '&db=gene&term=' + gene_term + "&retmax=" + str(retmax);
#
#  # Find count. 
#  soup = soupify(init_url);
#  count = int(soup.find('count').string);
#
#  for retstart in xrange(0,count,retmax):
#    base_url = ESEARCH_URL + "&db=gene&term=%s&retstart=%s&retmax=%s" % (gene_term,retstart,retmax);
#    soup = soupify( base_url );
#    for id in soup.findAll('id'):
#      uids.append( str(id.string) );
#
#  return uids;

def symb2uid(gene_list):
  uids = list();

  # Fix: if gene_list is empty, and we submit the query.. it will attempt to return 
  # every Homo sapien gene known to man. Ouch. 
  if len(gene_list) == 0: 
    return uids;

  for chunk in subsets(gene_list,25): 
    gene_term = "+OR+".join([i + "[SYMB]" for i in chunk]) + "+AND+Homo+sapiens[orgn]";

    init_url = ESEARCH_URL + '&db=gene&retmax=200&term=' + gene_term;
    soup = soupify(init_url);
    for id in soup.findAll('id'):
      uids.append( str(id.string) );

  return uids;

# Queries Entrez Gene via efetch, and returns a soup of gene information. 
def fetchGeneSoup(uid_list):
  term = ','.join(uid_list);
  final_url = EFETCH_URL + '&db=gene&retmode=xml&id=' + term;

  soup = soupify(final_url);
  return soup;

# Determines if a gene record is discontinued or not. 
def isDiscontinued(gene_root):
  status = None;
  result = gene_root.find('gene-track_status');
  if result != None:
    status = str(result['value']);

  # Note to self: if gene-track_status doesn't exist, the record is probably bad, 
  # so we consider it to be "discontinued." 

  if status != 'live':
    return True;
  else:
    return False;

# Returns OMIM soup for a list of UIDs. 
def fetchOMIMSoup(uid_list):
  try:
    term = ','.join(uid_list);
  except TypeError:
    print "OMIM soup error: ";
    print uid_list;
    raise;
  
  final_url = EFETCH_URL + '&db=omim&retmode=xml&id=' + term;

  soup = soupify(final_url);
  return soup

# Given a OMIM soup and a uid, returns the soup corresponding to only that UID. 
def findOMIMRoot(soup,uid):
  result = soup.find(lambda tag: tag.name == 'mim-entry_mimnumber' and tag.string == uid);
  if result != None:
    return result.parent;
  else:
    print "Error: OMIM UID %s not found in return soup from NCBI.." % str(uid);
    return None;

# Extract OMIM text from OMIM root. 
def extractOMIMText(root):
  text = None;
  result = root.find('mim-text_text');
  if result != None:
    text = result.string;

  text = re.sub(r"\({(\d+)}\)",r"(\1)",text);
  text = re.sub(r"{(\d+):(.+?)}",r"\2",text);

  return text;

# Given a gene soup and a gene, find the root for that gene in the soup. 
# If the gene UID does not exist in the soup, this method will throw an exception. 
def findUIDRoot(soup,uid):
  # Check for UID in soup. 
  hook = soup.find(lambda tag: tag.name == "gene-track_geneid" and tag.string == uid);
  if hook == None: 
    raise Exception, "UID %s not found." % uid;
  else:
    uid_root = hook.parent.parent.parent;
    uid_root.extract();
    return uid_root;

# Given a root, find the UID. 
def extractUID(root):
  uid = None;
  result = root.find('gene-track_geneid');
  if result != None:
    uid = result.string;

  return str(uid);

# Given a gene root, find KEGG pathways. 
def extractKEGG(root):
  pathways = dict();
  results = root.findAll(text=re.compile('KEGG pathway:'));
  for i in results:
    # Quick and easy way to get rid of "KEGG pathway:"
    pathway_name = str(i);
    pathway_name = pathway_name[14:];
    
    # Find the URL. 
    pathway_url = str(i.parent.parent.find('other-source_url').string);

    # Store in tree. 
    pathways[pathway_name] = pathway_url;

  return pathways;

# Given a gene root, find phenotypes. 
def extractPhenotypes(root):
  phenos = set();
  proots = root.findAll(lambda tag: tag.name == "gene-commentary_type" and tag['value'] == "phenotype");
  for i in proots:
    phenos.add(i.next.next.next.string);

  return phenos;

# Given a gene root, find GO terms. 
# Returns a set. 
def extractGO(root):
  go_terms = set();
  go_hooks = root.findAll(tag='dbtag_db',text='GO');
  if go_hooks != None:
    for hook in go_hooks: 
      term = hook.parent.parent.parent.parent.find('other-source_anchor');
      if term != None:
        go_terms.add(str(term.string));
  
  return go_terms;

# Given a gene root, find gene symbol. 
def extractSymb(root):
  symb = None;
  result = root.find('gene-ref_locus');
  if result != None:
    symb = result.string;
  else:
    symb = None;

  return str(symb);

# Given a gene root, find synonyms for the gene symbol. 
def extractSyns(root):
  syns = set();
  result = root.findAll('gene-ref_syn_e');
  for r in result:
    s = str(r.string);
    if s != "None":
      syns.add(s);

  return syns;

# Find the proper gene symbol from a root and a list
# of possible symbols. 
# Requires: set of symbols, root
def findCorrectSymbol(symbs,root):
  syns = extractSyns(root);
  root_symbol = extractSymb(root);
  
  if root_symbol != None:
    syns.add(root_symbol);
  
  inter = set.intersection(syns,symbs);
  if len(inter) == 0: 
    print "Possible symbols: " + str(symbs);
    print "Possible synonyms: " + str(syns);
    print "Root symbol: " + str(root_symbol);
    raise Exception, "Couldn't find proper symbol for given root - contact developer.";
  else:
    return inter.pop();

# Given a gene root, find chromosome position. 
def extractLocus(root):
  locus = None;
  result = root.find('gene-ref_maploc');
  if result != None:
    locus = result.string;
  
  return str(locus);

# Given a gene root, find summary string. 
def extractSummary(root):
  summary = None;
  result = root.find('entrezgene_summary');
  if result != None:
    summary = result.string;
  
  return str(summary);

# Given a gene root, find gene type. 
def extractGeneType(root):
  type = None;
  result = root.find('entrezgene_type');
  if result != None:
    type = result.string;
  
  return str(type);

# Given a gene root, return the OMIM ID. 
def extractOMIM(root):
  omim = None;
  mim_root = root.find(text='MIM')
  if mim_root != None:
    id_result = mim_root.parent.parent.find('object-id_id');
    if id_result != None:
      omim = id_result.string;
  
  return str(omim);

# Given a gene root, return the gene description, if one exists. 
def extractDesc(root):
  desc = None;
  result = root.find('gene-ref_desc');
  if result != None:
    desc = result.string;

  return str(desc);

# Given a gene root, find the gene type. 
def extractType(root):
  type = None;
  result = root.find('entrezgene_type');
  if result != None:
    type = result['value'];

  return str(type);

# Given a gene root, extract GeneRIFs. 
# Returns a dictionary, where: 
# { pubmed id } = [ list of generif texts ]
def extractGeneRIF(root):
  data = dict();

  generif_roots = [p.parent for p in root.findAll(lambda tag: re.compile("gene-commentary_type").search(tag.name) and tag.string == "18")];
  for p in generif_roots:
    # Get the text. 
    text = p.find("gene-commentary_text");
    if text == None:
      continue;
    else:
      text = str(text.string);

    # Get the pubmed id. 
    pmid = p.find("pubmedid");
    if pmid == None:
      continue;
    else:
      pmid = str(pmid.string);

    # Store it. 
    data.setdefault(pmid,[]).append(text);
  
  return data;
