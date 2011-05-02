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

from ncbi import *
from textwrap import *
from util import *
from datetime import date,datetime
import sys

class PubmedArticle:
  def __init__(self):
    self.title = None;
    self.authors = [];
    self.so = None;
    self.pubdate = None;
    self.journal = None;
    self.pmid = None;
  
  def setTitle(self,title):
    self.title = periodCheck(title);

  def setSO(self,so):
    self.so = periodCheck(so);

  def setJournal(self,journal):
    self.journal = periodCheck(journal);

  def setPMID(self,id):
    self.pmid = id;

  def setPubDate(self,pubdate):
    self.pubdate = pubdate;

  def addAuthor(self,author):
    self.authors.append(author);

  def setAuthors(self,authors):
    if isIterable(authors):
      self.authors = authors;
    else:
      raise ValueError, "Author list must be iterable.";

  def write(self,out=sys.stdout,wrap_width=80):
    try:
      first_author = self.authors[0] + " et al";
      outstr = "-- %s. \"%s\" %s %s PMID: %s" % (first_author,self.title,self.journal,self.so,self.pmid);
      outstr = outstr.encode("utf-8");
      for line in wrap(outstr,initial_indent="",subsequent_indent="   ",width=wrap_width):
        out.write(line);
        out.write("\n");
    except UnicodeDecodeError:
      print >> sys.stderr, "Error printing article %s" % self.pmid;
      raise;
    except UnicodeEncodeError:
      print >> sys.stderr, "Error printing article %s" % self.pmid;
      raise;
  
  def __str__(self):
    try:
      first_author = self.authors[0] + " et al";
      outstr = "-- %s. \"%s\" %s %s PMID: %s" % (first_author,self.title,self.journal,self.so,self.pmid);
      outstr = outstr.encode("utf-8");
    except:
      print >> sys.stderr, "Error generating string for article %s." % str(self.pmid);
      raise;
    
    return outstr;

class Pubmed:
  _PMID = {};
  _BAD_PMIDS = set();
 
  @staticmethod
  def valueOf(pubmed_id):
    article = Pubmed._PMID.get(pubmed_id);

    if article != None:
      return article;
    elif pubmed_id in PubMed._BAD_PMIDS:
      article = PubmedArticle();
      article.setPMID(pubmed_id);
      return article;
    else:
      article = PubmedArticle();
      article.setPMID(pubmed_id);
      Pubmed._PMID[pubmed_id] = article;

    return article;

  @staticmethod
  def populate(pubmed_ids):
    if not isIterable(pubmed_ids):
      raise ValueError, "Must provide list of pubmed IDs to populate.";
      
    # Kill IDs that we have already seen before. 
    pubmed_ids = filter(lambda x: not Pubmed._PMID.has_key(x),pubmed_ids);

    # Kill IDs that we've seen before, but they were malformed for some reason. 
    pubmed_ids = filter(lambda x: x not in Pubmed._BAD_PMIDS,pubmed_ids);

    # Get the soup for all the pubmed ids. 
    if len(pubmed_ids) > 0:
      for chunk in subsets(pubmed_ids,50):
        soup = fetchPubmedSoup(chunk);
        for pmid in chunk:
          article = PubmedArticle();
          article.setPMID( pmid );
          root = findPMIDRoot(soup,pmid);

          if root:
            article.setTitle( extractTitle(root) );
            article.addAuthor( extractFirstAuthor(root) );
            article.setJournal( extractJournal(root) );
            article.setSO( extractSO(root) );
            article.setPubDate( extractDate(root) );

            Pubmed._PMID[pmid] = article;
          
          else:
            print >> sys.stderr, "Warning: could not find Pubmed root object for PMID %s.." % str(pmid);
            Pubmed._BAD_PMIDS.add(pmid);
  
  @staticmethod
  def writeAll(pubmed_ids,out=sys.stdout):
    # Terminal width. 
    term_width = terminal_size()[0];
    if out != sys.stdout: 
      term_width = 80;

    # Kill bad IDs. 
    pubmed_ids = filter(lambda x: x not in Pubmed._BAD_PMIDS,pubmed_ids);
    
    Pubmed.populate(pubmed_ids);
    articles = [Pubmed.valueOf(i) for i in pubmed_ids];
    for article in sorted(articles,key = lambda x: x.pubdate,reverse=True):
      Pubmed.valueOf(article.pmid).write(out,wrap_width=term_width);
      
  @staticmethod
  def writeReST(pubmed_ids,search_terms=[],out=sys.stdout):
    # Kill bad IDs. 
    pubmed_ids = filter(lambda x: x not in Pubmed._BAD_PMIDS,pubmed_ids);
    
    # Populate just in case. 
    Pubmed.populate(pubmed_ids);
    
    # Write articles in ReST format. 
    articles = [Pubmed.valueOf(i) for i in pubmed_ids];
    for article in sorted(articles,key = lambda x: x.pubdate,reverse=True):
      first_author = article.authors[0] + " et al";
      title_marked = markup_term(article.title,search_terms);
      fields = [unicode(i) for i in [first_author,title_marked,article.journal,article.so,article.pmid]];
      fields = map(lambda x: x.encode("utf-8","replace"),fields);
      outstr = "%s. \"%s\" %s %s `PMID %s`_" % tuple(fields);
      print >> out, "* %s" % outstr;
    
    print >> out, "";
    
    for id in pubmed_ids:
      print >> out, ".. _PMID %s: %s" % (str(id),PUBMED_WEB_URL + str(id));
      
    print >> out, "";
