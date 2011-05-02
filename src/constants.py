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

from util import convertFlank

# Default distance to use, in megabases. 
DEFAULT_DISTANCE = convertFlank('250kb');

# Default number of genes to retrieve within a given distance (if any exist, that is.) 
DEFAULT_NUMGENES = 9999;

# Default number of pubmed articles to retrieve. 
DEFAULT_PAPERNUM = 10;

SCANDB_WEB_URL = "http://www.scandb.org/";
SCANDB_PAPER_URL = "http://www.ncbi.nlm.nih.gov/pubmed/19933162";

# NCBI taxonomy ID for humans. 
TAXON_ID_HUMAN = '9606';

# NCBI URLs. 
OMIM_WEB_URL = "http://www.ncbi.nlm.nih.gov/omim/";
PUBMED_WEB_URL = "http://www.ncbi.nlm.nih.gov/pubmed/";

# Time to wait before a connection to MiMI, ScanDB, PubMed, etc. will timeout. 
CON_TIMEOUT = 10; # seconds

# SQLite table names. 
SQLITE_SNP_POS = "snp_pos";
SQLITE_REFFLAT = "refFlat";
SQLITE_REFSNP_TRANS = "refsnp_trans";
