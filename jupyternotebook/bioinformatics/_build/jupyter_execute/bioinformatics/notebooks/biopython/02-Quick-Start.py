#!/usr/bin/env python
# coding: utf-8

# # Chapter 2: Quick Start

# In[1]:


import IPython
print( IPython.sys_info() )


# In[2]:


import Bio
print(Bio.__version__)


# Following the tutorial here:
# 
# [Biopython Tutorial](http://biopython.org/DIST/docs/tutorial/Tutorial.pdf)

# # 2.2 Working with sequences

# In[4]:


from Bio.Seq import Seq
my_seq = Seq("AGTACACTGGT")
my_seq


# In[11]:


my_seq.complement()


# In[12]:


my_seq.reverse_complement()


# ## 2.4.1 FASTA and Bio.SeqIO

# ["NCBI Search"](https://www.ncbi.nlm.nih.gov/nuccore/?term=Cypripedioideae)

# ## Download FASTA files

# In[42]:


import wget
import tempfile

path = tempfile.mkdtemp()
url = 'https://raw.githubusercontent.com/biopython/biopython/master/Doc/examples/ls_orchid.fasta'
wget.download(url, out=path)
url = 'https://raw.githubusercontent.com/biopython/biopython/master/Doc/examples/ls_orchid.gbk'
wget.download(url, out=path)

import os

os.listdir(path)


# In[46]:


from pathlib import Path
fasta_path = Path(path + "/ls_orchid.fasta")
gbk_path = Path(path + "/ls_orchid.gbk")
fasta_path


# ## Review

# In[49]:


from Bio import SeqIO
for seq_record in SeqIO.parse( fasta_path, "fasta" ):
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))

