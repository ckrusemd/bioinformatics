#!/usr/bin/env python
# coding: utf-8

# # Chapter 4: Sequence annotation objects

# * Chapter 3 introduced the sequence classes. Immediately above the Seq class is the Sequence Record or **SeqRecord class**, defined in the Bio.SeqRecord module. 
# * This class allows higher level features such as identiers and features (as SeqFeature objects) to be associated with the sequence, and is used throughout the sequence input/output interface Bio.SeqIO described fully in Chapter 5.

# ## 4.1 The SeqRecord object

# The SeqRecord (Sequence Record) class is dened in the Bio.SeqRecord module. This class allows higher level features such as identiers and features to be associated with a sequence (see Chapter 3), and is the basic data type for the Bio.SeqIO sequence input/output interface (see Chapter 5).
# 
# ** The SeqRecord class itself is quite simple, and oers the following information as attributes:**
# * **.seq** { The sequence itself, typically a Seq object.
# * **.id** { The primary ID used to identify the sequence { a string. In most cases this is something like an accession number.
# * **.name** { A \common" name/id for the sequence { a string. In some cases this will be the same as the accession number, but it could also be a clone name. I think of this as being analogous to the LOCUS id in a GenBank record.
# * **.description** { A human readable description or expressive name for the sequence { a string.
# * **.letter_annotations** { Holds per-letter-annotations using a (restricted) dictionary of additional information about the letters in the sequence. The keys are the name of the information, and the information is contained in the value as a Python sequence (i.e. a list, tuple or string) with the same length as the sequence itself. This is often used for quality scores (e.g. Section 20.1.6) or secondary structure information (e.g. from Stockholm/PFAM alignment les).
# * **.annotations** { A dictionary of additional information about the sequence. The keys are the name of
# the information, and the information is contained in the value. This allows the addition of more
# \unstructured" information to the sequence.
# * **.features** { A list of SeqFeature objects with more structured information about the features on a sequence
# (e.g. position of genes on a genome, or domains on a protein sequence). The structure of sequence
# features is described below in Section 4.3.
# * **.dbxrefs** - A list of database cross-references as strings.

# ## 4.2 Creating a SeqRecord
# 
# Most commonly, you would download a sequence and read it. It is possible to generate your own sequence with meta-data.

# In[1]:


from Bio.Seq import Seq
simple_seq = Seq("GATC")
simple_seq


# In[2]:


from Bio.SeqRecord import SeqRecord
simple_seq_r = SeqRecord(simple_seq)
simple_seq_r


# In[3]:


simple_seq_r.id = "AC12345"
simple_seq_r.description = "Made up sequence I wish I could write a paper about"
simple_seq_r.annotations["evidence"] = "There is evidence"
simple_seq_r.annotations["author"] = "CKS"
simple_seq_r


# Working with _per-letter-annotations_ is similar, letter_annotations is a dictionary like attribute which will let you assign any Python sequence (i.e. a string, list or tuple) which has the same length as the sequence:

# In[4]:


simple_seq_r.letter_annotations["phred_quality"] = [40, 40, 38, 30]
print(simple_seq_r.letter_annotations)


# ## 4.2.2 SeqRecord objects from FASTA files

# * This example uses a fairly large FASTA file containing the whole sequence for Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, originally downloaded from the NCBI. 
# * This file is included with the Biopython unit tests under the GenBank folder, or online [NC_005816.fna](https://raw.githubusercontent.com/biopython/biopython/master/Tests/GenBank/NC_005816.fna) from our website.

# In[5]:


import wget
import tempfile
path = tempfile.mkdtemp()
url = "https://raw.githubusercontent.com/biopython/biopython/master/Tests/GenBank/NC_005816.fna"
wget.download(url, out=path)

from pathlib import Path
fasta_path = Path(path + "/NC_005816.fna")


# In[6]:


from Bio import SeqIO
record = SeqIO.read(fasta_path, "fasta")
print( len(record.seq) )
print( record.name )
print( record.description )
print( record.id )
record.seq


# ## 4.2.3 SeqRecord objects from GenBank files

# In[7]:


import wget
import tempfile
path = tempfile.mkdtemp()
url = "https://raw.githubusercontent.com/biopython/biopython/master/Tests/GenBank/NC_005816.gb"
wget.download(url, out=path)

from pathlib import Path
genbank_path = Path(path + "/NC_005816.gb")


# In[8]:


from Bio import SeqIO
record = SeqIO.read(genbank_path, "genbank")
record


# In[9]:


print( len(record.seq) )
print( record.name )
print( record.description )
print( record.id )
print( record.dbxrefs )
record.annotations


# # 4.3 Feature, location and position objects

# ## 4.3.1 SeqFeature objects
# 
# * The key idea about each SeqFeature object is to describe a region on a parent sequence, typically a SeqRecord object. 
# * That region is described with a location object, typically a range between two positions (see Section 4.3.2 below).
# * The attributes of a SeqFeature are:
# 
# * **.type** { This is a textual description of the type of feature (for instance, this will be something like `CDS'
# or `gene').
# 
# * **.location** { The location of the SeqFeature on the sequence that you are dealing with, see Section 4.3.2 below. The SeqFeature delegates much of its functionality to the location object, and includes a number of shortcut attributes for properties of the location:
# * **- .ref** { shorthand for .location.ref { any (dierent) reference sequence the location is referring to. Usually just None.
# * **- .ref_db** { shorthand for .location.ref_db { species the database any identier in .ref refers to. Usually just None.
# * ** - .strand** { shorthand for .location.strand { the strand on the sequence that the feature is located on. For double stranded nucleotide sequence this may either be 1 for the top strand, 􀀀1 for the bottom strand, 0 if the strand is important but is unknown, or None if it doesn't matter. This is None for proteins, or single stranded sequences.
# 
# * **.qualifiers** { This is a Python dictionary of additional information about the feature. The key is some kind of terse one-word description of what the information contained in the value is about, and the value is the actual information. For example, a common key for a qualier might be \evidence" and the value  might be \computational (non-experimental)." This is just a way to let the person who is looking at the feature know that it has not be experimentally (i. e. in a wet lab) conrmed. Note that other the
# value will be a list of strings (even when there is only one string). This is a reection of the feature
# tables in GenBank/EMBL les.
# * **.sub_features** { This used to be used to represent features with complicated locations like `joins' in Gen-
# Bank/EMBL les. This has been deprecated with the introduction of the CompoundLocation object,
# and should now be ignored.
# 
# 

# ## 4.3.2 Positions and locations

# * **position** { This refers to a single position on a sequence, which may be fuzzy or not. For instance, 5, 20,
# <100 and >200 are all positions.
# * **location** { A location is region of sequence bounded by some positions. For instance 5..20 (i. e. 5 to 20) is
# a location.

# ### 4.3.2.1 FeatureLocation object
# Unless you work with eukaryotic genes, most SeqFeature locations are extremely simple - you just need
# start and end coordinates and a strand. That's essentially all the basic FeatureLocation object does.
# In practise of course, things can be more complicated. First of all we have to handle compound locations
# made up of several regions. Secondly, the positions themselves may be fuzzy (inexact).
# ### 4.3.2.2 CompoundLocation object
# Biopython 1.62 introduced the CompoundLocation as part of a restructuring of how complex locations made
# up of multiple regions are represented. The main usage is for handling `join' locations in EMBL/GenBank
# les.
# ### 4.3.2.3 Fuzzy Positions
# So far we've only used simple positions. One complication in dealing with feature locations comes in the
# positions themselves. In biology many times things aren't entirely certain (as much as us wet lab biologists
# try to make them certain!). For instance, you might do a dinucleotide priming experiment and discover that
# the start of mRNA transcript starts at one of two sites. This is very useful information, but the complication
# comes in how to represent this as a position. To help us deal with this, we have the concept of fuzzy positions.
# **Basically there are several types of fuzzy positions, so we have ve classes do deal with them:**
# * ExactPosition { As its name suggests, this class represents a position which is specied as exact along
# the sequence. This is represented as just a number, and you can get the position by looking at the
# position attribute of the object.
# * BeforePosition { This class represents a fuzzy position that occurs prior to some specied site. In Gen-
# Bank/EMBL notation, this is represented as something like <13
# , signifying that the real position is
# located somewhere less than 13. To get the specied upper boundary, look at the position attribute
# of the object.
# * AfterPosition { Contrary to BeforePosition, this class represents a position that occurs after some speci
# ed site. This is represented in GenBank as >13
# , and like BeforePosition, you get the boundary
# number by looking at the position attribute of the object.
# * WithinPosition { Occasionally used for GenBank/EMBL locations, this class models a position which
# occurs somewhere between two specied nucleotides. In GenBank/EMBL notation, this would be
# represented as `(1.5)', to represent that the position is somewhere within the range 1 to 5. To get the
# information in this class you have to look at two attributes. The position attribute species the lower
# boundary of the range we are looking at, so in our example case this would be one. The extension
# attribute species the range to the higher boundary, so in this case it would be 4. So object.position
# is the lower boundary and object.position + object.extension is the upper boundary.
# * OneOfPosition { Occasionally used for GenBank/EMBL locations, this class deals with a position where
# several possible values exist, for instance you could use this if the start codon was unclear and there
# where two candidates for the start of the gene. Alternatively, that might be handled explicitly as two
# related gene features.
# * UnknownPosition { This class deals with a position of unknown location. This is not used in Gen-
# Bank/EMBL, but corresponds to the `?' feature coordinate used in UniProt.

# ### 4.3.2.4 Location testing

# In[10]:


from Bio import SeqIO
my_snp = 4350

for feature in record.features:
    if my_snp in feature:
        print("%s %s" % (feature.type, feature.qualifiers.get("db_xref")))


# In[11]:


from Bio import SeqIO
my_snp = 1888
for feature in record.features:
    if my_snp in feature:
        print("%s %s" % (feature.type, feature.qualifiers.get("db_xref")))


# ## 4.3.3 Sequence described by a feature or location

# * **SeqFeature** is a description of how to get a location from a parent sequence

# In[12]:


from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
parent_sequence = Seq("ACCGAGACGGCAAAGGCTAGCATAGGTATGAGACTTCCTTCCTGCCAGTGCTGAGGAACTGGGAGCCTAC")
feature = SeqFeature(FeatureLocation(5, 18), type="gene", strand=-1) # <- notice that it does not reference the parent_sequence


# In[13]:


feature_seq = seq[feature.location.start:feature.location.end]
print(feature_seq)


# In[14]:


feature_seq = seq[feature.location.start:feature.location.end].reverse_complement()
print(feature_seq)


# Conversely, the **SeqFeature** also has the extract functionality where a sequence is referenced.

# In[15]:


feature_seq = feature.extract(seq)
print(feature_seq)


# # 4.4 Comparison (and differences)

# In[16]:


from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
record1 = SeqRecord(Seq("ACGT"), id="test")
record2 = SeqRecord(Seq("ACGT"), id="test")


# ```
# record1 == record2
# ```
# SeqRecord comparison is deliberately not implemented. Explicitly compare the attributes of interest.

# In[17]:


print( record1.id )
print( record2.id )


# In[18]:


print( record1.seq )
print( record2.seq )


# # 4.5 References

# ````
# Bio.SeqFeature.Reference class
# 
# * journal
# * title 
# * authors
# * medline_id 
# * pubmed_id 
# * comment 

# # 4.6 The format method

# The format method of **SeqRecord** allows you to print in a format acceptable by SeqIO.

# In[19]:


from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
record = SeqRecord(
    Seq(
        "MMYQQGCFAGGTVLRLAKDLAENNRGARVLVVCSEITAVTFRGPSETHLDSMVGQALFGD"
        "GAGAVIVGSDPDLSVERPLYELVWTGATLLPDSEGAIDGHLREVGLTFHLLKDVPGLISK"
        "NIEKSLKEAFTPLGISDWNSTFWIAHPGGPAILDQVEAKLGLKEEKMRATREVLSEYGNM"
        "SSAC"
        ),
        id="gi|14150838|gb|AAK54648.1|AF376133_1",
        description="chalcone synthase [Cucumis sativus]"
        )


# In[20]:


print(record.format("fasta"))


# # 4.7 Slicing a SeqRecord

# In[21]:


from Bio import SeqIO
record = SeqIO.read(genbank_path, "genbank")
record


# In[22]:


len(record)


# In[23]:


len(record.features)


# Selecting a specific gene that has already been sequenced.

# In[24]:


print(record.features[21])


# Reversely, select a parent record and get the associated features.

# In[25]:


sub_record = record[4300:4800]
print( len( sub_record ) )
sub_record


# In[26]:


sub_record.features


# In[27]:


print( sub_record.features[0] )


# In[28]:


print( sub_record.features[1] )


# In[29]:


sub_record.description


# In[30]:


sub_record.name


# # 4.8 Adding SeqRecord objects
# 
# You can add SeqRecord objects together, giving a new SeqRecord

# In[31]:


from Bio import SeqIO
record = SeqIO.read(genbank_path, "genbank")
record


# In[32]:


record


# In[33]:


record[2000:]


# In[34]:


record[:2000]


# In[35]:


len(record.features)


# ## Shifting

# In[36]:


shifted = record[2000:] + record[:2000]
shifted


# In[37]:


len(shifted.features)


# # 4.9 Reverse-complementing SeqRecord objects

# In[ ]:




