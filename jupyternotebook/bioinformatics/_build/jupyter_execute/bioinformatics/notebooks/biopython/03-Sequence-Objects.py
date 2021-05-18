#!/usr/bin/env python
# coding: utf-8

# # Chapter 3: Sequence objects

# * Sequences are essentially strings of letters like AGTACACTGGT, which seems very natural since this is the most common way that sequences are seen in biological le formats
# * The most important dierence between Seq objects and standard Python strings is they have dierent methods. Although the Seq object supports many of the same methods as a plain string, its translate() method diers by doing biological translation, and there are also additional biologically relevant methods like reverse_complement().

# ## 3.1 Sequences act like strings

# In[1]:


from Bio.Seq import Seq
my_seq = Seq("GATCG")
for index, letter in enumerate(my_seq):
    print("%i %s" % (index, letter))

print( "First element: " + my_seq[0] )
print( "Total elements: " + str( len( my_seq ) ) )


# ### Count elements

# In[2]:


from Bio.Seq import Seq
my_seq = Seq("GATCGATGGGCCTATATAGGATCGAAAATCGC")


# In[3]:


len(my_seq)


# In[4]:


my_seq.count("G")


# In[5]:


my_seq.count("GA")


# #### GC%

# In[6]:


100 * float(my_seq.count("G") + my_seq.count("C")) / len(my_seq)


# In[7]:


from Bio.SeqUtils import GC
GC( my_seq )


# ## 3.2 Slicing a sequence

# In[8]:


my_seq[4:12]


# In[9]:


my_seq[0::2]


# In[10]:


# Reverse
print( my_seq )
print( my_seq[::-1] )


# In[11]:


# Back to string:
str(my_seq)


# # 3.4 Concatenating or adding sequences

# In[12]:


from Bio.Seq import Seq
seq_1 = Seq("ATATAGG")
seq_2 = Seq("ACGT")
seq_1 + seq_2


# ## Loop-based concat

# In[13]:


from Bio.Seq import Seq
list_of_seqs = [Seq("ACGT"), Seq("AACC"), Seq("GGTT")]
concatenated = Seq("")
for s in list_of_seqs:
    concatenated += s
concatenated


# ## Join-based concat

# ## TODO: Not working but code shows this example?

# In[14]:


from Bio.Seq import Seq
concatenated = Seq('NNNNN').join([Seq("AAA")])
concatenated


# In[15]:


from Bio.Seq import Seq
contigs = [Seq("ATG"), Seq("ATCCCG"), Seq("TTGCA")]
spacer = Seq("N"*10)
spacer.join(contigs)


# ## 3.6 Nucleotide sequences and (reverse) complements

# * DNA: A=T, G≡C
# * RNA: A=U, G≡C
# 
# * Nucleic acid sequence of bases that can form a double- stranded structure by matching base pairs. For example, the **complementary sequence** to C-A-T-G (where each letter stands for one of the bases in DNA) is G-T-A-C
# * The **reverse complement** of a DNA sequence is formed by reversing the letters, interchanging A and T and interchanging C and G. Thus the reverse complement of ACCTGAG is CTCAGGT.

# In[16]:


from Bio.Seq import Seq
my_seq = Seq("GATCGATGGGCCTATATAGGATCGAAAATCGC")
print( my_seq )
print( my_seq.complement() )
print( my_seq.reverse_complement() )


# # 3.7 Transcription

# In[17]:


from Bio.Seq import Seq
coding_dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
template_dna = coding_dna.reverse_complement()
print("5' - " + coding_dna + " - 3'")
print("3' - " + template_dna[::-1] + " - 5'")


# ## Transcribe

# ### DNA Coding Strand transcription

# "As you can see, all this does is to replace T by U."

# In[18]:


coding_dna.transcribe()


# In[19]:


coding_dna


# ### Template Strand transcription

# In[20]:


template_dna.reverse_complement().transcribe()


# ### mRNA to DNA
# 
# "Again, this is a simple U -> T substitution:"

# In[21]:


from Bio.Seq import Seq
messenger_rna = Seq("AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG")
print( messenger_rna )
print( messenger_rna.back_transcribe() )


# # 3.8 Translation

# Translation is done from NCBI coding schemes.
# 
# * [NCBI Coding Scheme](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)

# ## Translate from mRNA

# In[22]:


from Bio.Seq import Seq
messenger_rna = Seq("AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG")
print( messenger_rna )
print( messenger_rna.translate() )


# ## Translate from coding DNA

# In[23]:


from Bio.Seq import Seq
coding_dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
print( coding_dna )
print( coding_dna.translate() )


# ### Take into account the DNA type

# In[24]:


print( coding_dna.translate(table="Vertebrate Mitochondrial") )


# ### NCBI table number (instead of using Vertebrate Mitochondrial)

# In[25]:


print( coding_dna.translate(table=2) )


# ### To Stop Codon

# In[26]:


print( coding_dna.translate(to_stop=True) )


# In[27]:


print( coding_dna.translate(to_stop=True, table="Vertebrate Mitochondrial") )


# ## Complete coding sequence CDS

# Telling BioPython to use a CDS means the coding becomes **Methionine** instead of **Valine**.

# In[28]:


from Bio.Seq import Seq
gene = Seq( "GTGAAAAAGATGCAATCTATCGTACTCGCACTTTCCCTGGTTCTGGTCGCTCCCATGGCA"
            "GCACAGGCTGCGGAAATTACGTTAGTCCCGTCAGTAAAATTACAGATAGGCGATCGTGAT"
            "AATCGTGGCTATTACTGGGATGGAGGTCACTGGCGCGACCACGGCTGGTGGAAACAACAT"
            "TATGAATGGCGAGGCAATCGCTGGCACCTACACGGACCGCCGCCACCGCCGCGCCACCAT"
            "AAGAAAGCTCCTCATGATCATCACGGCGGTCATGGTCCAGGCAAACATCACCGCTAA")


# In[29]:


gene.translate(table="Bacterial")


# In[30]:


gene.translate(table="Bacterial", cds = True)


# ### And it will throw an exception if wrong.
# 
# "Sequence length 296 is not a multiple of three"

# In[31]:


gene = Seq( "GTGAAAAAGATGCAATCTATCGTACTCGCACTTTCCCTGGTTCTGGTCGCTCCCATGGCA"
            "GCACAGGCTGCGGAAATTACGTTAGTCCCGTCACTAAAATTACAGATAGGCGATCGTGAT"
            "AATCGTGGCTATTACTGGGATGGAGGTCACTGGCGCGACCACGGCTGGTGGAAACAACAT"
            "TATGAATGGCGAGGCAGTCGCTGGCACCTACACGGACCGCCGCCACCGCCCGCCACCAT"
            "AAGAAAGCTCCTCATGATCATCACGGCGGTCATGGTCCAGGCAAACATCACCGCTAA")
gene.translate(table="Bacterial", cds = True)


# # 3.9 Translation Tables

# * Internally these use codon table objects derived from the NCBI information at [gc.prt](ftp://ftp.ncbi.nlm.nih.gov/entrez/misc/data/)
# * Also shown on https: [wprintgc.cgi](//www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi )

# In[32]:


from Bio.Data import CodonTable
standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
mito_table = CodonTable.unambiguous_dna_by_name["Vertebrate Mitochondrial"]


# In[33]:


print(standard_table)


# In[34]:


print(mito_table)


# ## Example: Find Codons

# In[35]:


standard_table.stop_codons


# In[36]:


standard_table.start_codons


# In[37]:


print( mito_table.forward_table["ACG"] )
print( mito_table.forward_table["CCG"] )


# # 3.10 Comparing Seq objects

# As of Biopython 1.65, sequence comparison only looks at the sequence and compares like the Python string.

# # 3.11 MutableSeq objects

# In[38]:


from Bio.Seq import Seq
my_seq = Seq("GCCATTGTAATGGGCCGCTGAAAGGGTGCCCGA")
my_seq[5] = "G" # Throws exception


# In[39]:


mutable_seq = my_seq.tomutable()
mutable_seq


# In[40]:


mutable_seq[0] = "C"
mutable_seq


# In[41]:


mutable_seq.remove("T")
mutable_seq


# ## Back to immutable

# In[42]:


immutable_again_seq = mutable_seq.toseq()
immutable_again_seq


# In[43]:


immutable_again_seq[5] = "G"


# # 3.12 UnknownSeq objects
# 
# A (more) memory efficient way of holding a sequence of unknown letters.

# In[44]:


from Bio.Seq import UnknownSeq
unk = UnknownSeq(20)
print(unk)


# * DNA is commonly labeled with N
# * Protein is commonly labeled with X

# In[45]:


unk_dna = UnknownSeq(20,character="N")
unk_dna


# In[46]:


unk_protein = UnknownSeq(20,character="X")
unk_protein


# In[ ]:




