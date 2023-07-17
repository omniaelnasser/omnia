
####Reference sequences for H1N1#######
#specify the path of sequence_IDs file on PC
Reference_IDs_path = 'D:/NEW_LAB/PROGRAMMING/project/Reference_IDs.txt'
try:                                      ##Error handling
    with open(Reference_IDs_path,'w') as Reference_IDs:      # open the file in write mode
        Reference_IDs.write('WGD05271\nWAB60225\nWDR22980\nWCJ17069\nUUB84697\nWDR22978\nWDR22947\nUZY70773\nWEI47079\nWDR22981')
except FileNotFoundError:
    print("File not exist")
##Download sequences from online NCBI database
from Bio import Entrez          #import Entrez to access NCBI database
from Bio import SeqIO
from Bio.Seq import Seq         ##import Seq for dealing with sequences
from Bio.SeqRecord import SeqRecord    ##Seqrecord for each sequence and related informations
Entrez.email = 's.Gamal2132@nu.edu.eg'
db = 'protein'       #Access protein database in NCBI
sequences_list = []  ## create a list to get a final list for all seq records
Reference_IDs = open(Reference_IDs_path,'r')   #open the IDs file in read mode
# specify the path of Refernce_sequence file on PC
Reference_sequences_path = 'D:/NEW_LAB/PROGRAMMING/project/Reference_sequences.txt'
#looping for each  ID in the file to download its sequence
for y in Reference_IDs:
    with Entrez.efetch(db=db, id=y, rettype='fasta', retmode='text') as handle:
        seq_record = SeqIO.read(handle, 'fasta')
        print(seq_record)
        sequences_list.append(seq_record)
###Create a single variable for each SeqRecord in the list
rec1 = (sequences_list[0])
rec2 = (sequences_list[1])
rec3 = (sequences_list[2])
rec4 = (sequences_list[3])
rec5 = (sequences_list[4])
rec6 = (sequences_list[5])
rec7 = (sequences_list[6])
rec8 = (sequences_list[7])
rec9 = (sequences_list[8])
rec10 = (sequences_list[9])
##create a new list for all variables
rec_list =[rec1,rec2,rec3,rec4,rec5,rec6,rec7,rec8,rec9,rec10]
#write all sequence records in Reference_sequence
SeqIO.write(sequences=rec_list, handle=Reference_sequences_path, format='fasta')
###################################################################

######Case sequences for H1N1
#specify the path of case_IDs file on PC
case_IDs_path = 'D:/NEW_LAB/PROGRAMMING/project/case_IDs.txt'
case_IDs = open(case_IDs_path,'w')    # open the file in write mode
#write accession No of the sequences
case_IDs.write('AJM70759\nANH22085\nWHU33647\nWCB23209\nWEI46848\nWEY08902\nWEY08900\nULA27163\nUXL02131\nBDO47307')
case_IDs.close()     ##close the file

##Download sequences from online NCBI database
from Bio import Entrez      #import Entrez to access NCBI database
from Bio import SeqIO
from Bio.Seq import Seq      ##import Seq for dealing with sequences
from Bio.SeqRecord import SeqRecord     ##Seqrecord for each sequence and related informations
Entrez.email = 's.Gamal2132@nu.edu.eg'
db = 'protein'         #Access protein database in NCBI
case_list = []         ## create a list to get a final list for all seq records
case_IDs = open(case_IDs_path,'r')     #open the case_IDs file in read mode
# specify the path of case_sequence file on PC
case_sequences_path = 'D:/NEW_LAB/PROGRAMMING/project/case_sequences.txt'
#looping for each  ID in the file to download its sequence
for z in case_IDs:
    with Entrez.efetch(db=db, id=z, rettype='fasta', retmode='text') as handle:
        set_record = SeqIO.read(handle, 'fasta')
        print(set_record)
        case_list.append(set_record)
###Create a single variable for each SeqRecord in the list
rec1 = (case_list[0])
rec2 = (case_list[1])
rec3 = (case_list[2])
rec4 = (case_list[3])
rec5 = (case_list[4])
rec6 = (case_list[5])
rec7 = (case_list[6])
rec8 = (case_list[7])
rec9 = (case_list[8])
rec10 = (case_list[9])
##create a new list for all variables
rec_list =[rec1,rec2,rec3,rec4,rec5,rec6,rec7,rec8,rec9,rec10]
#write all sequence records in case_sequences file
SeqIO.write(sequences=rec_list, handle=case_sequences_path, format='fasta')
###############################################################################
# Align the above Reference sequences and save their alignment in a file

##import clustalomegacommandline
from Bio.Align.Applications import ClustalOmegaCommandline
###specify the path of clustalomega
ClustalOmega = r"C:\Users\SOFTZONE\Downloads\Compressed\clustal-omega-1.2.2-win64\clustalo.exe"

###specify the path of infile (Reference_sequences) and outfile (Alignment_references)
file_in = r"D:/NEW_LAB/PROGRAMMING/project/Reference_sequences.txt"
file_out = r"D:/NEW_LAB/PROGRAMMING/project/Alignment_references.txt"

##running alignment and choose clustal format for the output result and force writing in the output file
clustalomega_cline = ClustalOmegaCommandline(ClustalOmega, infile=file_in, outfile=file_out, outfmt="clustal", force=True)
print(clustalomega_cline())

###Extraction of consensus sequence
from Bio import AlignIO      #dealing with alignment
from Bio.Align import AlignInfo     #extract information from alignment
#create a consensus from the MSA of reference sequences
alignment = AlignIO.read(handle= file_out, format="clustal")   #read alignment file
summary_align = AlignInfo.SummaryInfo(alignment)
consensus = summary_align.dumb_consensus()    # Seq object contain consensus sequence
#############################################################################

# Build a phylogentic tree for these the refernce sequences.
from Bio import AlignIO
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor

#read the alignment file
align = AlignIO.read(handle= file_out, format="clustal")
print(align)

#calculate distance matrix
calculator = DistanceCalculator('identity')
distancematrix = calculator.get_distance(align)
print(distancematrix)
## Create a DistanceTreeConstructor object
constructor = DistanceTreeConstructor()
#building phylogenetic tree using neighbor joining method
tree = constructor.nj(distancematrix)
Phylo.draw_ascii(tree)
from matplotlib import pylab
A = tree.common_ancestor({"name":"WDR22947.1"}, {"name":"WGD05271.1"})
A.color = "red"
B = tree.common_ancestor({"name":"WDR22980.1"}, {"name":"WDR22978.1"})
B.color = "blue"
tree.rooted = True
Phylo.draw(tree)
###############################################################
##preparing the file of MSA for case sequences
##specify the file path on PC
consensus_path = 'D:/NEW_LAB/PROGRAMMING/project/consensus_case_sequences.txt'
from Bio.Seq import Seq         ## dealing with sequences
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord   ##Seqrecord for each sequence and related informations
##create a seqrecord from consensus
consensus_record = SeqRecord(consensus, id='Reference sequence', description='Consensus sequence NA Influenza virus A')
##create a list of consensus record and other 10 case records created before
MSA_records =[consensus_record,rec1,rec2,rec3,rec4,rec5,rec6,rec7,rec8,rec9,rec10]
## writing 11 seqrecords in consensus_case_sequences
SeqIO.write(sequences=MSA_records, handle=consensus_path, format='fasta')
############################################################################


####Multiple sequence alignmet for consensus and case sequences
###specify the path of infile (consensus_case_sequences) and outfile (MSA_consensus_case_sequences)
file_sequences = r"D:/NEW_LAB/PROGRAMMING/project/consensus_case_sequences.txt"
file_MSA = r"D:/NEW_LAB/PROGRAMMING/project/MSA_consensus_case_sequences.txt"

##running alignment and choose clustal format for the output result and force writing in the output file
clustalomega_cline = ClustalOmegaCommandline(ClustalOmega, infile=file_sequences, outfile=file_MSA, outfmt="clustal", force=True)
print(clustalomega_cline())
#################################################################

##Phylogenetic tree
#read the alignment file
align_MSA = AlignIO.read(handle= file_MSA, format="clustal")
print(align_MSA)

#calculate distance matrix
calculator = DistanceCalculator('identity')
distancematrix = calculator.get_distance(align_MSA)
print(distancematrix)
## Create a DistanceTreeConstructor object
constructor = DistanceTreeConstructor()
#building phylogenetic tree using neighbor joining method
tree = constructor.nj(distancematrix)
Phylo.draw_ascii(tree)
from matplotlib import pylab
tree.rooted = True
Phylo.draw(tree)
########################################################################################################################
# The Percentage of Chemical Constituent from a separate Function File
# By Ahmed Mohamed Aziza
import Chemical_Constituents as cc                      # Import the function file in the code
file_input = './consensus_case_sequences.txt'           # Determine the relative Path of the desired file
cc.chemical_constituents(file_input)                    # Call the Function and the Resulted output is generated
########################################################################################################################
# The Dissimilar Regions in Sequences
# by Ahmed Mohamed Aziza
import dissimilar_regions as dsr                        # Import the function file
file_input = './MSA_consensus_case_sequences.txt'       # Determine the relative Path of the desired file
dsr.dissimilar_region(file_input)                       # Call the Function and the Resulted output is generateds
