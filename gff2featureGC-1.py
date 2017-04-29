#! /usr/bin/env python3.6


import sys
import collections



#function goes at the begining. a function to clear up a DNA sequence

def clean_seq(input_seq):
    clean = input_seq.upper()
    clean = clean.replace('N','')
    return clean
#clean variable is a private and is not accessible outside the function

#Function to calculate lenght and composition of sequence
def nuc_freq(sequence,base1,base2,sig_digs=2):
    #calculate the length of the sequence
    length=len(sequence)

    #genome covered with the feature
    genome_cover=(len(sequence)/len(genome))*100

    #count the number of this nucleotide

    count_of_base1=sequence.count(base1)

    count_of_base2=sequence.count(base2)

    # calculate the gc content
    gc_content=((count_of_base1+count_of_base2)/length)*100

    #return the frequencey and the lenght
    return (length,genome_cover,round(gc_content, sig_digs))
    #return()  #this will return the both variables

#Function to generate complement sequence
def reverse_compl(dna_seq):
    replacement1=dna_seq.replace('A','t')
    replacement2=replacement1.replace('T','a')
    replacement3=replacement2.replace('C','g')
    replacement4=replacement3.replace('G','c')
    return (replacement4.upper())





usage = sys.argv[0] + ": watermelon.fsa watermelon.gff"

if len(sys.argv) < 3:
    print(usage)
    sys.exit("\nThis script requires both watermelon FSA file and a watermelon GFF file\n")

watermelon_gff = sys.argv[1]
watermelon_fsa= sys.argv[2]

#print(gff + "\n" + genome)





gff_file="watermelon.gff"
fsa_file="watermelon.fsa"


#open the files 
gff_in=open(gff_file,'r')
fsa_in=open(fsa_file,'r')



#Creating dictionary
    
#key=feature type, value=concatenation of all sequences of that type-not useful for anything other than calculating AT/GC Content
feature_sequences={}

#key exon, value=sequence
exon_sequences={}

#key is gene and value is concatenation of all exon sequences in a gene
gene_sequences={}


#Declare a variable
genome=''

line_number=0



for line in fsa_in:
    #print(str(line_number) + ":" + line)

    line=line.rstrip('\n')
    
    if line_number > 0:
        genome+=line
        
    line_number+=1


#did we get the genome correctly
#print(len(genome))

#close the file fsa
fsa_in.close()


cds     = ''
trna    = ''
rrna    = ''
intron  = ''
misc    = ''
repeats = ''



for line in gff_in:
    line=line.rstrip('\n')
    #types=line.split('type ')
   # other_type=types[len(types)-1]
   # print(other_type)

    fields=line.split('\t')
    type=fields[2]
    strand=fields[6]

        #print(gene[0])
    start=int(fields[3])
    
    end=int(fields[4])
    #print(type, "\t", start,"\t", end)


    #extract feature from the genome

    fragment=genome[start-1:end]

    fragment=clean_seq(fragment)


        
    if type in feature_sequences:
        feature_sequences[type] +=fragment
    else:
        feature_sequences[type]=fragment

    if type=='CDS':
        gene_feature=fields[8]
        gene1=gene_feature.split(';')
        gene=gene1[0]
        gene_sequence=genome[start-1:end]

        if strand=='-':
            #print("Before ")
            complement_sequence=reverse_compl(gene_sequence)
            exon_sequences[gene]=complement_sequence
            #print("After ")
        else:
            exon_sequences[gene]=gene_sequence
        
    #print (clean)
    #print(gene[0])
    
#close the GFF file
gff_in.close()
    

        

# order the exon sequences
ordered_exons_sequences = collections.OrderedDict(sorted(exon_sequences.items()))
    
for exon, seq in ordered_exons_sequences.items():
    print(">", exon,"\n",seq)  #print out the exon and the sequences


    
for feature, sequence in feature_sequences.items():

    #calculate the nucleotide composition for each feature
    (feature_length,cover, feature_comp)=nuc_freq(sequence,base1='C',base2='G',sig_digs=2)
     
    print(feature.ljust(20), str(feature_length), "(%1.1f" %cover, "%)", "\t",str(feature_comp)+"%")


    




#function for reverse complement
#function for G+C
#print/store/build cds for each gene
#capture the gene name and _ + if _ call reverse function
