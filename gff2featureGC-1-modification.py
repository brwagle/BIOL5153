#! /usr/bin/env python3.6

import sys

#function goes at the begining. a function to clear up a DNA sequence

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



def clean_seq(input_seq):
    clean = input_seq.upper()
    clean = clean.replace('N','')
    return clean
#clean variable is a private and is not accessible outside the function


def nuc_freq(sequence,base,sig_digs=2):
    #calculate the length of the sequence
    length=len(sequence)

    #count the number of this nucleotide

    count_of_base=sequence.count(base)

    # calculate the base frequency
    freq_of_base=count_of_base/length
    #return the frequencey and the lenght

    return (length,round(freq_of_base, sig_digs))


    #return()  #this will return the both variables

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
    start=int(fields[3])
    
    end=int(fields[4])
    #print(type, "\t", start,"\t", end)


    #extract feature from the genome

    fragment=genome[start-1:end]

    fragment=clean_seq(fragment)
    
    #print (clean)

    
    #sys.exti()

    if type=='CDS':
        cds+=fragment
        

    if type=='intron':
        intron+=fragment

    if type=='misc_feature':
        misc+=fragment

    if type=='repeat_region':
        repeats+=fragment
    if type=='rRNA':
        rrna+=fragment
    if type=='tRNA':
        trna+=fragment




            
#loop over the 4 nucleotide
        
types=[cds,intron,misc,repeats,rrna,trna]

#creating a dictionary to store keys and values

dict={

'exon':types[0],'intron':types[1],'misc_feature':types[2],'repeat_region':types[3],'rRNA':types[4],'tRNA':types[5]

 }



#Modification of assignement 7 to calculate nucleotide composition for each feature type

print("length and nucleotide composition of each feature type")
for feature_type, seq in dict.items():
    for nucleotide in ('A','C','G','T'):
        #calculate the nucleotide composition for each feature
        (feature_length,feature_comp)=nuc_freq(seq,base=nucleotide,sig_digs=2)
       
        print(feature_type.ljust(20)+ str(feature_length)+"\t"+nucleotide+ " "+str(feature_comp))





#Assignemnt 7 without modifications

gc_content_cds=((cds.count('G')+cds.count('C'))/len(cds))*100
gc_content_intron=((intron.count('G')+intron.count('C'))/len(intron))*100
gc_content_misc=((misc.count('G')+misc.count('C'))/len(misc))*100
gc_content_repeats=((repeats.count('G')+repeats.count('C'))/len(repeats))*100
gc_content_rrna=((rrna.count('G')+rrna.count('C'))/len(rrna))*100
gc_content_trna=((trna.count('G')+trna.count('C'))/len(trna))*100



# genome covered by each feature
length_cds=len(cds)
#print(length_cds)
genome_cds=(length_cds/len(genome))*100

length_intron=len(intron)
genome_intron=(length_intron/len(genome))*100

length_misc=len(misc)
genome_misc=(length_misc/len(genome))*100

length_repeats=len(repeats)
genome_repeats=(length_repeats/len(genome))*100

length_rrna=len(rrna)
genome_rrna=(length_rrna/len(genome))*100

length_trna=len(trna)
genome_trna=(length_trna/len(genome))*100
    

'''
#I tried these codes and didn't get good alighment, so I switch to different styles
print("exon \t %7d (%2.1f) \t %2.2f" % (length_cds,genome_cds,gc_content_cds))
print("intron \t %7d (%2.1f) \t %2.2f" % (length_intron,genome_intron,gc_content_intron))
print("misc_feature \t %7d (%2.1f) \t %2.2f" % (length_misc,genome_misc,gc_content_misc))
print("rrna\t %7d (%2.1f) \t %2.2f" % (length_rrna,genome_rrna,gc_content_rrna))
print("repeats_region\t %7d (%2.1f) \t %2.2f" % (length_repeats,genome_repeats,gc_content_repeats))
print("trna\t %7d (%2.1f) \t %2.2f" % (length_trna,genome_trna,gc_content_trna))
'''


#printing the results of assignment 7
for args in (('exon',length_cds,genome_cds,gc_content_cds), ('intron',length_intron,genome_intron,gc_content_intron),('misc_feature',length_misc,genome_misc,gc_content_misc),('rrna',length_rrna,genome_rrna,gc_content_rrna),('repeat_region',length_repeats,genome_repeats,gc_content_repeats),('trna',length_trna,genome_trna,gc_content_trna)):
    print ('{0:<20} {1:<6}({2:<4.1f}%) \t {3:<.2f}'.format(*args))


#close the GFF file
gff_in.close()
    
