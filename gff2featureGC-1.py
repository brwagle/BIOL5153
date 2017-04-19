#! /usr/bin/env python3.6

#Assignment 7

#Reading the watermelon.fsa file
my_file=open("watermelon.fsa", 'r')
my_dna=my_file.read()

#total lenght of genome in watermelon.fsa file
genome_length=len(my_dna)

#reading the watermelon.gff file
my_file2=open("watermelon.gff", 'r')

#Using grep and cut command in terminal start and stop position of each feature were redirected to a separate text file which are used for further processing


#exon length and GC content calculation
exon_locations=open("exon.txt")
coding_sequence=""
for line in exon_locations:
    positions=line.split('\t')
    start=int(positions[0])
    stop=int(positions[1])
    exon=my_dna[start:stop]
    coding_sequence=coding_sequence + exon


length_coding=len(coding_sequence)

#percentage of genome covered by exon
percent_of_genome_exon=(length_coding/genome_length)*100

#Calculation of GC content of exon
gc_content_exon=((coding_sequence.count('G') + coding_sequence.count("C"))/length_coding)*100


#intron length and GC content calculation

intron_locations=open("intron.txt")
noncoding_sequence=""
for line in intron_locations:
    positions=line.split('\t')
    start=int(positions[0])
    stop=int(positions[1])
    intron=my_dna[start:stop]
    noncoding_sequence=noncoding_sequence + intron

# Calculation of GC content of intron
gc_content_intron=((noncoding_sequence.count('G') + noncoding_sequence.count("C"))/len(noncoding_sequence))*100

length_noncoding=len(noncoding_sequence)

# percentage of genome covered by introns
percent_of_genome_intron=(length_noncoding/genome_length)*100


#misc-feature length and GC content calculation

misc_locations=open("misc_feature.txt")
misc_sequence=""
for line in misc_locations:
    positions=line.split('\t')
    start=int(positions[0])
    stop=int(positions[1])
    misc=my_dna[start:stop]
    misc_sequence=misc_sequence + misc

# Calculation of GC content in misc_feature sequence
gc_content_misc=((misc_sequence.count('G') + misc_sequence.count("C"))/len(misc_sequence))*100


length_misc=len(misc_sequence)

# percentage of genome covered by misc_feature
percent_of_genome_misc=(length_misc/genome_length)*100


#rRNA length and GC content calculation

rRNA_locations=open("rRNA.txt")
rRNA_sequence=""
for line in rRNA_locations:
    positions=line.split('\t')
    start=int(positions[0])
    stop=int(positions[1])
    rRNA=my_dna[start:stop]
    rRNA_sequence=rRNA_sequence + rRNA

#Calculation of GC content in rRNA
gc_content_rRNA=((rRNA_sequence.count('G') + rRNA_sequence.count("C"))/len(rRNA_sequence))*100

length_rRNA=len(rRNA_sequence)

# percentage of genome covered by rRNA
percent_of_genome_rRNA=(length_rRNA/genome_length)*100



#repeat region length and GC content calculation

repeat_locations=open("repeat.txt")
repeat_sequence=""
for line in repeat_locations:
    positions=line.split('\t')
    start=int(positions[0])
    stop=int(positions[1])
    repeat=my_dna[start:stop]
    repeat_sequence=repeat_sequence + repeat


#Calculation of GC content in repeat region
gc_content_repeat=((repeat_sequence.count('G') + repeat_sequence.count("C"))/len(repeat_sequence))*100


length_repeat=len(repeat_sequence)

# percentage of genome covered by repeat region
percent_of_genome_repeat=(length_repeat/genome_length)*100




#tRNA length and GC content calculation

tRNA_locations=open("tRNA.txt")
tRNA_sequence=""
for line in tRNA_locations:
    positions=line.split('\t')
    start=int(positions[0])
    stop=int(positions[1])
    tRNA=my_dna[start:stop]
    tRNA_sequence=tRNA_sequence + tRNA


# Calculation of GC content in tRNA
gc_content_tRNA=((tRNA_sequence.count('G') + tRNA_sequence.count("C"))/len(tRNA_sequence))*100


length_tRNA=len(tRNA_sequence)

# percentage of genome covered by rRNA
percent_of_genome_tRNA=(length_tRNA/genome_length)*100




print("exon ",length_coding,"(",round(percent_of_genome_exon, 2),"%)",round(gc_content_exon, 2))
print("intron ",length_noncoding,"(",round(percent_of_genome_intron, 2),"%)",round(gc_content_intron, 2))
print("misc_feature ",length_misc,"(",round(percent_of_genome_misc, 2),"%)",round(gc_content_misc, 2))
print("rRNA ",length_rRNA,"(",round(percent_of_genome_rRNA, 2),"%)",round(gc_content_rRNA, 2))
print("repeat_region ",length_repeat,"(",round(percent_of_genome_repeat, 2),"%)",round(gc_content_repeat, 2))
print("tRNA ",length_tRNA,"(",round(percent_of_genome_tRNA, 2),"%)",round(gc_content_tRNA, 2))
