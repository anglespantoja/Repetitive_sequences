##### USER DEFYINED PARAMETERS ####
lenght_of_crRNA_with_PAM = 25
lenght_of_crRNA = lenght_of_crRNA_with_PAM - 3 #This has to be changed manually
PAMS = ["TTT"]
#Specify number of most repetitive sequences to be determined
REP = 100

####### DEFINYING FUNCTIONS ######

#Inverse complement
def InversComplement(input):
    output = ''
    for nucleotide in input:
        nucleotide = nucleotide.upper()

        if nucleotide == 'A':
            output += 'T'
        elif nucleotide == 'T':
            output += 'A'
        elif nucleotide == 'G':
            output += 'C'
        else:
            output += 'G'

    return(output[::-1]) # (::-1) - Will deliver reverse sequence

##### OPEN, READING, MODIFYING THE GENOME AND GENERATING COMPLEMENTARY STRAND #####

genome1 = open('SHORT_K12.fna', 'r')
next(genome1) #Jumps the first line which should be the FASTA title
genome1 = genome1.read()
genome1 = genome1.upper().replace("\n","")

genome1 = genome1 + genome1[:(lenght_of_crRNA_with_PAM-1)]
lenght = len(genome1)
print("Lenght of genome 1 is: " + str(lenght))
complementary_genome1 = InversComplement(genome1)

#### FOR SECOND GENOME ####

genome2 = open('SHORT_NISSLE.fna', 'r')
next(genome2) #Jumps the first line which should be the FASTA title
genome2 = genome2.read()
genome2 = genome2.upper().replace("\n","")

genome2 = genome2 + genome2[:(lenght_of_crRNA_with_PAM-1)]
lenght = len(genome2)
print("Lenght of genome 2 is: " + str(lenght))

complementary_genome2 = InversComplement(genome2)

####### FOR FIRST GENOME #########

import re
targets1 = []
for PAM in PAMS:
    new_targets1 = re.findall(rf'(?=({PAM}.{{22}}))', genome1)
    new_targets_complementary1 = re.findall(rf'(?=({PAM}.{{22}}))', complementary_genome1)
    targets1 = targets1 + new_targets1 + new_targets_complementary1

####### FOR SECOND GENOME ########## 

targets2 = []
for PAM in PAMS:
    new_targets2 = re.findall(rf'(?=({PAM}.{{22}}))', genome2)
    new_targets_complementary2 = re.findall(rf'(?=({PAM}.{{22}}))', complementary_genome2)
    targets2 = targets2 + new_targets2 + new_targets_complementary2

###### FROM LISTS, FIND DIFFERENT #####

matched_targets_in_genome1 = []
for target2 in targets2:
    print(str((targets2.index(target2)/len(targets2))*100))
    not_matched = "NO"
    memory = "NO"
    for target1 in targets1:
        if re.search(target2,target1):
            memory = "YES"
        else:
            not_matched = "YES"
    if not_matched == "YES" and memory == "NO":
        matched_targets_in_genome1.append(target2)

##### DETERMINATION OF MOST REPETITIVE crRNA TARGETED SEQUENCES

from collections import Counter
most_common_unique_targets_genome1 = [sequence for sequence in Counter(matched_targets_in_genome1).most_common(REP)]
print(most_common_unique_targets_genome1)

