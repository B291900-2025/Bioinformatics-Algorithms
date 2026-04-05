
import smith_waterman_p as SW
import  process_fasta_file as pff
import math
import random
import programme_settings as ps
import textwrap

ps.read()

randomloops =int(ps.settings["BUILD_EXPECT"]["random_loops"])
alphabet =ps.settings["BUILD_EXPECT"]["alphabet"]

print("*********************************************************")
print("    Running Random Trial for SmithWaterman               ")
print("    Random loops: %s                                " % randomloops)
print("    Building Expect Scores and Updating Main Settings     ")
print("**********************************************************")

counts_sizes = list()
count_bases = {}
count_bases_valid = {}


def count_the_bases():
    #Just for fun... count the bases
    global count_bases


    for i in alphabet:
        print(i + " "+str(ord(i)-ord('A')))
        #min is 'A' max is 'Z' (index 0 to 24)
        #other way chr(65) which is A
        count_bases[i] = 0



#build a function to call
def count_residues(seq, csizes = counts_sizes, cbases =count_bases):

    csizes.append(len(seq))


    for v in seq.upper():
        cbases[v] = cbases[v]+1

#        if v =='Z':
#            v =0

    return 0

def run_trial():
    global counts_sizes
    global count_bases
    global count_bases_valid

    count_the_bases()

    #feed into the fasta file processing
    ffile =ps.settings["DEFAULT"]["database"]
    pff.process_fasta_file(ffile, count_residues, 1)


    print("Count sizes sze: ",len(counts_sizes))
    print("Raw Base Counts from the Database: ",count_bases)

    #add up all counts
    sum =0
    for b in count_bases.values():
        sum += b
    print("Sum Bases: ",sum)

    valid_residues = str(ps.settings["BLAST"]["valid_residues"])

    #delete elements from the dictionary that are not valid bases
    all_residues =""
    for z in count_bases.keys():#nope
        if valid_residues.find(z)==-1:
            continue
        count_bases_valid[z] = count_bases[z]
        # express in per k of bases
        count_bases_valid[z] = math.floor(count_bases[z]/1000)

        all_residues += z * count_bases_valid[z]

    print("Count Bases : ",textwrap.fill(str(count_bases_valid),80))
    print("All residues in bases:(first 200) ",all_residues[1:100] )
    print("All residues in bases:(last 200) ",all_residues[len(all_residues)-100:] )


    #random shuffle the container all_residues and counts_sizes
    #select elements randomised and then random numbers...
    all_residues =''.join(random.sample(all_residues,len(all_residues)))
    random.shuffle(counts_sizes)

    print("Count Bases : ",textwrap.fill(str(count_bases_valid)))
    print("Random All residues in bases:(first 200) ",all_residues[1:100] )
    print("Random All residues in bases:(last 200) ",all_residues[len(all_residues)-100:] )
    print("Counts sizes [1:100] ",counts_sizes[1:100])

    #we don't want to overwrite the stored settings unless we allow this
    if ps.settings["BUILD_EXPECT"]["update_settings"]=="False":
        print("Update settings set to False, no settings saved!")
        return
    else:
        print("Updating settings...")

    #OK open a file to write SW outputs from randomised data
    with open(ps.settings["BUILD_EXPECT"]["random_file"], "a") as output_file:

        for counter in range(0,randomloops):
            seq1 =[]
            seq2 =[]
            s1 = counts_sizes[random.randint(0, len(counts_sizes)-1)]
            s2 = counts_sizes[random.randint(0, len(counts_sizes)-1)]

            for p in range(s1):
                seq1.append(all_residues[random.randint(0, len(all_residues)-1)])
            for p in range(s2):
                seq2.append(all_residues[random.randint(0, len(all_residues)-1)])

            seq1 = ''.join(seq1)
            seq2 = ''.join(seq2)
            #print("Seq1: "+seq1)
            #print("Seq2: "+seq2)

            if counter%1000 ==0:
                print("#"+str(counter))

            output_file.write(str(SW.perform_smith_waterman(seq1,seq2))+"\n")

#Run the trail
run_trial()




