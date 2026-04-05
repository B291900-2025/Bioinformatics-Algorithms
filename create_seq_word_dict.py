############################################################################
#   Use to create the data structure from words of the query seq in BLAST  #
############################################################################
import blosum as bl
import math
from collections import defaultdict
import textwrap

import programme_settings

programme_settings.read()

dist = bl.BLOSUM(int(programme_settings.settings["DEFAULT"]["blosum"]))

tscore = int(programme_settings.settings["BLAST"]["tscore"])
tscore_sim = int(programme_settings.settings["BLAST"]["tscore_sim"])

#list of valid residues
residues =programme_settings.settings["BLAST"]["valid_residues"]

def create_word_dict(seq):
    global tscore


    word_dictionary =defaultdict(list)
    word_size =int(programme_settings.settings["BLAST"]["word_size"])
    seq=seq.upper()

    len_seq=len(seq)/word_size
    len_seq =math.floor(len_seq) * word_size



    for i in range(0, len_seq):
        #print(i)

        # indexes
        current_word_pos = i

        s = seq[i:i+word_size]
        if len(s) % word_size !=0:
            continue

        #score the words and see if they are suitable
        score=0
        for j in range(0,word_size):
            score = score + dist[s[j]][s[j]]
        if score < tscore:
            continue

        word_dictionary[s].append(current_word_pos)

        expand_to_similar(word_dictionary,seq,current_word_pos,word_size)

        #Expand words to mis=1


    return word_dictionary

def expand_to_similar(word_dictionary, seq, current_word_pos, word_size):

    global dist
    global tscore
    global tscore_sim
    global residues

#Expand words to mis=1
    for i in range(0, word_size):
        s = seq[current_word_pos:current_word_pos + word_size]
        if len(s) % word_size != 0:
            continue
        #base score-these bases
        for k in range(0,len(residues)):
            s2 = list(s)
            s2[i]=residues[k]
            s2 = "".join(s2)

            #skip matching sequence
            if s2 ==s :
                continue


            score = 0
            for j in range(0, word_size):
                score = score + dist[s[j]][s2[j]]
            if score < tscore_sim:
                continue
            #Otherwise add the score to the
            word_dictionary[s2].append(current_word_pos)






def test_dict():
#build the dictionary with some test data
    d = "pwnaaplhnfgedflqpyvqlqqnfsasdlevnleatreshahfstpqalelflnysvtp"
    td = create_word_dict(d)

    string_wrap = textwrap.fill(str(td), 80)
    print(string_wrap)
    #print(td)
    print("Dist contains for QSYV: lookup", td.get("QSYV"))

#test_dict()