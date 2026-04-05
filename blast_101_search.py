print("**********************************************************************************")
print("*                                                                                *")
print("*                     Blast 101 Search Basic Blast Programme                     *")
print("*                                   Beta version!                                *")
print("**********************************************************************************")
# Simon Tomlinson Bioinformatics Algorithms 2025
#  You need to install the required packages eg blosum

import smith_waterman_p as SW
import time
from operator import itemgetter
import  process_fasta_file as pff
import create_seq_word_dict as tdict
from collections import defaultdict
import blosum as bl
from print_logger import *
import programme_settings
import calc_bit_and_evalues as cbe

#These are used for timing of various events
aligner_timer_secs ={}
programme_settings.read()

# matrix for blosum default = blosum62
ddist = bl.BLOSUM(int(programme_settings.settings["BLAST"]["blosum"]))

max_extension_length =int(programme_settings.settings["BLAST"]["max_extension_length"])

qsequence = programme_settings.settings["DEFAULT"]["query_sequence"]
qsequence =qsequence.upper()

#this can now become a BLAST object that ranks the seqs based upon a simple score
#SW alignment only needs to be performed at the end
#UP TO Here

query_sequence = tdict.create_word_dict(qsequence)

aligntime = 0
max_scores = int(programme_settings.settings["BLAST"]["max_scores"])
max_alignments = int(programme_settings.settings["BLAST"]["max_alignments"])
word_size =int(programme_settings.settings["BLAST"]["word_size"])
lib_seqs =int(programme_settings.settings["DEFAULT"]["current_library_seq_count"])
lib_seq_length = int(programme_settings.settings["DEFAULT"]["current_library_size"])

def extend_diagonal(pos_s0_s1,s0,s1):
 #   SW.dist is the current distance matrix
    global dist
    global word_size
    #after that step back or forward depending upon which is the best
    left_pos = pos_s0_s1
    right_pos = pos_s0_s1

#extend until the edge is reached or currentscore==0
# first word residues are from the matching triplet so must exist
# they are also to the right of the co-ordinate

    p1 = right_pos[0]
    p2 = right_pos[1]

    cscore = 0
    for i in range(0,word_size):
        right_pos = (p1 +i,p2+i)
        r1 = s0[right_pos[0]]
        r2 = s1[right_pos[1]]
        cscore += ddist[r1][r2]

    bescore = cscore

#after that adjust the position pointers to extend the bestpos by one step
#if scores for left right extension are equal add both

    extending = True

    while extending:
        tmp_right_pos=(right_pos[0] +1,right_pos[1]+1)
        tmp_left_pos =(left_pos[0] -1,left_pos[1]-1)

        if tmp_right_pos[0] >=len(s0) or tmp_right_pos[1] >=len(s1):
            tmp_right_pos = None
        if tmp_left_pos[0] <0 or tmp_left_pos[1] <0:
            tmp_left_pos = None
        if tmp_left_pos  is None and tmp_right_pos is None:
            extending = False #no further extensions possible
            continue
        if tmp_left_pos is not None:
            tmp_left_score = SW.dist[s0[tmp_left_pos[0]]] [s1[tmp_left_pos[1]]]
        else:
            tmp_left_score =None
        if tmp_right_pos is not None:
            tmp_right_score = SW.dist[s0[tmp_right_pos[0]]] [s1[tmp_right_pos[1]]]
        else:
            tmp_right_score =None

        #block the final extension of the
        if tmp_right_score is None:
            cscore += tmp_left_score
            left_pos = tmp_left_pos
        elif tmp_left_score is None:
            cscore += tmp_right_score
            right_pos = tmp_right_pos
        else:
            if tmp_right_score >= tmp_left_score:  #slight bias here
                cscore+= tmp_right_score
                right_pos = tmp_right_pos
            else:
                cscore += tmp_left_score
                left_pos = tmp_left_pos
        if cscore <= int(programme_settings.settings["BLAST"]["min_extension_score"]):
            extending = False
        if cscore > bescore:
            bescore = cscore

    return bescore


#Performs the core part of the BLAST search
def process_blast(myline_database):
    #word_size =int(programme_settings.settings["BLAST"]["word_size"]) #word size for searching
    #max_extension_length =int(programme_settings.settings["BLAST"]["max_extension_length"])
    count_matches =0

    global query_sequence
    global qsequence
    global aligntime


    t3 = time.time()

   #add some BLAST  stuff here
    res = []
    for i in range(0, len(myline_database)):

        s = myline_database[i:i + word_size]
        if len(s) % word_size != 0:
            continue

        if len(query_sequence[s]) > 0: #OK found a match
            count_matches +=1

            for j in range(0,len(query_sequence[s])):
                res.append((query_sequence[s][j],i))
    #this gives a pair of matches  j indexed is query match i is target match
    #sort based upon query position
    # get smallest first!

    #sanity check
    if count_matches <2:
        return 0
    res.sort(key=itemgetter(0,1), reverse=False)

    #There is a need to avoid synonym matches which is why this code runs here!!

    res_store =defaultdict(list)
    cnt =0

    res_store[cnt].append(res[0])

    for el in range(1,len(res)):
        if res[el][0]==res_store[cnt][0]:
            res_store[cnt].append(res[el])
        else:
            cnt+=1
            res_store[cnt].append(res[el])

    #then search first A  and last B in contiguous order
    #we now perform an all-against all cfn using res_store
    #
    length = len(res_store)

    #can't have any pairs with 1 match
    if length < 2:
        return 0

    bestscore = 0

    #print(length)
    for i in range(0, length):
        for j in range(0, length):

            if i == j:
                continue
            len_i =len(res_store[i])
            len_j = len(res_store[j])

            for k in range(0,len_i):
                for m in range(0,len_j):

                    #keep if in the right order
                    if (res_store[i][k][0] > res_store[j][m][0]) or (res_store[i][k][1] > res_store[j][m][1]):
                        continue
                    #keep if on same diagonal
                    if (res_store[j][m][0] - res_store[i][k][0]) != (res_store[j][m][1] - res_store[i][k][1]):
                        continue

                    #keep less than 40
                    if (res_store[j][m][0] - res_store[i][k][0]) > max_extension_length:
                        continue

                    #score largest contiguous section in the correct order
                    #score the segment (without gaps) using the scoring matrix
                    #return the score...
                    currentscore = extend_diagonal(res_store[i][k],qsequence,myline_database)
                    if currentscore > bestscore:
                        bestscore = currentscore

                    currentscore = extend_diagonal(res_store[j][k],qsequence,myline_database)
                    if currentscore > bestscore:
                        bestscore = currentscore



    t4 = time.time()
    aligntime =aligntime +(t4-t3)
    return bestscore

#process the file
def process_fasta_file():
    global aligner_timer_secs

    res = pff.process_fasta_file(programme_settings.settings["DEFAULT"]["database"], process_blast, max_scores,aligner_timer_secs)


    #change the order to print the best result first!
    res.sort(key=itemgetter(2),reverse = True)

    print(" HSP Scores: ", len(res)) #tmp test
    #Processing HSPs with SW alignment
    return res

#Split here...
def print_final_results(res):
    global t1, t2, t3, t4, t0
    global qsequence
    #create a new object to store the final results -res 2
    res2 =[]

    for val in res:
        swscore  = SW.perform_smith_waterman(qsequence, val[1].upper(), False, False)
        res2.append((val[0],val[1],val[2],val[3],swscore))


    res2.sort(key=itemgetter(4),reverse = True)

    #Completed ... Outputting Final Blast Scores
    print = logger

    print("###############################################################")
    print("BLAST Hits  for SwissProt Database Search")
    print("Max Number of Results to return: ", max_scores)
    print("Matrix: Blosum%s, Gap Weight: %s" % (programme_settings.settings["DEFAULT"]["blosum"],programme_settings.settings["DEFAULT"]["seq_gap"]) )
    print("Query Sequence:")
    print(qsequence)
    print("###############################################################")

    #grap the parameters to use for the stats
    k = float(programme_settings.settings["DEFAULT"]["k"])
    scale = float(programme_settings.settings["DEFAULT"]["lam"])
    #len = len(qsequence)


    #write the raw data to a file
    with open("logs\BLsearch.csv", 'w') as outputfile:
        outputfile.write("Name\tSWScore\tExscore\tBitScore\tExpect\tSeqIndex"+os.linesep)

        for i in res2:
            outputfile.write("\""+str(i[0])[1:]+"\""+"\t")
            outputfile.write(str(i[4]) + "\t")
            outputfile.write(str(i[2]) + "\t")
            #bitscore and expect
            outputfile.write(cbe.get_bit_score_s(i[4],k,scale) + "\t")
            outputfile.write(cbe.get_expect_s(i[4],k,scale) + "\t")
            outputfile.write(str(i[3]) + os.linesep)

    #bit of code duplication here -but needed a different formatting
    print("Name\tSWScore\tExScore\tBitScore\tExpect\tSeqIndex" + os.linesep)

    for i in res2:
        print("\"" + str(i[0])[1:40] + "...\"" + "\t", end='')
        print(str(i[4]) + "\t", end='')
        print(str(i[2]) + "\t", end='' )
        print(cbe.get_bit_score_s(i[4], k,scale) + "\t",end='')
        print(cbe.get_expect_s(i[4], k,scale) + "\t",end='')
        print(str(i[3]) + os.linesep)
        #raw_score,k, scale,number_seqs

    loop =0
    for val in res2:
        print(os.linesep + "Result number : ", loop +1, end ="")
        print(" Sequence File Index: ", val[3])
        print(" Name: ", val[0], end ="")
        print(" Score: ", val[2])
        print("SWScore: ",val[4], end ="")
        print(" Bit Score: ",cbe.get_bit_score_s(i[4], k,scale),end ="")
        print(" Expect: " ,cbe.get_expect_s(i[4],k,scale))
        #results_count = results_count +1
        SW.perform_smith_waterman(qsequence,val[1].upper(),False,True)
        loop+=1
        if loop >=max_alignments:
            break

    # print out the final counts
    print_timer()

    print("~~~~~~~~~Finished~~~~~~~~~")

def blast101_run():
    init_print_timer()

    t1 = time.time()
    #process the FASTA files
    res =process_fasta_file()

    t2 = time.time()
    aligner_timer_secs["Elapsedtime"] += t2-t1

    # print out the final results
    print_final_results(res)


#this just makes a neat printout of the timer
def print_timer():

    global aligner_timer_secs

    print("***********FINAL TIMING RESULTS ************")
    print("\nElement\tTime(secs)")

    for k,v in aligner_timer_secs.items():
        vv = float(v)

        print(f"{k:}\t{vv:.2f} secs")

    if aligner_timer_secs["Elapsedtime"] !=0:
        print("Percentage time taken BLAST searching: %.2f %%" % ((aligner_timer_secs["Fastatime"] / aligner_timer_secs["Elapsedtime"]) * 100))
        print("********************************************")

def init_print_timer():
    aligner_timer_secs["Fastatime"] =0
    aligner_timer_secs["Elapsedtime"] = 0


blast101_run()