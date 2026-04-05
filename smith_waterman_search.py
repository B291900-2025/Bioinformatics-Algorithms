########################################################################
#         Search using SW                                              #
#         Super SLOW compared to BLAST                                 #
########################################################################
import smith_waterman_p as SW
import time
import os
from operator import itemgetter
import  process_fasta_file as pff
from print_logger import *
import calc_bit_and_evalues as cbe
import programme_settings

programme_settings.read()

aligntime = 0
max_scores = int(programme_settings.settings["SWSEARCH"]["max_sw_scores"])

query_sequence = programme_settings.settings["DEFAULT"]["query_sequence"]
query_sequence =query_sequence.upper()
database = programme_settings.settings["DEFAULT"]["database"]

def processSW(myline_database):
    global aligntime
    global query_sequence
    global settings
    q_seq=query_sequence.upper()
    db_entry = myline_database.upper()
    t3 = time.time()
    res = SW.perform_smith_waterman(q_seq,db_entry,False,False)
    t4 = time.time()
    aligntime =aligntime +(t4-t3)
    return res

#just buildit

t0 = time.time()

print =logger

res = pff.process_fasta_file(database, processSW, max_scores)

#printing final output by iterating the scores list
results_count = 1
t1 = time.time()
#change the order to print the best result first!
res.sort(key=itemgetter(2),reverse = True)

print("###############################################################")
print("SW Alignment Results for SwissProt Database Search")
print("Max Number of Results to return: ", max_scores)
print("Matrix: Blosum%s, Gap Weight: %s" % (programme_settings.settings["DEFAULT"]["blosum"],programme_settings.settings["DEFAULT"]["seq_gap"]))
print("Query Sequence:")
print(query_sequence)
print("###############################################################\n")

#grap the parameters to use for the stats
lib_size = int(programme_settings.settings["DEFAULT"]["current_library_size"])
K = float(programme_settings.settings["DEFAULT"]["K"])
lam = float(programme_settings.settings["DEFAULT"]["lam"])
len = len(query_sequence)

#write the raw data to a file
#res ==i  is (current_header, myline, score,fileindex)
with open("logs\SWsearch.csv", 'w') as outputfile:
    outputfile.write("Name\tSWScore\tBitScore\tE-value\tSeqIndex\n")

    for i in res:
        outputfile.write("\""+str(i[0])[1:]+"\""+"\t")
        outputfile.write(str(i[2]) + "\t")
        outputfile.write(cbe.get_bit_score_s(i[2], K, lam) + "\t")
        outputfile.write(cbe.get_expect_s(i[2], K, lam) + "\t")
        outputfile.write(str(i[3]) + os.linesep)

#bit of code duplication here -but needed a different formatting
print("Name\tSWScore\tBitScore\tE-value\tSeqIndex")

for i in res:
    print("\"" + str(i[0])[1:40] + "...\"" + "\t", end='')
    print(str(i[2]) + "\t", end='')
    print(cbe.get_bit_score_s(i[2], K, lam)+ "\t", end='')
    print(cbe.get_expect_s(i[2], K, lam)+ "\t", end= '')
    print(str(i[3]) + "\n")

for val in res:
    print(os.linesep + "Result number : ", results_count)
    print("Sequence File Index: ",val[3])
    print("Name: ", val[0])
    #print(val[1])
    print("Score: %s, BitScore %s, e-value %s" % (val[2],cbe.get_bit_score_s(val[2], K, lam),cbe.get_expect_s(val[2], K, lam)) )
    results_count = results_count +1
    SW.perform_smith_waterman(query_sequence,val[1].upper(),False,True)


print("Time taken: {} secs".format(t1-t0))
print("Percentage time taken SW: ",aligntime/(t1-t0)*100,"%")

print("~~~~~~~~~Finished~~~~~~~~~")

