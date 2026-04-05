############################################################
#  Takes a FASTA file                                      #
#  Takes a callback function                               #
#  Call function like processSW(query_sequence, myline)    #
#  Find the SW optimum sequence                            #
############################################################


from operator import itemgetter
import time

bestscore = 0
bestheader = ""
bestbody = ""
aligntime = 0

global aligner_timer_secs

#OK add tuple for storing scores
res = []

#add two seqs
def add_score(nsssc, max_scores):
    global res

    res.append(nsssc)
    res.sort(key=itemgetter(2),reverse = False) #ascending scores
    if len(res) > max_scores:
        res.pop(0)


def process_fasta_file(fastafile, processfunction, max_scores,ats =None):

    f1 = open(fastafile, 'r')
    f1lines = f1.readlines()
    read_body = False
    count = 1
    myline = ""
    fileindex =0

    global bestscore
    global bestheader
    global bestbody
    global  aligntime

    current_header=""

    for line in f1lines:
        line = line.rstrip()
        if line.startswith('>'):
            count = count + 1
            if count %1000==0:
                print("#",count)

            if read_body == True:

                # we are ending the processing of a seq
                t3=time.time()
                score = processfunction( myline)
                t4 =time.time()
                if ats is not None:
                    # store the score
                    ats["Fastatime"] += t4 - t3

                if score > bestscore:
                    bs = (current_header, myline, score,fileindex)
                    add_score(bs, max_scores)
                    fileindex+=1
                read_body = False

        if read_body == False:
            myline = ""
            current_header = line
            read_body = True
        else:
            if len(line) > 0:
                myline = myline + line
    f1.close()

    if  read_body == True:
        t3 =time.time()
        score = processfunction(myline)
        t4 = time.time()
        if ats is not None:
            # store the score
            ats["Fastatime"] += t4 - t3



        if score > bestscore:
            bs = (current_header, myline, score,fileindex)
            add_score(bs, max_scores)
            fileindex+=1

        print("Processing completed: ",count)

    return res