print("**********************************************************************************")
print("*                                                                                *")
print("*                         BLAST BioPython Search Example                         *")
print("*                                   Beta version!                                *")
print("**********************************************************************************")


#This code uses biopython
from Bio.Blast import NCBIWWW
import programme_settings as ps
import time

#dbs to use ftp://ftp.ncbi.nlm.nih.gov/blast/db/
#docs https://biopython.org/docs/latest/Tutorial/
########################################################################
#                        BioPython BLAST Example                       #
########################################################################
ps.read()
#make a toy file to test
myseq = ps.settings["DEFAULT"]["query_sequence"]
with open("logs/toyfasta.fasta", "w") as testfile:
    testfile.write(">testfile\n")
    testfile.write(myseq +"\n")
print("Writing Testfile for BioPython...")

print("Starting BLAST....")
t0 =time.time()
record = open("logs/toyfasta.fasta",encoding='utf8').read()
results_stream = NCBIWWW.qblast("blastp","swissprot",record,format_type = "Text")
data = results_stream.read()
print("Completed reading....")
print(data)
t1 = time.time()
print("Elapsed time: %dsecs"%(t1-t0))

#help(NCBIWWW.qblast)