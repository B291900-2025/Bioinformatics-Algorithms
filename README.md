# Blast101_code
This archive contains the main code for the Bioinformatics Algorithms BLAST class that implements BLAST-like and Smith Waterman search in Python.

This is code is designed for teaching and uses a limited subset of Python- that which we teach in introductory Python courses.  

The class ICA extends this code to add tests so that the validity of the code can be checked.  This is 'beta code' designed for teaching, not for use outside the class!

The code is designed to run inside the PyCharm IDE.   The code has multiple 'entry point' for the code execution.  To run a single part of the code, select the file of interest
from within PyCharm, right click and choose Run.  This will run the selected code.

biopython_eg.py- this code runs an example BLAST job to the NCBI service and retrieves and displays the result.  This is how you might choose to run BLAST from the NCBI
using your own script.

blast_101_search.py- this code runs the local BLAST 101 search using the parameters specified in the settings.ini file

smith_waterman_search.py - performs a SW search of the current database- NB this code can be slow for large databases!!

settings.ini - contains the settings used by the programmes.

logs- this contains logs of the various runs performed including some example SW search results (that take quite a long time to compute)

Other files/folders can also be run or viewed as required- these are just the main entry points for the  Blast101 code...

Simon Tomlinson March 2026





