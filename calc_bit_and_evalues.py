###############################################################################
#   Stats calculations                                                        #
#   Notice here this code uses Pandas to read the data                        #
#   This code then uses Seaborn for some axis (and hidden histogram)          #
#   Numpy for the numerical part                                              #
#   Scipy for the scientific computing part                                   #
#   Matplotlib for plotting the final results                                 #
#   All the tools mentioned in the Python Primer week 0!                      #
###############################################################################

import math
import programme_settings as ps
#see https://stackoverflow.com/questions/72985838/scipy-gumbel-fit-does-not-fit-what-is-the-correct-shape-of-the-array-datafra
#see https://www.ncbi.nlm.nih.gov/BLAST/tutorial/Altschul-1.html
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import scipy.stats as ss


ps.read()

def build_fit():
    print("**************************************************************************")
    print("*               Fit Gumbel Distribution from Simulation Data             *")
    print("**************************************************************************")


    data = pd.read_csv(ps.settings["BUILD_EXPECT"]["random_file"])
    values = data.iloc[:, 0]
    sns.histplot(values, kde=True, stat='density')
    #plt.show()
    x = np.linspace(0,max(values),101)
    #bins = np.histogram(values, bins=x, density=True)
    fit = ss.gumbel_r.fit(values.to_numpy())

    sigma =fit[0] #loc AKA sigma K
    mu =fit[1] #scale lambda

    if ps.settings["BUILD_EXPECT"]["update_settings"] == "True":
        ps.settings["DEFAULT"]["K"] = str(sigma)
        ps.settings["DEFAULT"]["lam"] = str(mu)
        ps.write()

    dist = ss.gumbel_r(loc=sigma,scale=mu)
    plt.plot(x, dist.pdf(x),color ="red")

    if ps.settings["BUILD_EXPECT"]["show_fit"]=="True":
        plt.show()


#K and lam are from the Gumbel fit
#raw score comes from the search
#https://www.cs.cmu.edu/~durand/03-711/2015/Lectures/DBSearchingAndBlast-statistics.pdf
#score from https://academic.oup.com/bioinformatics/article/40/3/btae097/7629128
def get_bit_score(raw_score, k, scale):
    #Calculate bit score from E value

    #total count of bases in the database
    n = int(ps.settings["DEFAULT"]["current_library_size"])

    #count of seqs
    m = int(ps.settings["DEFAULT"]["current_library_seq_count"])

    e = get_expect(raw_score,k,scale)

    #average seq size
    avg = n/m

    # bs =log2(space/evalue)
    #problem that e is sometimes exactly zero
    if(e ==0):
        bs = float('nan')
        return bs

    bs = math.log2((n * avg)/e)

    return bs

#here just use the p-value from the dist to fit to the data
#
def get_p_value(raw_score, k, scale):
    p =1-ss.gumbel_r.cdf(raw_score, loc=k, scale=scale)
    return p

def get_expect(raw_score,k, scale):
    p = get_p_value(raw_score, k, scale)
    if p==0:
        return 0

    e = 1-math.exp(-1*p)
    return e

#just write a nice string representation for easy output
def get_bit_score_s(raw_score, k, scale, precision =0):

    s =get_bit_score(raw_score, k, scale)
    f = "{:."+str(precision)+"f}"
    s = f.format(s)
    return s

#just write a nice string representation for easy output
def get_expect_s(raw_score,k, scale,precision = 5):
    s = get_expect(raw_score,k,scale)

    if s<0.0001:
        f = "{:." + str(precision) + "E}"
        s = f.format(s)
    else:
        f = "{:." + str(precision) + "f}"
        s = f.format(s)
    return s




#test method based upon current parameters
def test():
    #from nanog search & uniprot_sprot database
    raw_score =300
    k = 6.23037
    scale = 29.01160

    print("Expect score (Nanog):", get_expect(raw_score,k,scale))
    print("Bit score (Nanog):", get_bit_score(raw_score,k,scale))
    print("Bit score (Nanog) 0dp",get_bit_score_s(raw_score,k,scale))
    print("Expect (Nanog) 5dp", get_expect_s(raw_score,k,scale,precision =5))

    print("Expect score (low):", get_expect(40,k,scale))
    print("Bit score (low):", get_bit_score(40,k,scale))
    print("Expect (low) 5dp", get_expect_s(40, k, scale,precision= 5))

#run the code
#build_fit()
#test()