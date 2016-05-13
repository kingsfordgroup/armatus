#!/usr/env python

import numpy as np
import seaborn as sb
import matplotlib.pyplot as plt
import argparse
import math
from scipy.sparse import coo_matrix



def plotall(datamat,domains1,domains2,bounds,legendname1,legendname2,outputname):
    """ Show heatmap of Hi-C data along with any domain sets given
    :param datamat: Hi-C data matrix as numpy array
    :param domains1: nx2 list of domains (optional, use [] to just see heatmap of Hi-C matrix)
    :param domains2: nx2 list of domains (optional, use [] if no second set of domains)
    :param bounds: (x,y) to view only bins between x and y (optional, use () to see entire chromosome)
    :param legendname1: legend label for first set of domains
    :param legendname2: legend label for second set of domains
    :param outputname: filename of image to be saved (optional - use [] to view instead of save)
    :return: either show image (if outputname == []) or save image as outputname
    """

    if bounds == (): # plot full Hi-C matrix with all TADs
        logdata = np.ma.log(datamat)
        logdata = logdata.filled(0)
        labelspacing = int(math.floor(round(len(logdata),-int(math.floor(math.log10(len(logdata)))))/10))
        ax = sb.heatmap(logdata,cbar=False,xticklabels=labelspacing,yticklabels=labelspacing)


        if domains1 != []:
            for interval in domains1: # plot outline of each domain
                plt.plot((interval[0]-1,interval[1]),(len(logdata)-interval[0]+1,len(logdata)-interval[0]+1),'g')
                plt.plot((interval[1],interval[1]),(len(logdata)-interval[0]+1,len(logdata)-interval[1]),'g')
                dom1Artist = plt.Line2D((0,1),(0,0), color='green', linestyle='solid')
        if domains2 != []:
            for interval in domains2:
                plt.plot((interval[0]-1,interval[1]),(len(logdata)-interval[1],len(logdata)-interval[1]),'b')
                plt.plot((interval[0]-1,interval[0]-1),(len(logdata)-interval[0]+1,len(logdata)-interval[1]),'b')
                dom2Artist = plt.Line2D((0,1),(0,0), color='blue', linestyle='solid')
    else: # show only range of matrix between bounds
        logdata = np.ma.log(datamat[bounds[0]:bounds[1],bounds[0]:bounds[1]])
        logdata = logdata.filled(0)
        labelspacing = int(math.floor(round(len(logdata),-int(math.floor(math.log10(len(logdata)))))/10))
        ax = sb.heatmap(logdata,cbar=False,xticklabels=labelspacing,yticklabels=labelspacing)

        if domains1 != []:
            for interval in domains1:
                if interval[0] >= bounds[0] and interval[1] <= bounds[1]:
                    interval -= bounds[0]
                    plt.plot((interval[0]-1,interval[1]),(len(logdata)-interval[0]+1,len(logdata)-interval[0]+1),'g')
                    plt.plot((interval[1],interval[1]),(len(logdata)-interval[0]+1,len(logdata)-interval[1]),'g')
            dom1Artist = plt.Line2D((0,1),(0,0), color='green', linestyle='solid')
        if domains2 != []:
            for interval in domains2:
                if interval[0] >= bounds[0] and interval[1] <= bounds[1]:
                    interval -= bounds[0]
                    plt.plot((interval[0]-1,interval[1]),(len(logdata)-interval[1],len(logdata)-interval[1]),'b')
                    plt.plot((interval[0]-1,interval[0]-1),(len(logdata)-interval[0]+1,len(logdata)-interval[1]),'b')
            dom2Artist = plt.Line2D((0,1),(0,0), color='blue', linestyle='solid')

    if legendname1 and legendname2:
        legend = ax.legend([dom1Artist,dom2Artist], [legendname1, legendname2],frameon = 1)
        legendframe = legend.get_frame()
        legendframe.set_facecolor('white')
        legendframe.set_edgecolor('black')
    elif legendname1:
        legend = ax.legend([dom1Artist],[legendname1])
        legendframe = legend.get_frame()
        legendframe.set_facecolor('white')
        legendframe.set_edgecolor('black')

    # save image to file if filename was given, .png is default if no extension given
    if outputname:
        plt.savefig(outputname)
    else: # just display image
        plt.show()



def parseRaoFormat(datamat,res):
    """ turn sparse Rao data format into dense matrix for heatmap
    :param datamat: Hi-C data in sparse format as numpy array (n x 3)
    :param res: resolution of data
    :return: dense Hi-C data matrix
    """

    datamat[:,0:2] = datamat[:,0:2]/res
    datamat = coo_matrix((datamat[:,2], (datamat[:,0],datamat[:,1]) ))
    datamat = datamat.todense()
    if datamat.shape[0] > datamat.shape[1]:
        # add column(s) of zeros to make square matrix
        ncols = datamat.shape[0] - datamat.shape[1]
        sqmat = np.zeros((datamat.shape[0],datamat.shape[0]))
        sqmat[:,:-1*ncols] = datamat
        datamat = sqmat
    elif datamat.shape[1] > datamat.shape[0]:
        # add row(s) of zeros to make square matrix
        nrows = datamat.shape[1] - datamat.shape[0]
        sqmat = np.zeros((datamat.shape[1],datamat.shape[1]))
        sqmat[:-1*nrows,:] = datamat
        datamat = sqmat
    datamat = datamat + np.transpose(datamat) - np.diagonal(datamat)*np.identity(len(datamat))
    return datamat



def main(datafile, res, domainfile1, domainfile2, domainres1, domainres2, windowbounds, legendname1, legendname2, outputname):

    datamat = np.genfromtxt(datafile,delimiter='\t')
    if datamat.shape[1] == 3: # Rao format
        datamat = parseRaoFormat(datamat, res)
    else: # remove any NaNs
        datamat = datamat[~np.isnan(datamat)]
        datamat = np.reshape(datamat,(np.sqrt(len(datamat)), np.sqrt(len(datamat))))
    if domainfile1:
        domains1 = np.genfromtxt(domainfile1,delimiter='\t')
        domains1 = domains1[~np.isnan(domains1)]/domainres1
        domains1 = np.reshape(domains1,(len(domains1)/2,2))
    else:
        domains1 = []
    if domainfile2:
        domains2 = np.genfromtxt(domainfile2,delimiter='\t')
        domains2 = domains2[~np.isnan(domains2)]/domainres2
        domains2 = np.reshape(domains2,(len(domains2)/2,2))
    else:
        domains2 = []
    if windowbounds:
        bounds = (int(windowbounds[0]),int(windowbounds[1]))
    else:
        bounds = ()
    if not legendname1: # make filenames legend entry, if none is given
        legendname1 = domainfile1
    if not legendname2:
        legendname2 = domainfile2

    plotall(datamat,domains1,domains2,bounds,legendname1,legendname2,outputname)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Domain visualization tool for Hi-C data.')
    parser.add_argument('-i', metavar='inputFile', help='raw Hi-C data filename (tab-delimited text file of Hi-C data or Rao format)')
    parser.add_argument('-r', metavar='Resolution', default = [], type=int, help='Hi-C Resolution (only needed if using Rao data format)')
    parser.add_argument('-b', metavar=('startBound','endBound'), nargs=2, default=(), help='Bounds for viewing window (optional)')
    parser.add_argument('-d1', metavar='domainFile1', default=[], help='TAD file')
    parser.add_argument('-d2', metavar='domainFile2', default=[], help='second TAD file (optional)')
    parser.add_argument('-dr1', metavar='domainResolution1', type=int, default=1, help='Resolution of domains in domainFile1')
    parser.add_argument('-dr2', metavar='domainResolution2', type=int, default=1, help='Resolution of domains in domainFile2')
    parser.add_argument('-l1', metavar='legendName1', default=[], type=str, help='Legend name for first set of domains')
    parser.add_argument('-l2', metavar='legendName2', default=[], type=str, help='Legend name for second set of domains')
    parser.add_argument('-o', metavar='outputFile', default=[], type=str, help='Filename for saved image file')
    args = parser.parse_args()

    main(args.i, args.r, args.d1, args.d2, args.dr1, args.dr2, args.b, args.l1, args.l2, args.o)

