import numpy as np
from Bio import SeqIO
import random
from collections import defaultdict
from Levenshtein import distance
import plotly.plotly as py
import plotly.figure_factory as ff
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as sch
import scipy.spatial as scs
import math

# Making the dataset
lists = []

for record in SeqIO.parse("DNASequences.fasta", "fasta"):
  lists.append(record)

# Loading the proximity matrix
f = open("prox_first.bin", 'rb')
d = np.load(f)
f.close()
print("Distances Matrix Loaded")

# Functions
def min_in_matrix(X):
    min = math.inf
    min_i = 0
    min_j = 0
    for i in range(len(X)):
        for j in range(len(X[i])):
            if i != j and min > X[i][j]:
                min = X[i][j]
                min_i = i
                min_j = j
    return min_i,min_j,min

def check_inf(X):
    for i in range(len(X)):
        for j in range(len(X[i])):
            if X[i][j]!=math.inf:
                return 1
    return 0

def marked(Z):
    for i in range(len(Z)):
        if Z[i]==0:
            return 1
    return 0

def complete_link(X):
    Z = [] # New Linkage matrix
    iterator=0
    clusters = {} # keeps track of the clusters
    count_clusters = {}
    for i in range(len(X)):
        count_clusters[i] = 1
    while check_inf(X):
        min_i,min_j,min = min_in_matrix(X)
        # Adding to the Linkage matrix
        if min_i in clusters and min_j in clusters:
            Z.append([clusters[min_i],clusters[min_j],min,count_clusters[clusters[min_i]]+count_clusters[clusters[min_j]]])
        elif min_i in clusters:
            a = count_clusters[clusters[min_i]]+1
            Z.append([clusters[min_i],min_j,min,a])
        elif min_j in clusters:
            a = 1+count_clusters[clusters[min_j]]
            Z.append([min_i,clusters[min_j],min,a])
        else:
            Z.append([min_i,min_j,min,2])

        # updating the clusters
        clusters[min_j] = len(X)+iterator
        iterator+=1
        count_clusters[clusters[min_j]] = 1
        # Updating the clusters counts
        if min_i in clusters and min_j in clusters:
            count_clusters[clusters[min_j]] = count_clusters[clusters[min_i]] + count_clusters[clusters[min_j]]
        elif min_i in clusters:
            count_clusters[clusters[min_j]] = count_clusters[clusters[min_i]] + 1
        elif min_j in clusters:
            count_clusters[clusters[min_j]] = count_clusters[clusters[min_j]] + 1
        else:
            count_clusters[clusters[min_j]] = 2
        # Updating the rest of the array
        for i in range(len(X)):
            if i!=min_i and i!=min_j:
                a = X[min_j][i]
                b = X[min_i][i]
                if a >= b:
                    X[min_j][i] = X[i][min_j] =  a
                else:
                    X[min_j][i] = X[i][min_j] = b
        # removing one of the data points
        for i in  range(len(X)):
            X[min_i][i]=X[i][min_i] = math.inf

    return X,Z

for i in range(len(d)):
    d[i][i] = math.inf

# Complete_Link
d_complete, Zdash_complete = complete_link(d)
Zdash_complete = np.array(Zdash_complete,dtype='float64')
sch.dendrogram(Zdash_complete)
plt.savefig('complete.png')
plt.show()
print("Complete Generated")
