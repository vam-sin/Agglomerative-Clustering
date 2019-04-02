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

def average_link(X):
    iterator=0
    Z = [] # New Linkage matrix
    clusters = {} # keeps track of the clusters
    count_clusters = {}
    for j in range(len(X)):
        count_clusters[j] = 1
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
        # Updating the rest of the array
        for i in range(len(X)):
            if i !=min_i and i!=min_j:
                dist = 0
                a = 0 # Total points in both clusters combined
                # Find a
                if min_i in clusters and min_j in clusters:
                    a = count_clusters[clusters[min_i]] + count_clusters[clusters[min_j]]
                elif min_i in clusters:
                    a = count_clusters[clusters[min_i]] + count_clusters[min_j]
                elif min_j in clusters:
                    a = count_clusters[min_i] + count_clusters[clusters[min_j]]
                else:
                    a = 2
                # Find total distance
                # if check_inf_row(X[i]):
                if min_i in clusters:
                    dist+=X[i][min_i]*count_clusters[clusters[min_i]]
                else:
                    dist+=X[i][min_i]
                if min_j in clusters:
                    dist+=X[min_j][i]*count_clusters[clusters[min_j]]
                else:
                    dist+=X[min_j][i]
                dist = dist/float(a)
                X[min_j][i] = X[i][min_j] =  float(dist)
        # Updating the clusters counts and clusters
        clusters[min_j] = len(X)+iterator
        iterator+=1
        count_clusters[clusters[min_j]] = 1
        if min_i in clusters and min_j in clusters:
            count_clusters[clusters[min_j]] = count_clusters[clusters[min_i]] + count_clusters[clusters[min_j]]
        elif min_i in clusters:
            count_clusters[clusters[min_j]] = count_clusters[clusters[min_i]] + 1
        elif min_j in clusters:
            count_clusters[clusters[min_j]] = count_clusters[clusters[min_j]] + 1
        else:
            count_clusters[clusters[min_j]] = 2
        # removing one of the data points
        for i in  range(len(X)):
            X[min_i][i]=X[i][min_i] = math.inf
    return X,Z

for i in range(len(d)):
    d[i][i] = math.inf

# average_link
d_average, Zdash_average = average_link(d)
Zdash_average = np.array(Zdash_average,dtype='float64')
sch.dendrogram(Zdash_average)
plt.savefig('average.png')
plt.show()
print("Average Generated")
