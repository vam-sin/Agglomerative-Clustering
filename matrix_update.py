import numpy as np
from Levenshtein import distance
import scipy.cluster.hierarchy as sch
import scipy.spatial as scs
import matplotlib.pyplot as plt

def min_in_matrix(X):
    min = 1000000
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
            if X[i][j]!=1000000:
                return 1
    return 0

def marked(Z):
    for i in range(len(Z)):
        if Z[i]==0:
            return 1
    return 0

def single_link(X):
    Z = [] # New Linkage matrix
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
        clusters[min_j] = len(X)+min_j-1
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
                if a <= b:
                    X[min_j][i] = X[i][min_j] =  a
                else:
                    X[min_j][i] = X[i][min_j] = b
        # removing one of the data points
        for i in  range(len(X)):
            X[min_i][i]=X[i][min_i] = 1000000

    return X,Z

def complete_link(X):
    Z = [] # New Linkage matrix
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
        clusters[min_j] = len(X)+min_j-1
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
            X[min_i][i]=X[i][min_i] = 1000000

    return X,Z

def average_link(X):
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
        print(X)
        # Updating the rest of the array
        for i in range(len(X)):
            if i !=min_i and i!=min_j:
                dist = 0
                a = 0 # Total points in both clusters combined
                # Find a
                print(min_i,min_j,i)
                print(clusters)
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
                print(a,dist)
                dist = dist/float(a)
                print(dist)
                X[min_j][i] = X[i][min_j] =  float(dist)
                print(X)
        # Updating the clusters counts
        # updating the clusters
        clusters[min_j] = len(X)+min_j-1
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
            X[min_i][i]=X[i][min_i] = 1000000
    return X,Z

trial = np.array([[0, 1, 2, 3],
                 [1, 0, 4, 5],
                 [2, 4, 0, 6],
                 [3, 5, 6, 0]])

trial = np.array(trial,dtype='float64')

for i in range(len(trial)):
    trial[i][i] = 1000000

# Single_Link
# trial_single, Zdash_single = single_link(trial)
# Zdash_single = np.array(Zdash_single,dtype='float64')
# print(Zdash_single)

# complete_link
# trial_complete, Zdash_complete = complete_link(trial)
# Zdash_complete = np.array(Zdash_complete,dtype='float64')
# print(Zdash_complete)

# Average link
trial_average, Zdash_average = average_link(trial)
Zdash_average = np.array(Zdash_average,dtype='float64')
print(Zdash_average)

Z = np.array([[0.,1,1,2],
             [2.,4,2,3],
             [3.,5,3,4]])
print(Z)

sch.dendrogram(Z)
plt.show()

# sch.dendrogram(Zdash_single)
# plt.show()

# sch.dendrogram(Zdash_complete)
# plt.show()

sch.dendrogram(Zdash_average)
plt.show()
