#!/usr/bin/env python
# coding: utf-8

# In[172]:

import csv 
import scipy as sc
import networkx as nx
import math as m
import numpy as np
import random
import matplotlib.pyplot as plt
import math as m
import mystic as my  
from mystic.constraints import integers
from mystic.solvers import diffev
from mystic.math import almostEqual 
from mystic.ensemble import BuckshotSolver
from mystic.ensemble import LatticeSolver
from mystic.ensemble import buckshot,lattice
from mystic.differential_evolution import DifferentialEvolutionSolver2
from mystic.differential_evolution import diffev2
from mystic.solvers import fmin_powell

# In[82]:


G = nx.barabasi_albert_graph(30,4)

# In[83]:


def pssn(lam,k):
    return (lam**k)*np.exp(-lam)/m.factorial(k)

def pssn_change(lam,k,q):
    if k ==0: 
        return (1-pssn(lam,0))*(1-q) + pssn(lam,0)
    else:
        return pssn(lam,k)*q 

def pssn_modified(lambd,k,f,q):
    if k == 0:
        return f*pssn_change(lambd,0,q) + 1 - f
    else:
        return f*pssn_change(lambd,k,q)
    
def edge_in_path(edge,path):
    if(edge[0] in path):
        index = path.index(edge[0])
        if(index + 1 != len(path) and path[index + 1] == edge[1]):
            return True 
        else:
            return False
    else:
        return False
def find_list_for_edge(edg): 
    listed = []
    for i in G.nodes():
        for j in G.nodes():
            temp = list(nx.all_shortest_paths(G,i,j))
            count = 0   
            if len(temp) == 1:
                if edge_in_path(edg,temp[0]):
                    count = 1
            else:
                for k in temp:
                    if edge_in_path(edg,k):
                        count += 1/len(temp)
            listed.append(count)
    return listed

def remove_zero(listed):
    result = []
    for i in listed:
        if i != 0:
            result.append(i)
    return result


# In[170]:


q = 0.6


# In[85]:


def gen_till_length(ld,f,length,q):
    dist_returned = []
    for k in range(0,length):
        dist_returned.append(pssn_modified(ld,k,f,q))
    return dist_returned


# In[86]:


def find_final_dist(edg,lmd,q):
    final = []
    zero_removed = remove_zero(find_list_for_edge(edg))
    result = np.convolve(gen_till_length(lmd,zero_removed[0],60,q),gen_till_length(lmd,zero_removed[1],60,q))
    for i in range(2,len(zero_removed)):
        result = np.convolve(result,gen_till_length(lmd,zero_removed[i],60,q))
    return result
    

# In[87]:

 
list_of_centrality = nx.edge_betweenness_centrality(G)
use_this = sorted(list_of_centrality, key = list_of_centrality.get, reverse = True)
t1 = find_final_dist(use_this[0],4,0.6)[0:180]
t2 = find_final_dist(use_this[1],4,0.6)[0:180]
t3 = find_final_dist(use_this[2],4,0.6)[0:180]


# In[88]:


plt.plot(np.linspace(0,179,180),t1, 'ro' ,label = "maximum betw")
plt.plot(np.linspace(0,179,180),t2, 'bo' ,label = "2nd maximum betw")
plt.plot(np.linspace(0,179,180),t3, 'go', label = "3rd maximum betw")
#plt.title("Packet Flow Probability Distribution at links with highest betweenness")
plt.xlabel("# of packets")
plt.ylabel("Probability")
plt.legend(loc = 1)
plt.text(150,0.005, '(d)', fontsize = 20)
plt.savefig('Final_dist_link_25.pdf')
plt.show()


# In[89]:


def find_clip(listed,clip):
    for j in range(0,len(listed)):
        if sum(listed[j:]) <= clip: 
            return j 
            break 

def pmf_maxima(edg,lmd,q):
    edge_pmf=find_final_dist(edg,lmd,q)
    for j in range(1,pmf_length(edg,lmd,q)):
        if(edge_pmf[j]>edge_pmf[j+1]):
            return j   

def pmf_length(edg,lmd,q):
    num=0
    for i in find_final_dist(edg,lmd,q):
      num+=1
    return num       

# In[173]:


#here we find the link capacities of all the edges
capacity_edges = {}
edge_enum={}
clip = 0.15 
total_cap=1000 
net_cap=total_cap-(total_cap % nx.number_of_edges(G))
lamd = 4
a=0 
for e in G.edges():
    edge_enum[a]=e
    a+=1

#here we find loss of an edge e of the graph 
def find_loss(e,lmd,q,k):  #here k is capacity of edge e not the no. of packets 
    edge_pmf=find_final_dist(e,lmd,q) 
    l_e=sum(edge_pmf[k+1:])
    return l_e   

def objective(x):
    s=0 
    for i in range(len(x)):
        s+=find_loss(edge_enum[i],lamd,q,x[i]) 
    return s 

added=lambda x: [i for i in x]        
cons=lambda x: my.constraints.impose_sum(total_cap,added(x)) 
round=np.round 
 
constr=my.constraints.and_(cons,round) 

bounds=[(0,0)]*a
for i in range(nx.number_of_edges(G)):
    c=pmf_length(edge_enum[i],lamd,q)
    bounds[i]=(0,c-1) 
x0=[0]*nx.number_of_edges(G)
for i in range(a):
    x0[i]=m.floor(net_cap/nx.number_of_edges(G)) 
c=total_cap % nx.number_of_edges(G) 
while c>0:
    n=random.randrange(0,a)
    x0[n]+=1
    c-=1              
print(x0) 
s=objective(x0)
print(s)
r=constr(x0)
print(r)
#result=diffev2(objective,x0=x0,bounds=bounds,constraints=constr,npop=104,gtol=50,maxiter=100000,disp=True,full_output=True)
res=fmin_powell(objective,x0=x0,constraints=constr,bounds=bounds) 
print(res)    
'''def get_E_prime(E_set):
    E_prime = []
    for e in E_set: 
        E_prime.append(e)
        E_prime.append((e[1],e[0])) 
    return E_prime

for e in get_E_prime(G.edges()):
    temp_list = find_final_dist(e,lamd,q)
    clipper = find_clip(temp_list,clip) 
    if(clipper != 0):
        capacity_edges[e] =  clipper
'''