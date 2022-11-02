#!/usr/bin/env python
# coding: utf-8

# In[172]:

import csv 
import scipy as sc
from scipy.optimize import minimize 
import networkx as nx
import math as m
import numpy as np
import random
import matplotlib.pyplot as plt
import math as m
from coopr.pyomo import *
from pyomo.opt import SolverFactory
import pyomo.environ as pyo 


# In[82]:


G = nx.watts_strogatz_graph(30,4,0.6)

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
        if sum(listed[j:]) <= clip: #edge_clips[e] maybe  
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
total_cap=5600
net_cap=total_cap-(total_cap % nx.number_of_edges(G)) 
print(net_cap)
lamd = 6
#q=0.6 
a=0 
for e in G.edges():
    edge_enum[a]=e
    a+=1

#here we find loss of an edge e of the graph 
def find_loss(e,lmd,q,k):  #here k is capacity of edge e not the packet no.
    edge_pmf=find_final_dist(e,lmd,q) 
    l_e=sum(edge_pmf[k+1:])
    return l_e    

#initial guess on capacity vector
x0={}
for i in range(a):
    x0[i]=m.floor(net_cap/nx.number_of_edges(G))       
#x0[a-1]=m.floor(net_cap/nx.number_of_edges(G)) + total_cap % nx.number_of_edges(G)  
c=total_cap % nx.number_of_edges(G)
while c>0:
    n=random.randrange(0,a)
    x0[n]+=1
    c-=1


#Define the optimization model in Pyomo
model=ConcreteModel() 
key_list=list(edge_enum.keys())
edge_list=list(edge_enum.values())
model.links=Set(initialize=key_list) 
lb={} 
for j in key_list:
    lb[j]=pmf_maxima(edge_enum[j],lamd,q)   
print(lb)    
ub={}
for j in key_list:
    ub[j]=pmf_length(edge_enum[j],lamd,q)     
print('\n',ub)    
def fb(model,i):
    return (lb[i],ub[i]) 
model.vars=Var(model.links,initialize=x0,domain=PositiveIntegers,bounds=fb) 
def conrule(model): 
    s=0
    for e in model.links:
       s+= model.vars[e]  
    return(s==net_cap)
model.cons=Constraint(rule=conrule)
def objrule(model):   
    s=0
    for e in model.links:
        s+=find_loss(edge_enum[e],lamd,q,pyo.value(model.vars[e])) 
    return s 
model.obj=Objective(rule=objrule,sense=minimize)      
opt=SolverFactory('ipopt') 
res=opt.solve(model)

#Display results
p=0
cap_vec=[]
cent_vec=[]
for e in model.links:
    print(" ",pyo.value(model.vars[e]))
    cap_vec.append(pyo.value(model.vars[e]))
    cent_vec.append(nx.edge_betweenness_centrality(G)[edge_enum[e]]) 
    p+=pyo.value(model.vars[e])
print(p)  
print(np.corrcoef(cap_vec,cent_vec))   

'''
def get_E_prime(E_set):
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


# In[174]:


capacity_edges



# In[92]:


list_of_centrality = nx.edge_betweenness_centrality(G)
use_this = sorted(list_of_centrality, key = list_of_centrality.get, reverse = True)
edge = use_this[0]
X = np.linspace(0,179,180)
Y = find_final_dist(edge, lamd,q)[0:180]
plt.plot(X,Y,'bo')
plt.axvline(x = capacity_edges[edge], ymax = 0.4, ls = '--')
plt.fill_between(X[capacity_edges[edge]:179],Y[capacity_edges[edge]:179], facecolor = 'yellow')
plt.xlabel("# of packets")
plt.ylabel("Probability")
plt.text(150,0.005, '(c)', fontsize = 20)
plt.savefig('Max_evaluate025_cap.pdf')
plt.show()


# In[32]:


#simulation part 


# In[93]:


def break_paths_into_edges(path):
    final_list = []
    for k in range(0,len(path)-1):
        final_list.append((path[k],path[k+1]))
    return final_list


# In[183]:


itr_val = 90
dict_of_dicts = {}

for i in range(0,itr_val):
    test_dict_in_itr = {}
    
    #initialize this to be zeros
    for e in get_E_prime(G.edges()):
        test_dict_in_itr[e] = 0
    
    pkts_and_destinations = {}
    for j in G.nodes():
        j_copy = [j]
        dest = random.choice(list(set(G.nodes())-set(j_copy)))
        if(np.random.uniform(0,1) <= q):
            no_of_packets = np.random.poisson(lamd)
        else:
            no_of_packets = 0
        pkts_and_destinations[j] = (no_of_packets, dest, nx.single_source_shortest_path(G,j)) 
    
    for k in range(0,len(G.nodes())):
        for m in range(0,len(G.nodes())):
            if k != m:
                for e in break_paths_into_edges(pkts_and_destinations[k][2][m]):
                    test_dict_in_itr[e] += pkts_and_destinations[k][0]
    
    dict_of_dicts[i] = test_dict_in_itr
        


# In[184]:


#now_check_theory_with_simulation:
fraction_of_success = {}
for e in get_E_prime(G.edges()):
    fraction_of_success[e] = 0 

for j in range(0,itr_val):
    for e in get_E_prime(G.edges()):
        if capacity_edges[e] >= dict_of_dicts[j][e]: 
            fraction_of_success[e] += 1/itr_val  
        else:
            fraction_of_success[e] += 0


# In[185]:


fraction_of_success


# In[188]:


a_s = list(fraction_of_success.values())
weights = np.ones_like(a_s)/float(len(a_s))
plt.hist(a_s,weights=weights)
plt.xlabel("C")
plt.ylabel("Fraction of links free for C fraction of time windows")
plt.text(0.25,0.36,'(d)',fontsize = '20')
plt.savefig('Histogram_06_90.pdf')
plt.show()


# In[154]:


GMOSvC = {}
GMOSvC_intervals = {}


# In[189]:


a = np.histogram(list(fraction_of_success.values()))[0].tolist()
b = np.histogram(list(fraction_of_success.values()))[1].tolist()
a = [x / sum(a) for x in a]
GMOSvC[(q,itr_val)] = a
GMOSvC_intervals[(q,itr_val)] = b


# In[213]:


plt.plot(GMOSvC_intervals[(1,30)][0:10],GMOSvC[(1,30)], 'ro' ,label = "(q,#TF) = (1,30)")
plt.plot(GMOSvC_intervals[(0.6,30)][0:10],GMOSvC[(0.6,30)],'go' ,label = "(q,#TF) = (0.6,30)")
plt.xlabel('C')
plt.ylabel('$g(\lambda,q)$')
plt.legend(loc = 2)
plt.text(0.4,0.5,'(e)',fontsize = '20')
plt.savefig('GMOSvC_FIG1.pdf')
plt.show()


# In[214]:


plt.plot(GMOSvC_intervals[(1,90)][0:10],GMOSvC[(1,90)], 'ro' ,label = "(q,#TF) = (1,90)")
plt.plot(GMOSvC_intervals[(0.6,90)][0:10],GMOSvC[(0.6,90)], 'go',label = "(q,#TF) = (0.6,90)")
plt.xlabel('C')
plt.ylabel('$g(\lambda,q)$')
plt.legend(loc = 2)
plt.text(0.4,0.5,'(f)',fontsize = '20')
plt.savefig('GMOSvC_FIG2.pdf')
plt.show()


# In[22]:


# To check if variance changes with changing of q. FOR CONSTANT CLIP
Q_list = np.linspace(0.1,1,10)
Avg_SD_list = []
for take in Q_list:
    q = take
    a = 0
    for e in get_E_prime(G.edges()):
        temp_list = find_final_dist(e,lamd,q) 
        a = a + np.sqrt(np.var(temp_list))
    a = a/len(G.edges())
    Avg_SD_list.append(a)


# In[29]:


plt.plot(Q_list,Avg_SD_list,'bo')
plt.xlabel('q')
plt.ylabel('Standard Deviation of PMF')
plt.text(0.85,0.03,'(f)',fontsize = 15)
plt.savefig('SDvsq.pdf')
plt.show()
'''

# In[ ]:




