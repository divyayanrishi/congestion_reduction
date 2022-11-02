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


G = nx.erdos_renyi_graph(30,0.2,seed=None)

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
total_cap=5700
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

#Model for optimization using Pyomo
model=ConcreteModel() 
key_list=list(edge_enum.keys())
edge_list=list(edge_enum.values())
model.links=Set(initialize=key_list) 
lb={} 
for j in key_list:
    lb[j]=pmf_maxima(edge_enum[j],lamd,q)   
#print(lb)    
ub={}
for j in key_list:
    ub[j]=pmf_length(edge_enum[j],lamd,q)     
#print('\n',ub)    
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
model.obj=Objective(rule=objrule)      
opt=SolverFactory('ipopt') 
res=opt.solve(model)

#Display optimized edge capacities
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