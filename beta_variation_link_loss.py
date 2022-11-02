import csv 
import scipy as sc
from scipy.optimize import minimize 
import networkx as nx
import math as m
import numpy as np
import random
import matplotlib.pyplot as plt
import math as m
from random import choice 
import pandas as pd 

G=nx.erdos_renyi_graph(30,0.4)  
print(nx.number_of_edges(G))
nx.draw(G)
plt.savefig(r'D:\Research\Complex Networks\cap_cent_scatterPlots&graphs\erd_ren_beta(30,0.4)_1.png') 

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

q=0.6    

def gen_till_length(ld,f,length,q):
    dist_returned = []
    for k in range(0,length):
        dist_returned.append(pssn_modified(ld,k,f,q))
    return dist_returned
def find_final_dist(edg,lmd,q):
    final = []
    zero_removed = remove_zero(find_list_for_edge(edg))
    result = np.convolve(gen_till_length(lmd,zero_removed[0],60,q),gen_till_length(lmd,zero_removed[1],60,q))
    for i in range(2,len(zero_removed)):
        result = np.convolve(result,gen_till_length(lmd,zero_removed[i],60,q))
    return result
list_of_centrality = nx.edge_betweenness_centrality(G)
use_this = sorted(list_of_centrality, key = list_of_centrality.get, reverse = True)
t1 = find_final_dist(use_this[0],4,0.6)[0:180]
t2 = find_final_dist(use_this[1],4,0.6)[0:180]
t3 = find_final_dist(use_this[2],4,0.6)[0:180]

def find_clip(listed,clip):
    for j in range(0,len(listed)):
        if sum(listed[j:]) <= clip: #edge_clips[e] maybe  
            return j
            break 

capacity_edges = {}
edge_enum={}
clip = 0.15
net_cap=4180
lamd = 4
a=0
for e in G.edges(): 
    edge_enum[a]=e
    a+=1

#here we find loss of an edge e of the graph 
def find_loss(e,lmd,q,k): 
    edge_pmf=find_final_dist(e,lmd,q) 
    l_e=sum(edge_pmf[(m.floor(k)+1):len(edge_pmf)])
    return l_e       
edge_caps={}
list_cent=[] 
avg_edg_cap=[]
cap_sum=0
'''for i in range(10,91,10):
    cap_sum=0
    i_list=[]
    for e in G.edges():
        list_t=find_final_dist(e,lamd,q)
        k=find_clip(list_t,(1-(i/100)))
        i_list.append(k)    
    edge_caps[i/100]=i_list     
prob_list=list(edge_caps.keys()) 

for i in prob_list.keys():
    cap_sum=sum(prob_list[i])
    avg_edg_cap.append(cap_sum) 
'''
for e in G.edges():
    list_cent.append(nx.edge_betweenness_centrality(G)[e])
cent_sum=sum(list_cent[:]) 

edge_caps['edge_cent']=list_cent 
#edge_caps['avg_cap']=avg_edg_cap  
#df=pd.DataFrame(edge_caps,columns=list(edge_caps.keys()))    
#print(df) 
#df.to_excel(r'D:\Research\Complex Networks\Bar_alb(30,4)q0.6l4.xlsx',sheet_name='Sheet1',index=False)  
beta=list(range(1,10)) 
beta_loss=[] 
for i in range(len(beta)): 
    s=0
    for e in G.edges(): 
        c=m.floor(1+beta[i]*100*nx.edge_betweenness_centrality(G)[e]/cent_sum) 
        s=s+find_loss(e,lamd,q,c)   
    beta_loss.append(s) 
print(beta_loss)
plt.clf()     
plt.scatter(beta,beta_loss)
plt.xlabel('beta')
plt.ylabel('loss value')  
plt.show() 
beta_dict={}
beta_dict['beta']=beta
beta_dict['loss value']=beta_loss
df=pd.DataFrame(beta_dict,columns=list(beta_dict.keys()))
df.to_excel(r'D:\Research\Complex Networks\Dist_sheets\beta_loss.xlsx',sheet_name='Sheet1')
plt.xlabel('beta')
plt.ylabel('loss value')  
plt.show()      

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