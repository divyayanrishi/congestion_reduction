
# coding: utf-8

# In[1]:

import networkx as nx 
from networkx import Graph 
import random
import matplotlib.pyplot as plt
def newman_square(n,m,p,seed=None):
    print("Executing ",p)
    G=nx.grid_graph(dim=[n,m],periodic=False)  
    pos=nx.spring_layout(G)
    nlist=list(G.nodes())
    nlist2=nlist
    elist=20*list(G.edges())
    c=0   
    for j in range(20):
        i=1
        for u in nlist:
            if random.random()<p:
                w=random.choice(list(nlist))
                while w==u or G.has_edge(u,w):
                    w=random.choice(list(nlist))
                    if G.degree(u)>=n-1:
                        break    
                G.add_edge(u,w) 
                elist.append((u,w))   
        c=c+nx.diameter(G)
        #print("%d" % (len(elist)/20),end=' ')
    d=c/20
    nx.write_gexf(G,"D:\Research\Complex Networks\\newman_data_old.gexf")                
        #nx.draw_networkx_nodes(G,pos,nodelist=nlist,with_labels=True)
        #nx.draw_networkx_edges(G,pos,edgelist=elist)
        #plt.show()
    return d,len(elist)/20
#Main function
list_prob=[]
for k in range(0,95,1):
    list_prob.append(k/100)
list_diam=[]
list_edge=[]
a=0
b=0
for i in range(0,len(list_prob)):
    a,b=newman_square(10,10,list_prob[i],seed='None')
    list_diam.append(a)
    list_edge.append(b)
print(list_prob)
print(list_diam)    
print(list_edge)
plt.plot(list_prob,list_diam,scalex=True,scaley=True)
plt.xlabel('Wiring Probability')
plt.ylabel('Network Diameter')
plt.show()
plt.plot(list_prob,list_edge,scalex=True,scaley=True)
plt.xlabel('Wiring Probability')
plt.ylabel('No. of edges')
plt.show()


# In[ ]:



