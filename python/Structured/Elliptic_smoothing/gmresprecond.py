# Autogenerated with SMOP 
from smop.core import *
# 

    
@function
def gmresprecond(A=None,b=None,invM=None,x0=None,m=None,*args,**kwargs):
    varargin = gmresprecond.varargin
    nargin = gmresprecond.nargin

    ##
#        Project: Fluid - structure interaction on deformable surfaces
#         Author: Luca Di Stasio
#    Institution: ETH Zrich
#                 Institute for Building Materials
# Research group: Computational Physics for Engineering Materials
#        Version: 0.1
#  Creation date: May 21st, 2014
#    Last update: May 21st, 2014
    
    #    Description: 
#          Input: 
#         Output:
    
    ##
    
    N=size(A,1)
    v=zeros(N,m + 1)
    h=zeros(m + 1,m)
    r0=b - dot(A,x0)
    beta=sqrt(sum(r0 ** 2))
    v[:,1]=r0 / beta
    j=1
    loop=1
    while loop and j <= m:

        wj=dot(A,(dot(invM,v[:,j])))
        for i in arange(1,j).reshape(-1):
            h[i,j]=dot(wj.T,v[:,i])
            wj=wj - dot(h[i,j],v[:,i])
        h[j + 1,j]=sqrt(sum(wj ** 2))
        if h[j + 1,j] == 0:
            loop=0
        else:
            v[:,j + 1]=wj / h[j + 1,j]
            if j == m:
                loop=0
            else:
                j=j + 1

    
    if j != m:
        v=v[:,1:j]
        h=h[1:j + 1,1:j]
    
    r=zeros(j + 1,1)
    r[1,1]=beta
    for k in arange(1,j).reshape(-1):
        P=eye(j + 1)
        s=h[k + 1,k] / (sqrt(h[k + 1,k] ** 2 + h[k,k] ** 2))
        c=h[k,k] / (sqrt(h[k + 1,k] ** 2 + h[k,k] ** 2))
        P[k,k]=c
        P[k + 1,k + 1]=c
        P[k,k + 1]=s
        P[k + 1,k]=- s
        h=dot(P,h)
        r=dot(P,r)
    
    y=backsolve(h[1:m,:],r[1:m,:])
    x=x0 + dot(v[:,1:m],y)
    return x