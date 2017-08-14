# Autogenerated with SMOP 
from smop.core import *
# 

    
@function
def generatecells2D(Nx=None,Ny=None,*args,**kwargs):
    varargin = generatecells2D.varargin
    nargin = generatecells2D.nargin

    ##
#        Project: Fluid - structure interaction on deformable surfaces
#         Author: Luca Di Stasio
#    Institution: ETH Zrich
#                 Institute for Building Materials
# Research group: Computational Physics for Engineering Materials
#        Version: 0.1
#  Creation date: July 30th, 2014
#    Last update: July 30th, 2014
    
    #    Description: 
#          Input: 
#         Output:
    
    ##
    
    cells=zeros(dot((Nx - 1),(Ny - 1)),5)
    for j in arange(1,Ny - 1).reshape(-1):
        for k in arange(1,Nx - 1).reshape(-1):
            nodeindex1=k + dot((j - 1),Nx)
            nodeindex2=k + dot((j - 1),Nx) + Nx
            nodeindex3=k + dot((j - 1),Nx) + Nx + 1
            nodeindex4=k + dot((j - 1),Nx) + 1
            cellindex=k + dot((j - 1),(Nx - 1))
            cells[cellindex,1]=4
            cells[cellindex,2]=nodeindex1
            cells[cellindex,3]=nodeindex4
            cells[cellindex,4]=nodeindex3
            cells[cellindex,5]=nodeindex2
    
    return cells