# Autogenerated with SMOP 
from smop.core import *
# 

    
@function
def R(lattice=None,*args,**kwargs):
    varargin = R.varargin
    nargin = R.nargin

    ##
#        Project: Fluid-Structure Interaction on Deformable Surfaces
#         Author: Luca Di Stasio
#    Institution: ETH Zrich
#                 Institute for Building Materials
# Research group: Computational Physics for Engineering Materials
#        Version: 0.1
#  Creation date: July 24th, 2014
#    Last update: July 24th, 2014
    
    #          Input: meshed computational domain compdomain
#         Output: mesh in the physical domain
    
    ##
    
    Rvec=zeros(size(lattice,1),1)
    return Rvec