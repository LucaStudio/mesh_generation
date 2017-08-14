# Autogenerated with SMOP 
from smop.core import *
# 

    
@function
def sphericalcoordinates(r=None,phi=None,theta=None,*args,**kwargs):
    varargin = sphericalcoordinates.varargin
    nargin = sphericalcoordinates.nargin

    ##
#        Project: Fluid - structure interaction on deformable surfaces
#         Author: Luca Di Stasio
#    Institution: ETH Zrich
#                 Institute for Building Materials
# Research group: Computational Physics for Engineering Materials
#        Version: 0.1
#  Creation date: August 6th, 2014
#    Last update: August 6th, 2014
    
    #    Description: 
#          Input: 
#         Output:
    
    ##
    
    coords=matlabarray(cat(multiply(multiply(r,sin(phi)),cos(theta)),multiply(multiply(r,sin(phi)),sin(theta)),multiply(r,cos(phi))))
    return coords