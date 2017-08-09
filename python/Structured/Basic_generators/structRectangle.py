# Autogenerated with SMOP 
from smop.core import *
# 

    
@function
def structRectangle(x0=None,y0=None,lx=None,ly=None,Nx=None,Ny=None,*args,**kwargs):
    varargin = structRectangle.varargin
    nargin = structRectangle.nargin

    ##
#==============================================================================
# Copyright (c) 2016 Universit de Lorraine & Lule tekniska universitet
# Author: Luca Di Stasio <luca.distasio@gmail.com>
#                        <luca.distasio@ingpec.eu>
    
    # Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# 
# 
# Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
# Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in
# the documentation and/or other materials provided with the distribution
# Neither the name of the Universit de Lorraine or Lule tekniska universitet
# nor the names of its contributors may be used to endorse or promote products
# derived from this software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#==============================================================================
    
    #  DESCRIPTION
#  
#  A function to generate meshed rectangles
    
    #  Input: x0 - scalar - x-coordinate of center
#         y0 - scalar - y-coordinate of center
#         lx - scalar - Side length in x-direction
#         ly - scalar - Side length in y-direction
#         Nx - scalar - Number of ELEMENTS in x-direction
#         Ny - scalar - Number of ELEMENTS in y-direction
#  Output: mesh - (Nx+1)*(Ny+1) x 2 matrix - mesh nodes ordered through helical
#          indexing
    
    ##
    
    mesh=zeros(dot((Nx + 1),(Ny + 1)),2)
    deltax=lx / Nx
    deltay=ly / Ny
    xs=(arange((x0 - dot(0.5,lx)),(x0 + dot(0.5,lx)),deltax)).T
    for j in arange(1,Ny + 1).reshape(-1):
        mesh[dot((j - 1),(Nx + 1)) + 1:dot(j,(Nx + 1)),1]=xs
        mesh[dot((j - 1),(Nx + 1)) + 1:dot(j,(Nx + 1)),2]=dot(((y0 - dot(0.5,ly)) + dot(deltay,(j - 1))),ones((Nx + 1),1))
    
    return mesh