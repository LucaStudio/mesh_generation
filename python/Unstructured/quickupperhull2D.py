# Autogenerated with SMOP 
from smop.core import *
# 

    
@function
def quickupperhull2D(p=None,printflag=None,*args,**kwargs):
    varargin = quickupperhull2D.varargin
    nargin = quickupperhull2D.nargin

    ##
#        Project: Fluid-Structure Interaction on Deformable Surfaces
#         Author: Luca Di Stasio
#    Institution: ETH Zrich
#                 Institute for Building Materials
# Research group: Computational Physics for Engineering Materials
#        Version: 0.1
#  Creation date: April 23rd, 2014
#    Last update: April 23rd, 2014
    
    #          Input: N x 3 vector points in plane (indices + coordinates)
#         Output: M x 1 vector of indeces of points belonging to the upper hull
    
    ##
    
    N=size(p,1)
    A,Ia=mintournamenttree(p,true,2,nargout=2)
    B,Ib=maxtournamenttree(p,true,2,nargout=2)
    pold=copy(p)
    p=matlabarray([])
    for i in arange(1,N).reshape(-1):
        if i != Ia and i != Ib:
            p=matlabarray(cat([p],[pold[i,:]]))
    
    p,Ip=mintournamenttree(p,false,2,nargout=2)
    clear('Ip')
    #------------------------------ UPPER HULL --------------------------------
    
    Cindex=0
    ABCarea=0
    tocheckupper=matlabarray([])
    for i in arange(1,N - 2).reshape(-1):
        area=orient2D(A[1,2:end()].T,B[1,2:end()].T,p[i,2:end()].T)
        if area > 0 and area > ABCarea:
            if Cindex != 0:
                tocheckupper=matlabarray(cat([tocheckupper],[Cindex]))
            Cindex=copy(i)
            ABCarea=copy(area)
        else:
            if area > 0:
                tocheckupper=matlabarray(cat([tocheckupper],[i]))
    
    C=p[Cindex,:]
    if printflag:
        figure()
        plot(p[:,2],p[:,3],'*')
        hold('on')
        for i in arange(1,size(tocheckupper,1)).reshape(-1):
            plot(p[tocheckupper[i,1],2],p[tocheckupper[i,1],3],'gd','LineWidth',2)
            hold('on')
        plot(A[:,2],A[:,3],'rd','LineWidth',2)
        hold('on')
        plot(B[:,2],B[:,3],'rd','LineWidth',2)
        hold('on')
        plot(C[:,2],C[:,3],'kd','LineWidth',2)
        hold('on')
        plot(cat([A[:,2]],[B[:,2]]),cat([A[:,3]],[B[:,3]]),'--k','LineWidth',2)
        hold('on')
        grid('on')
        pause
    
    tocheckupperAC=matlabarray([])
    tocheckupperCB=matlabarray([])
    for i in arange(1,size(tocheckupper,1)).reshape(-1):
        area=orient2D(A[1,2:end()].T,B[1,2:end()].T,p[tocheckupper[i],2:end()].T)
        areaAC=orient2D(A[1,2:end()].T,C[1,2:end()].T,p[tocheckupper[i],2:end()].T)
        if (areaAC > 0 and area > 0) or (orient2D(C[1,2:end()].T,B[1,2:end()].T,p[tocheckupper[i],2:end()].T) > 0 and area > 0):
            if areaAC > 0:
                tocheckupperAC=matlabarray(cat([tocheckupperAC],[tocheckupper[i]]))
            else:
                tocheckupperCB=matlabarray(cat([tocheckupperCB],[tocheckupper[i]]))
    
    if printflag:
        figure()
        plot(p[:,2],p[:,3],'*')
        hold('on')
        for i in arange(1,size(tocheckupperAC,1)).reshape(-1):
            plot(p[tocheckupperAC[i,1],2],p[tocheckupperAC[i,1],3],'gd','LineWidth',2)
            hold('on')
        plot(A[:,2],A[:,3],'rd','LineWidth',2)
        hold('on')
        plot(B[:,2],B[:,3],'rd','LineWidth',2)
        hold('on')
        plot(C[:,2],C[:,3],'kd','LineWidth',2)
        hold('on')
        plot(cat([A[:,2]],[B[:,2]]),cat([A[:,3]],[B[:,3]]),'--k','LineWidth',2)
        hold('on')
        plot(cat([A[:,2]],[C[:,2]]),cat([A[:,3]],[C[:,3]]),'--k','LineWidth',2)
        hold('on')
        grid('on')
        fprintf('There are %d points to check\\n',size(tocheckupperAC,1))
        pause
    
    tocheckupperAC=matlabarray(cat(tocheckupperAC,dot(A[1,1],ones(length(tocheckupperAC),1)),dot(A[1,2],ones(length(tocheckupperAC),1)),dot(A[1,3],ones(length(tocheckupperAC),1)),dot(C[1,1],ones(length(tocheckupperAC),1)),dot(C[1,2],ones(length(tocheckupperAC),1)),dot(C[1,3],ones(length(tocheckupperAC),1))))
    indexhullAC=matlabarray([])
    while logical_not(isempty(tocheckupperAC)):

        Xindex=0
        ACXarea=0
        localA=tocheckupperAC[1,2:4]
        localC=tocheckupperAC[1,5:7]
        currentindex=0
        for i in arange(1,size(tocheckupperAC,1)).reshape(-1):
            area=orient2D(localA[1,2:end()].T,localC[1,2:end()].T,p[tocheckupperAC[i],2:end()].T)
            if area > 0 and area > ACXarea:
                currentindex=copy(i)
                Xindex=tocheckupperAC[i,1]
                ACXarea=copy(area)
        tocheckupperAC=matlabarray(cat([tocheckupperAC[1:currentindex - 1,:]],[tocheckupperAC[currentindex + 1:end(),:]]))
        tempcheckupperAC=copy(tocheckupperAC)
        tocheckupperAC=matlabarray([])
        if Xindex != 0:
            X=p[Xindex,:]
            indexhullAC=matlabarray(cat([indexhullAC],[Xindex]))
            if printflag:
                plot(cat([X[:,2]],[localA[:,2]]),cat([X[:,3]],[localA[:,3]]),'--k','LineWidth',2)
                hold('on')
                plot(cat([localC[:,2]],[X[:,2]]),cat([localC[:,3]],[X[:,3]]),'--k','LineWidth',2)
                hold('on')
            for i in arange(1,size(tempcheckupperAC,1)).reshape(-1):
                if tempcheckupperAC[i,2] == localA[1,1] and tempcheckupperAC[i,5] == localC[1,1]:
                    if orient2D(localA[1,2:end()].T,X[1,2:end()].T,p[tempcheckupperAC[i,1],2:end()].T) > 0 or orient2D(X[1,2:end()].T,localC[1,2:end()].T,p[tempcheckupperAC[i,1],2:end()].T) > 0:
                        if orient2D(localA[1,2:end()].T,X[1,2:end()].T,p[tempcheckupperAC[i,1],2:end()].T) > 0:
                            tocheckupperAC=matlabarray(cat([tocheckupperAC],[tempcheckupperAC[i,1],tempcheckupperAC[i,2:4],p[Xindex,:]]))
                        else:
                            tocheckupperAC=matlabarray(cat([tocheckupperAC],[tempcheckupperAC[i,1],p[Xindex,:],tempcheckupperAC[i,5:7]]))
                else:
                    tocheckupperAC=matlabarray(cat([tocheckupperAC],[tempcheckupperAC[i,:]]))
        if printflag:
            plot(cat([localA[:,2]],[localC[:,2]]),cat([localA[:,3]],[localC[:,3]]),'--k','LineWidth',2)
            hold('on')
            for i in arange(1,size(tocheckupperAC,1)).reshape(-1):
                plot(p[tocheckupperAC[i,1],2],p[tocheckupperAC[i,1],3],'gd','LineWidth',2)
                hold('on')
            grid('on')
            xlabel('x')
            ylabel('y')
            title('Convex hull')
            fprintf('There are %d points to check\\n',size(tocheckupperAC,1))
            pause

    
    if printflag:
        figure()
        plot(p[:,2],p[:,3],'*')
        hold('on')
        for i in arange(1,size(tocheckupperCB,1)).reshape(-1):
            plot(p[tocheckupperCB[i,1],2],p[tocheckupperCB[i,1],3],'gd','LineWidth',2)
            hold('on')
        plot(A[:,2],A[:,3],'rd','LineWidth',2)
        hold('on')
        plot(B[:,2],B[:,3],'rd','LineWidth',2)
        hold('on')
        plot(C[:,2],C[:,3],'kd','LineWidth',2)
        hold('on')
        plot(cat([A[:,2]],[B[:,2]]),cat([A[:,3]],[B[:,3]]),'--k','LineWidth',2)
        hold('on')
        plot(cat([C[:,2]],[B[:,2]]),cat([C[:,3]],[B[:,3]]),'--k','LineWidth',2)
        hold('on')
        grid('on')
        fprintf('There are %d points to check\\n',size(tocheckupperCB,1))
        pause
    
    tocheckupperCB=matlabarray(cat(tocheckupperCB,dot(C[1,1],ones(length(tocheckupperCB),1)),dot(C[1,2],ones(length(tocheckupperCB),1)),dot(C[1,3],ones(length(tocheckupperCB),1)),dot(B[1,1],ones(length(tocheckupperCB),1)),dot(B[1,2],ones(length(tocheckupperCB),1)),dot(B[1,3],ones(length(tocheckupperCB),1))))
    indexhullCB=matlabarray([])
    while logical_not(isempty(tocheckupperCB)):

        Xindex=0
        CBXarea=0
        localC=tocheckupperCB[1,2:4]
        localB=tocheckupperCB[1,5:7]
        currentindex=0
        for i in arange(1,size(tocheckupperCB,1)).reshape(-1):
            area=orient2D(localC[1,2:end()].T,localB[1,2:end()].T,p[tocheckupperCB[i],2:end()].T)
            if area > 0 and area > CBXarea:
                currentindex=copy(i)
                Xindex=tocheckupperCB[i,1]
                CBXarea=copy(area)
        tocheckupperCB=matlabarray(cat([tocheckupperCB[1:currentindex - 1,:]],[tocheckupperCB[currentindex + 1:end(),:]]))
        tempcheckupperCB=copy(tocheckupperCB)
        tocheckupperCB=matlabarray([])
        if Xindex != 0:
            X=p[Xindex,:]
            indexhullCB=matlabarray(cat([indexhullCB],[Xindex]))
            if printflag:
                plot(cat([X[:,2]],[localC[:,2]]),cat([X[:,3]],[localC[:,3]]),'--k','LineWidth',2)
                hold('on')
                plot(cat([localB[:,2]],[X[:,2]]),cat([localB[:,3]],[X[:,3]]),'--k','LineWidth',2)
                hold('on')
            for i in arange(1,size(tempcheckupperCB,1)).reshape(-1):
                if tempcheckupperCB[i,2] == localC[1,1] and tempcheckupperCB[i,5] == localB[1,1]:
                    if orient2D(localC[1,2:end()].T,X[1,2:end()].T,p[tempcheckupperCB[i,1],2:end()].T) > 0 or orient2D(X[1,2:end()].T,localB[1,2:end()].T,p[tempcheckupperCB[i,1],2:end()].T) > 0:
                        if orient2D(localC[1,2:end()].T,X[1,2:end()].T,p[tempcheckupperCB[i],2:end()].T) > 0:
                            tocheckupperCB=matlabarray(cat([tocheckupperCB],[tempcheckupperCB[i,1],tempcheckupperCB[i,2:4],p[Xindex,:]]))
                        else:
                            tocheckupperCB=matlabarray(cat([tocheckupperCB],[tempcheckupperCB[i,1],p[Xindex,:],tempcheckupperCB[i,5:7]]))
                else:
                    tocheckupperCB=matlabarray(cat([tocheckupperCB],[tempcheckupperCB[i,:]]))
        if printflag:
            plot(cat([localC[:,2]],[localB[:,2]]),cat([localC[:,3]],[localB[:,3]]),'--k','LineWidth',2)
            hold('on')
            for i in arange(1,size(tocheckupperCB,1)).reshape(-1):
                plot(p[tocheckupperCB[i,1],2],p[tocheckupperCB[i,1],3],'gd','LineWidth',2)
                hold('on')
            grid('on')
            xlabel('x')
            ylabel('y')
            title('Convex hull')
            fprintf('There are %d points to check\\n',size(tocheckupperCB,1))
            pause

    
    upperhullAC=matlabarray([])
    for i in arange(1,size(indexhullAC,1)).reshape(-1):
        upperhullAC=matlabarray(cat([upperhullAC],[p[indexhullAC[i,1],:]]))
    
    upperhullCB=matlabarray([])
    for i in arange(1,size(indexhullCB,1)).reshape(-1):
        upperhullCB=matlabarray(cat([upperhullCB],[p[indexhullCB[i,1],:]]))
    
    if logical_not(isempty(upperhullAC)):
        upperhullAC,Iac=maxtournamenttree(upperhullAC,false,2,nargout=2)
    
    if logical_not(isempty(upperhullCB)):
        upperhullCB,Icb=maxtournamenttree(upperhullCB,false,2,nargout=2)
    
    #------------------------------ CONVEX HULL --------------------------------
    
    upperhull=matlabarray(cat([B],[upperhullCB],[C],[upperhullAC],[A]))
    return upperhull,A,B,C