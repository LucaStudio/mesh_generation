%%
%        Project: Fluid - structure interaction on deformable surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: March 28th, 2014
%    Last update: April 1st, 2014
%
%    Description: 
%          Input: 
%         Output:

clear all
close all
clc

D = 3;
dim1min = -5; 
dim1max = 5;
Ndim1 = 100;
dim2min = 0;
dim2max = 5;
Ndim2 = 100;
dim3min = 0;
dim3max = 10;
Ndim3 = 100;
Nlinesdim1 = 1;
Nlinesdim2 = 1;
Nlinesdim3 = 1;
edgeflag = 0;
faceflag = 0;
cellflag = 0;
printflag = 0;

mesh = generate_regularcuboids_mesh(D,dim1min,dim1max,Ndim1,dim2min,dim2max,Ndim2,dim3min,dim3max,Ndim3,Nlinesdim1,Nlinesdim2,Nlinesdim3,edgeflag,faceflag,cellflag,printflag);

changeddomain = changedomain(D,mesh,8,10,0,2*pi,0,10);

funcs = {'r*cos(theta)','r*sin(theta)','z'};
args = {'r','theta','z'};

analyticalmorph = analytictransformation(D,changeddomain,funcs,args);

Dshow = 3;
pointcolor = 'r';
pointdim = 2;
linecolor = 'b';
linedim = 1;
titlestring = 'Regular mesh';
xstring = 'x';
ystring = 'y';
zstring = 'z';
shownodes = true;
showlines = false;
shownodelabels = false;
showedgelabels = false;
showfacelabels = false;
showcelllabels = false;
nodelabelcolor = [0 0 0]; % black in rgb
edgelabelcolor = [1 0 0]; % red in rgb
facelabelcolor = [0 1 0]; % green in rgb
celllabelcolor = [0 0 1]; % blue in rgb
labelsize = 12;
xfigsize = 1000;
yfigsize = 1000;
xaxismin = -15;
xaxismax = 15;
yaxismin = -15;
yaxismax = 15;
zaxismin = -15;
zaxismax = 15;

f1 = show_mesh(Dshow,mesh,shownodes,showlines,shownodelabels,showedgelabels,showfacelabels,showcelllabels,pointcolor,pointdim,linecolor,linedim,nodelabelcolor,edgelabelcolor,facelabelcolor,celllabelcolor,labelsize,titlestring,xstring,ystring,zstring,xfigsize,yfigsize,xaxismin,xaxismax,yaxismin,yaxismax,zaxismin,zaxismax);

f3 = show_mesh(Dshow,changeddomain,shownodes,showlines,shownodelabels,showedgelabels,showfacelabels,showcelllabels,pointcolor,pointdim,linecolor,linedim,nodelabelcolor,edgelabelcolor,facelabelcolor,celllabelcolor,labelsize,titlestring,xstring,ystring,zstring,xfigsize,yfigsize,xaxismin,xaxismax,yaxismin,yaxismax,zaxismin,zaxismax);

f5 = show_mesh(Dshow,analyticalmorph,shownodes,showlines,shownodelabels,showedgelabels,showfacelabels,showcelllabels,pointcolor,pointdim,linecolor,linedim,nodelabelcolor,edgelabelcolor,facelabelcolor,celllabelcolor,labelsize,titlestring,xstring,ystring,zstring,xfigsize,yfigsize,xaxismin,xaxismax,yaxismin,yaxismax,zaxismin,zaxismax);


D = 2;
dim1min = -5; 
dim1max = 5;
Ndim1 = 100;
dim2min = 0;
dim2max = 5;
Ndim2 = 100;
dim3min = 0;
dim3max = 10;
Ndim3 = 100;
Nlinesdim1 = 1;
Nlinesdim2 = 1;
Nlinesdim3 = 1;
edgeflag = 0;
faceflag = 0;
cellflag = 0;
printflag = 0;

mesh = generate_regularcuboids_mesh(D,dim1min,dim1max,Ndim1,dim2min,dim2max,Ndim2,dim3min,dim3max,Ndim3,Nlinesdim1,Nlinesdim2,Nlinesdim3,edgeflag,faceflag,cellflag,printflag);

changeddomain = changedomain(D,mesh,8,10,0,2*pi,0,10);

funcs = {'r*cos(theta)','r*sin(theta)'};
args = {'r','theta'};

analyticalmorph = analytictransformation(D,changeddomain,funcs,args);

Dshow = 2;
pointcolor = 'r';
pointdim = 2;
linecolor = 'b';
linedim = 1;
titlestring = 'Regular mesh';
xstring = 'x';
ystring = 'y';
zstring = 'z';
shownodes = true;
showlines = false;
shownodelabels = false;
showedgelabels = false;
showfacelabels = false;
showcelllabels = false;
nodelabelcolor = [0 0 0]; % black in rgb
edgelabelcolor = [1 0 0]; % red in rgb
facelabelcolor = [0 1 0]; % green in rgb
celllabelcolor = [0 0 1]; % blue in rgb
labelsize = 12;
xfigsize = 1000;
yfigsize = 1000;
xaxismin = dim1min - 5;
xaxismax = dim1max + 5;
yaxismin = dim2min - 5;
yaxismax = dim2max + 5;
zaxismin = dim3min - 5;
zaxismax = dim3max + 5;

f2 = show_mesh(Dshow,mesh,shownodes,showlines,shownodelabels,showedgelabels,showfacelabels,showcelllabels,pointcolor,pointdim,linecolor,linedim,nodelabelcolor,edgelabelcolor,facelabelcolor,celllabelcolor,labelsize,titlestring,xstring,ystring,zstring,xfigsize,yfigsize,xaxismin,xaxismax,yaxismin,yaxismax,zaxismin,zaxismax);

f4 = show_mesh(Dshow,changeddomain,shownodes,showlines,shownodelabels,showedgelabels,showfacelabels,showcelllabels,pointcolor,pointdim,linecolor,linedim,nodelabelcolor,edgelabelcolor,facelabelcolor,celllabelcolor,labelsize,titlestring,xstring,ystring,zstring,xfigsize,yfigsize,xaxismin,xaxismax,yaxismin,yaxismax,zaxismin,zaxismax);

f6 = show_mesh(Dshow,analyticalmorph,shownodes,showlines,shownodelabels,showedgelabels,showfacelabels,showcelllabels,pointcolor,pointdim,linecolor,linedim,nodelabelcolor,edgelabelcolor,facelabelcolor,celllabelcolor,labelsize,titlestring,xstring,ystring,zstring,xfigsize,yfigsize,xaxismin,xaxismax,yaxismin,yaxismax,zaxismin,zaxismax);

