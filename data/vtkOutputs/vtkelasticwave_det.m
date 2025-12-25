% test_vtk_write_binary
clearvars

% x = vtx_coord(:,1);
% y = vtx_coord(:,2);
% z = vtx_coord(:,3);

nbparam = load('./NBparam.dat');
T = nbparam(1);
dt = nbparam(2);
nbcount = T/dt;
points = dlmread('./points.dat','');   %dlmread('points.dat','');
x = points(:,1);
y = points(:,2);
z = points(:,3);
np = length(x);

tetra = dlmread('./tetrahedrons4c.dat','');   %dlmread('points.dat','');


for j=1:nbcount

    
%%%% Mean %%%%%%%%%%%%%%%%%%%%%
    if(j < 10)
        str = sprintf('U_mean_000%d.dat',j);
    elseif (j < 100) 
        str = sprintf('U_mean_00%d.dat',j);
    else
        str = sprintf('U_mean_0%d.dat',j);
    end 
    
solution = dlmread(str,'');
U = solution(1:np);
V = solution(np+1:2*np);
W = solution(2*np+1:3*np);

% setup filename
str2 = sprintf('U_%d.vtk',j);
filename = str2;  %% Binary (use paraview)
data_title = 'displacement';

% organize data
% .type field must be 'scalar' or 'vector'
% number of vector components must equal dimension of the mesh
% data_struct.type = 'scalar';
% data_struct.name = 'displacement';
% data_struct.data = solution;

data_struct.type = 'vector';
data_struct.name = 'displacements';
data_struct.data = [U,V,W];   %% Need to provide nodal components

% if *all data* is oriented column-wise, i.e. lots of rows, few columns,
% then set this to *false*
flipped = false;
% otherwise, if you want to transpose the data, then set this to *true*

% write the file
tic
stat = vtk_writeFn(filename,data_title,points,tetra,data_struct,flipped);
toc



dispout = ['VTK output: ', filename];
disp(dispout)


%%%%%%%%%%%%%%%%%% SD %%%%%%%%%%%%%%%%


    if(j < 10)
        str = sprintf('sd_000%d.dat',j);
    elseif (j < 100) 
        str = sprintf('sd_00%d.dat',j);
    else
        str = sprintf('sd_0%d.dat',j);
    end 
    
solution = dlmread(str,'');
U = solution(1:np);
V = solution(np+1:2*np);
W = solution(2*np+1:3*np);

% setup filename
str2 = sprintf('sd_%d.vtk',j);
filename = str2;  %% Binary (use paraview)
data_title = 'sd';

% organize data
% .type field must be 'scalar' or 'vector'
% number of vector components must equal dimension of the mesh
% data_struct.type = 'scalar';
% data_struct.name = 'displacement';
% data_struct.data = solution;

data_struct.type = 'vector';
data_struct.name = 'sd';
data_struct.data = [U,V,W];   %% Need to provide nodal components

% if *all data* is oriented column-wise, i.e. lots of rows, few columns,
% then set this to *false*
flipped = false;
% otherwise, if you want to transpose the data, then set this to *true*

% write the file
tic
stat = vtk_writeFn(filename,data_title,points,tetra,data_struct,flipped);
toc



dispout = ['VTK output: ', filename];
disp(dispout)


%%%%%%%%%%% PCE Coefficients %%%%%%%%%%%%%%


for pc = 4:4:20
    
        
    pcname = sprintf('pce');
    
        if(pc < 10)
            str1 = sprintf('pce_0%d',pc);
        elseif(pc < 100) 
            str1 = sprintf('pce_%d',pc);
        end 
        


        if(j < 10)
            str = sprintf('%s_000%d.dat',str1,j);
        elseif (j < 100) 
            str = sprintf('%s_00%d.dat',str1,j);
        else
            str = sprintf('%s_0%d.dat',str1,j);
        end 

    solution = dlmread(str,'');
    U = solution(1:np);
    V = solution(np+1:2*np);
    W = solution(2*np+1:3*np);

    % setup filename
    str2 = sprintf('pce_%d_%d.vtk',pc,j);
    filename = str2;  %% Binary (use paraview)
    data_title = 'displacement';

    % organize data
    % .type field must be 'scalar' or 'vector'
    % number of vector components must equal dimension of the mesh
    % data_struct.type = 'scalar';
    % data_struct.name = 'displacement';
    % data_struct.data = solution;

    data_struct.type = 'vector';
    data_struct.name = 'displacements';
    data_struct.data = [U,V,W];   %% Need to provide nodal components

    % if *all data* is oriented column-wise, i.e. lots of rows, few columns,
    % then set this to *false*
    flipped = false;
    % otherwise, if you want to transpose the data, then set this to *true*

    % write the file
    tic
    stat = vtk_writeFn(filename,data_title,points,tetra,data_struct,flipped);
    toc



    dispout = ['VTK output: ', filename];
    disp(dispout)


end


end

