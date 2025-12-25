

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Domain Decomposition PCGM- % Acoustic Wave propagation solver %
% Direct Schur Complement - Fortran Decomposed Mesh
% Sudhi Sharma P V %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%% Mesh extraction from Preprocessor of Fortran %%%%%%%%%%%%%%%

meshpath = sprintf('./../../../../../data/meshData/');

meshfile = strcat(meshpath,'meshdim.dat');
meshdim = load(meshfile);
nP = meshdim(:,1);
nB = meshdim(:,4);
npart = meshdim(:,5);

pfile = strcat(meshpath,'points.dat');
points = load(pfile)';

tfile = strcat(meshpath,'triangles.dat');
triangles = load(tfile);
triangles = triangles(:,1:4)';

bnfile = strcat(meshpath,'boundary_nodes.dat');
nodes_BG = load(bnfile);
% n_interface = unique(n_interface);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


S_global = 0;

G_global = 0;

tol = 1*10^(-6);

%%%%%%%%%%%%%%%%%%%%%%
% New-Mark Beta Scheme

c = 1; %%% Wave Velocity

T = 1;

gamma = 0.5;

beta = 0.25;

deltaT = 0.01;

count = ceil(T/deltaT);

%%% Space coordinates %%%%

uini = zeros(nP,1);
p_n = zeros(nP,1);
a_n = zeros(nP,1);


x = points(1,:)';

y = points(2,:)';

m = 2;
n = 1;

% %% Initial Conditions %%%%%%%%%%
uini = sin(m.*pi.*x).*sin(n.*pi.*y);

% uini = cos(-x-y);

% uini = 1*exp(-50*((x-0.5).^2+(y-0.5).^2));


%   uini = zeros(1,length(p(1,:)));


% p_n = -sqrt(2)*sin(-x-y);

%  p_n(:,i) = zeros(length(p(1,:)),1);

% a_n(:,i) = zeros(length(p(1,:)),1);

%%% Gaussian Damped

% uini = 1*exp(-100*((x-0.5).^2+(y-0.5).^2)); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


U_wave = zeros(nP,count+1); %%% Total Solution at all time step

U = zeros(nP,1); %%%% Solution at one timestep


u_n = uini;

U_gamma_time = zeros(nB,count+1);

U_gamma = zeros(nB,1);

time = zeros(count+1,1);

%%% Iteration count and residual for Interface problem for all timesteps
iter = zeros(nB,count+1);
nres = zeros(nB,count+1);
%%%%%


NBcount =1;

U_wave(:,NBcount) = u_n;

while time<=T
    
time(NBcount+1)=time(NBcount)+deltaT;

NBcount = NBcount+1;

S_global = 0;

G_global = 0;


for i=1:npart


%%% Subdomain Level nodes    
n_sub = sprintf(strcat(meshpath,'nodes000%d.dat'),i);
nodes_Sub = load(n_sub);  
%%% Number of subdomain nodes
n_Sub = length(nodes_Sub);

 
%%% Boundary nodes of subdomain
nb_sub = sprintf(strcat(meshpath,'nbnodes000%d.dat'),i);
nodes_BSub = load(nb_sub);       
   
n_BSub = length(nodes_BSub);

n_ISub = n_Sub - n_BSub;

nodes_ISub = nodes_Sub(1:n_ISub,1);
    
% n_interior = nodes_ISub;

% Retrieving saved sub-domain matrices from FEniCS Assembly
path = sprintf('./data/Amats/subdom%d',i);

cd(path)

Mpath = sprintf('Mii%d',i);
Cpath = sprintf('Cii%d',i);
Ktpath = sprintf('K_Tii%d',i);
bpath = sprintf('bi%d',i);

fileA = matfile(Ktpath);
fileM = matfile(Mpath);
fileC = matfile(Cpath);
fileb = matfile(bpath);

K_Tii = fileA.K_Tii;
Mii = fileM.Mii;
Cii = fileC.Cii;
bi = fileb.bi;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Mpath = sprintf('Mig%d',i);
Cpath = sprintf('Cig%d',i);
Ktpath = sprintf('K_Tig%d',i);


fileA = matfile(Ktpath);
fileM = matfile(Mpath);
fileC = matfile(Cpath);

K_Tig = fileA.K_Tig;
Mig = fileM.Mig;
Cig = fileC.Cig;
Cgi = Cig';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Mpath = sprintf('Mgg%d',i);
Cpath = sprintf('Cgg%d',i);
Ktpath = sprintf('K_Tgg%d',i);
bpath = sprintf('bg%d',i);

fileA = matfile(Ktpath);
fileM = matfile(Mpath);
fileC = matfile(Cpath);
fileb = matfile(bpath);


K_Tgg = fileA.K_Tgg;
Mgg = fileM.Mgg;
Cgg = fileC.Cgg;
bg = fileb.bg;

K_Tgi = K_Tig';

Rpath = sprintf('R_%d',i);
fileR = matfile(Rpath);
R = double(fileR.R);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(NBcount==2)

%         un_sub = u_n(nodes_Sub);
%         pn_sub = p_n(nodes_Sub);
%         an_sub = a_n(nodes_Sub);
% 
%         un_interface = u_n(nodes_BG);
%         pn_interface = p_n(nodes_BG);
%         an_interface = a_n(nodes_BG);
% 
% 
%         un_interior = un_sub(1:n_ISub,1);
%         pn_interior = pn_sub(1:n_ISub,1);
%         an_interior = an_sub(1:n_ISub,1);
% 
%         un_sub_g = R*un_interface;
%         pn_sub_g = R*pn_interface;
%         an_sub_g = R*an_interface;
%         
%         save('un_interior','un_interior');
%         save('un_interface','un_sub_g');
%         
%         save('pn_interior','pn_interior');
%         save('pn_interface','pn_sub_g');
%         
%         save('an_interior','an_interior');
%         save('an_interface','an_sub_g');

    fileA = matfile('un_interior');
    un_interior = fileA.un_interior;
    
    fileA = matfile('un_interface');
    un_sub_g = fileA.un_sub_g;
    
    fileA = matfile('pn_interior');
    pn_interior = fileA.pn_interior;
    
    fileA = matfile('pn_interface');
    pn_sub_g = fileA.pn_sub_g;
    
    fileA = matfile('an_interior');
    an_interior = fileA.an_interior;
    
    fileA = matfile('an_interface');
    an_sub_g = fileA.an_sub_g;



        

else
    
    fileA = matfile('un_interior');
    un_interior = fileA.un_interior;
    
    fileA = matfile('un_interface');
    un_sub_g = fileA.un_sub_g;
    
    fileA = matfile('pn_interior');
    pn_interior = fileA.pn_interior;
    
    fileA = matfile('pn_interface');
    pn_sub_g = fileA.pn_sub_g;
    
    fileA = matfile('an_interior');
    an_interior = fileA.an_interior;
    
    fileA = matfile('an_interface');
    an_sub_g = fileA.an_sub_g;
       
%     un_sub_g = R*un_interface;
%     pn_sub_g = R*pn_interface;
%     an_sub_g = R*an_interface;
    
    
    

end



%%%%%% Assembling the Linear System to be solved for Implicit Method - for
%%%%%% Each Subdomain %%%%%% 

fun_Un_An_interior = ((un_interior + deltaT*pn_interior)/(deltaT^2 * beta)) ...
                       + ((1-2*beta)*an_interior)/(2*beta);

fun_Un_An_interface = ((un_sub_g + deltaT*pn_sub_g)/(deltaT^2 * beta)) + ((1-2*beta)*an_sub_g)/(2*beta);

Cmult_interior = (deltaT*gamma*fun_Un_An_interior) - pn_interior - deltaT*(1-gamma)*an_interior ;
 

Cmult_interface = (deltaT*gamma*fun_Un_An_interface) - pn_sub_g - deltaT*(1-gamma)*an_sub_g ;
 

M_n_interior = Mii * fun_Un_An_interior + Mig*fun_Un_An_interface;


M_n_interface = Mig' * fun_Un_An_interior + Mgg*fun_Un_An_interface;


C_n_interior = Cii * Cmult_interior + Cig* Cmult_interface;

C_n_interface = Cig' * Cmult_interior + Cgg* Cmult_interface;


b_t_interior = bi + M_n_interior + C_n_interior ;
b_t_interface = bg + M_n_interface + C_n_interface;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% K_transient u_n+1 = b_n


%%% u_new = K_transient\b_n will be solved by Two level DDM solver 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fprintf('Condition number of original matrix')
% cond(M_n)

S_local = K_Tgg - (K_Tgi * ((K_Tii)\ K_Tig));

G_local = b_t_interface - (K_Tgi * (K_Tii \ b_t_interior));

S_global = S_global + R' * S_local * R;

G_global = G_global + R' * G_local;



Biname = sprintf('Bt_interior_%d.mat',i);
Bgname = sprintf('Bt_interface_%d.mat',i);
save(Biname,'b_t_interior');
save(Bgname,'b_t_interface');


cd ../../../

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%% Starting Interface Solution %%%%%%%%%%%%%%%%%%%

% - S_global*U_gamma_time(:,NBcount-1);
% U_gamma_o = U_gamma_time(:,NBcount-1);
% U_gamma_o = zeros(length(node_interfacewhole),1);
% r_gamma_o = G_global;

U_gamma = S_global\G_global;

% [U_gamma,iternum,res] = DoPCGMLumped(npart,U_gamma_o,r_gamma_o,node_gamma,node_interior,R,D,tol,count,NBcount);

%  [U_gamma,iter,nres] = DoPCGM_NeumannTwoLevel(npart,node_c,node_r,...
%     node_corner,node_gamma,node_interior,G_global,U_gamma_o,r_gamma_o,Rc,Rr,Bc,R,D,tol,tolc);

U_gamma_t = U_gamma;
U_gamma_time(:,NBcount) = U_gamma;

% iter(:,NBcount) = iternum;
% 
% nres(:,NBcount) = res;

% % node_global = nonzeros(unique(GLnode));
% % 
% % for nr = 1:size(n_interface,1)
% %  U(node_global(n_interface(nr),:)) = U_gamma(nr);
% % end


%%- Scalar Assembly 
 nBG = length(U_gamma);
 for j = 1:nBG
   tempGB = nodes_BG(j);
   U(tempGB) = U_gamma(j);
 end


% Constructing the A and b Matrices after getting interface solution. Can
% be avoided by saving them and reloading at this stage.

for j=1:npart
  
    
%%% Subdomain Level nodes    
n_sub = sprintf(strcat(meshpath,'nodes000%d.dat'),j);
nodes_Sub = load(n_sub);   

n_Sub = length(nodes_Sub);

path = sprintf('./data/Amats/subdom%d',j);

cd(path)

% Retrieving saved sub-domain matrices
Ktpath = sprintf('K_Tii%d',j);
fileA = matfile(Ktpath);
K_Tii = fileA.K_Tii;

Ktpath = sprintf('K_Tig%d',j);
fileA = matfile(Ktpath);
K_Tig = fileA.K_Tig;

Biname = sprintf('Bt_interior_%d.mat',j);
Bgname = sprintf('Bt_interface_%d.mat',j);
filebi = matfile(Biname);
filebg = matfile(Bgname);

b_t_interior = filebi.b_t_interior;
b_t_interface = filebg.b_t_interface;

Rpath = sprintf('R_%d',j);
fileR = matfile(Rpath);
R = double(fileR.R);


U_interior  = K_Tii\(b_t_interior - (K_Tig * R * U_gamma));

U_interior_t = U_interior;

npi = length(U_interior);


%%- Scalar Assembly 
    for nn = 1:npi
       U(nodes_Sub(nn)) = U_interior(nn);
    end


% % % for ni = 1:size(n_interior,1)
% % %  U(GLnode(n_interior(ni),:,j)) = U_interior(ni);
% % % end

save('unew_interior','U_interior_t');

cd ../../../

end


%%%%%%%%%%%%%%%%% DDM solver ends %%%%%%%%%%%%%%%%%%%%%%%%%%

u_new = U;

U_wave(:,NBcount) = u_new;

% [u_n,p_n,a_n] = Update_acceleration(u_n,p_n,a_n,deltaT,beta,gamma,u_new);

Update_acceleration_Subdomain(npart,deltaT,beta,gamma,U_gamma_t,nodes_BG);

end


fprintf("Condition number of Schur Complement matrix")
cond(S_global)
%%%%%%%%%%%%%%%%%%%%%%%
% % % Analytic solution

nt = length(time);

SOL = zeros(length(x),nt); 
S = uini';

for i=1:nt
    SOL(:,i) = S.*cos(c*pi*(sqrt(m^2+n^2)*time(i)));
end

% %Absolute error
% E = abs(SOL-U_wave);
% max_E = max(max(E))


for j =1:nt
    
    
     s1 = pdesurf(points,triangles,U_wave(:,j));
     axis([0 1 0 1 -1 1]);
     zlim([-1 1])
     hold on
     colorbar
     pause(0.05);
     
         
    if(j<count)
         delete(s1)
    end

end



%%% Draw the time vaiation at a particular point

pos = 402;
figure(20)
plot(time,U_wave(pos,:));


figure(10)
plot(time,U_wave(pos,:),time,SOL(pos,:));

E_norm = zeros(nt,1);
for kk = 1:nt
E_norm(kk) = norm(U_wave(:,kk)-SOL(:,kk));
end
figure(12)
plot(time,E_norm);

E_norm(kk)




% 
% clear all;
% if exist('./data/Amats')
% [status, message, messageid] = rmdir('./data/Amats', 's')
% end
% delete *.mat  % Command deletes all files in folder with .mat extension








