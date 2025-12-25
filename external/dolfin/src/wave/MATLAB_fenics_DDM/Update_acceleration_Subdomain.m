%%%%% Function to Update Acceleration and Velocity %%%%%%%%

function Update_acceleration_Subdomain(npart,deltaT,beta,gamma,U_gamma_t,nodes_BG)

meshpath = sprintf('./../../../../../data/meshData/');

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


path = sprintf('./data/Amats/subdom%d',i);

cd(path)
    

Rpath = sprintf('R_%d',i);
fileR = matfile(Rpath);
R = double(fileR.R);



fileA = matfile('un_interior.mat');
un_interior = fileA.un_interior;

fileA = matfile('un_interface.mat');
un_sub_g = fileA.un_sub_g;


fileA = matfile('pn_interior.mat');
pn_interior = fileA.pn_interior;

fileA = matfile('pn_interface.mat');
pn_sub_g = fileA.pn_sub_g;


fileA = matfile('an_interior.mat');
an_interior = fileA.an_interior;

fileA = matfile('an_interface.mat');
an_sub_g = fileA.an_sub_g;

% u0_sub = u_0(nodes_Sub);
% p0_sub = p_0(nodes_Sub);
% a0_sub = a_0(nodes_Sub);
% 
% u0_interface = u_0(nodes_BG);
% p0_interface = p_0(nodes_BG);
% a0_interface = a_0(nodes_BG);
% 
% 
% u0_interior = u0_sub(1:n_ISub,1);
% p0_interior = p0_sub(1:n_ISub,1);
% a0_interior = a0_sub(1:n_ISub,1);
% 
% u0_sub_g = R*u0_interface;
% p0_sub_g = R*p0_interface;
% a0_sub_g = R*a0_interface;


unew_sub_g = R*U_gamma_t;


u_i = matfile('unew_interior.mat');

unew_interior = u_i.U_interior_t;


anew_interface = ((unew_sub_g - un_sub_g - deltaT*pn_sub_g)/(0.5*deltaT^2) - ((1-(2*beta))*an_sub_g))/(2*beta);

pnew_interface = pn_sub_g + deltaT * ((1-gamma)*an_sub_g + gamma*anew_interface);



anew_interior = ((unew_interior - un_interior - deltaT*pn_interior)/(0.5*deltaT^2) - ((1-(2*beta))*an_interior))/(2*beta);

pnew_interior = pn_interior + deltaT * ((1-gamma)*an_interior + gamma*anew_interior);


un_interior = unew_interior;
un_sub_g = unew_sub_g;

pn_interior = pnew_interior;
pn_sub_g = pnew_interface;

an_interior = anew_interior;
an_sub_g = anew_interface;



%%%%% Should keep the same name as the one used for saving old variables

save('un_interior','un_interior');
save('un_interface','un_sub_g');

save('an_interior','an_interior');
save('an_interface','an_sub_g');

save('pn_interior','pn_interior');
save('pn_interface','pn_sub_g');


cd ../../../

end



end
