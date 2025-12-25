%%%%% Function to Update Acceleration and Velocity %%%%%%%%

function [u_n,p_n,a_n] = Update_acceleration(u_0,p_0,a_0,deltaT,beta,gamma,u_new)


a_n = ((u_new - u_0 - deltaT*p_0)/(0.5*deltaT^2) - ((1-(2*beta))*a_0))/(2*beta);

p_n = p_0 + deltaT * ((1-gamma)*a_0 + gamma*a_n);

u_n = u_new;


end
