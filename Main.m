clc
clear all
close all

syms E

global mu
mu = 398600;


R0 = [ 1600     5310    3800];
V0 = [-7.350    0.4600  2.470];
t = 3200;

[R V] = rv_from_r0v0(R0, V0, t);
r=norm(R);
v=norm(V);
[h, e, i, omega, w, theta]= coe_from_rv(R,V,mu);

disp('A.Asgharpoor     email: A.Asgharpoor@ut.ac.ir')
disp('===================================================================================')


fprintf('\n Initial position vector (km):')
fprintf('\n r0 = (%g, %g, %g)\n', R0(1), R0(2), R0(3))
fprintf('\n Initial velocity vector (km/s):')
fprintf('\n v0 = (%g, %g, %g)', V0(1), V0(2), V0(3))
fprintf('\n\n Elapsed time = %g s\n',t)
fprintf('\n Final position vector (km):')
fprintf('\n r = (%g, %g, %g)\n', R(1), R(2), R(3))
fprintf('\n Final position:')
fprintf('\n r = %g Km ',r)
fprintf('\n')
fprintf('\n Final velocity vector (km/s):')
fprintf('\n v = (%g, %g, %g)', V(1), V(2), V(3))
fprintf('\n')
fprintf('\n Final velocity :')
fprintf('\n v = %g Km/s',v)
fprintf('\n')

a=h^2/ mu*(1-e^2);
eq = E == atan(sqrt( 1-e / 1+e ) * tand(theta/2)) /2;
E= vpasolve(eq);
M = E - e*sin(E);
n = sqrt(mu/a^3);
E_rad = M; 
dE = 99999;
eps = 1e-6; % [rad] control precision of Newton's method solution

while (abs(dE) > eps)
    dE = (E_rad - e * sin(E_rad) - M)/(1 - e * cos(E_rad));
    E_rad = E_rad -  dE;
end

p_m = a*(cos(E_rad) - e);
q_m = a*sqrt(1 - e^2)*sin(E_rad);

dMdt_rad_per_s = n;
dEdt_rad_per_s = dMdt_rad_per_s/(1 - e*cos(E_rad));
dpdt_m_per_s = -a*sin(E_rad)*dEdt_rad_per_s;
dqdt_m_per_s = a*cos(E_rad)*dEdt_rad_per_s*sqrt(1 - e^2);
ra_m = a*(1 + e);  % [m] apogee 
rp_m = a*(1 - e);  % [m] perigee
rEarth_m = 6378; % [m] Earth radius at equator
Evals = 0:0.01:360.0; % [deg] values of the eccentric anomaly around orbit 
pvals = a*(cosd(Evals)-e); % [m] orbit positions
qvals = a*sqrt(1 - e^2)*sind(Evals); % [m] orbit positions
figure('color','white','Renderer', 'painters', 'Position', [500 200 700 700])
tiledlayout(2,1)

nexttile

fill(rEarth_m.*cosd(Evals),rEarth_m.*sind(Evals),[0.75 1.00 0.75]);
hold on;

plot(pvals,qvals,'.b');
plot(p_m,q_m,'.r','MarkerSize',25);
p1_m = p_m+dpdt_m_per_s*300;
q1_m = q_m+dqdt_m_per_s*300;

plot([p_m p1_m],[q_m q1_m],'r');
theta = -atan2d(dqdt_m_per_s,dpdt_m_per_s);
plot(p1_m+[0 -0.1*rEarth_m*cosd(theta+30)],q1_m+[0  0.1*rEarth_m*sind(theta+30)],'r');
plot(p1_m+[0 -0.1*rEarth_m*cosd(theta-30)],q1_m+[0  0.1*rEarth_m*sind(theta-30)],'r');

% axes
plot([-1.25*rp_m 1.25*rp_m],[0 0],'k');
plot([0 0],[-1.25*ra_m 1.25*ra_m],'k');
text(1.25*rp_m, 0.1*rEarth_m, 'p');
plot(1.25*rp_m+[0 -0.1*rEarth_m*cosd(30)],[0 0.1*rEarth_m*sind(30)],'k');
plot(1.25*rp_m+[0 -0.1*rEarth_m*cosd(30)],[0 -0.1*rEarth_m*sind(30)],'k');
text(0.1*rEarth_m, 1.25*ra_m, 'q');
plot([0 -0.1*rEarth_m*sind(30)],1.25*ra_m+[0 -0.1*rEarth_m*cosd(30)],'k');
plot([0 +0.1*rEarth_m*sind(30)],1.25*ra_m+[0 -0.1*rEarth_m*cosd(30)],'k');

axis equal
axis tight
axis off
box on

title({'Orbit Palne'});

n_rad_per_s = sqrt(mu/a^3);  % [rad/s] mean motion
n_deg_per_s = rad2deg(n_rad_per_s); % [deg/s] mean motion
M_rad = M;
E_rad = M_rad; 
dE = 99999;
eps = 1e-6; % [rad] control precision of Newton's method solution
while (abs(dE) > eps)
    dE = (E_rad - e * sin(E_rad) - M_rad)/(1 - e * cos(E_rad));
    E_rad = E_rad -  dE;
end
p_m = a*(cos(E_rad) - e);
q_m = a*sqrt(1 - e^2)*sin(E_rad);

dMdt_rad_per_s = n_rad_per_s;
dEdt_rad_per_s = dMdt_rad_per_s/(1 - e*cos(E_rad));
dpdt_m_per_s = -a*sin(E_rad)*dEdt_rad_per_s;
dqdt_m_per_s = a*cos(E_rad)*dEdt_rad_per_s*sqrt(1 - e^2);
E_deg_epoch = rad2deg(E_rad);
Rz_Omega = [ ...
    [cosd(omega) sind(omega) 0]; ...
    [-sind(omega) cosd(omega) 0]; ...
    [0 0 1]];
Rx_i = [ ...
    [1 0 0]; ...
    [0 cosd(i) sind(i)]; ...
    [0 -sind(i) cosd(i)]];
Rz_omega = [ ...
    [cosd(w) sind(w) 0]; ...
    [-sind(w) cosd(w) 0]; ...
    [0 0 1]];

T= 2*pi*a^(3/2) / sqrt(mu);
Theta_deg = theta;
Rz_hour = [ [cosd(Theta_deg) sind(Theta_deg) 0]; 
    [-sind(Theta_deg) cosd(Theta_deg) 0]; 
    [0 0 1]];

r_pq = [p_m q_m 0]';


omega_deg = w;
Omega_deg = omega;
i_deg = i;


r_ECI = inv(Rz_Omega)*inv(Rx_i)*inv(Rz_omega)*r_pq;



Evals = 0:1:360.0; % [deg] values of the eccentric anomaly around orbit 
Orbit_p = a*(cosd(Evals)-e); % [m] orbit positions
Orbit_q = a*sqrt(1 - e^2)*sind(Evals); % [m] orbit positions
deltaT_s = ((Evals-E_deg_epoch) - e*sind(Evals-E_deg_epoch))/n_deg_per_s; % [s] time since epoch along orbit

Orbit_ECI = zeros(numel(deltaT_s),3);
Orbit_LLA = zeros(numel(deltaT_s),3);

for ipt = 1:size(Orbit_ECI,1)
    r_pq = [Orbit_p(ipt) Orbit_q(ipt) 0]';
    Orbit_ECI(ipt,:) = [inv(Rz_Omega)*inv(Rx_i)*inv(Rz_omega)*r_pq]'; %[Rz_Omega*Rx_i*Rz_omega*r_pq]';
    
end

nexttile

plot3(Orbit_ECI(:,1),Orbit_ECI(:,2),Orbit_ECI(:,3),'r','LineWidth',2)
hold on
[X, Y, Z] = sphere;
X = rEarth_m *X;
Y= rEarth_m*Y;
Z= rEarth_m * Z;
surf(X,Y,Z, 'edgecolor','none')
colormap(summer)
xlabel('ECI x [m]');
ylabel('ECI y [m]');
zlabel('ECI z [m]');
title('Satellite Orbit in ECI Coordinates');
grid on