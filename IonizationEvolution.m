clear clc

% This script contains functions for computing multilevel ionization
% evolution of high atomic number gas due to interation with high intenstity
% laser.  
% It includes algorithms for solving the complete multilevel ODE system 
% as well as an analytical solution using reduced 4-level ODEs.

% The functions below use krypton with 4J 2ps laser as an exmple. 
% Solutions for other gas mediums or laser parameters can also be obtianed
% by simply modifying the "ionized medium" in ionAtLocation() or onAtLocation_4level() 
% and the laser parameters in gaussian()

% Example: Kr at focal location of CO2 laser 
% full ODE solution
[Kr_fullODE,ne_fullODE,taxis] = ionAtLocation(0.1e-10,0e-6,1,0);

% reduced 4 level ODE solution
[Kr_reduce4ODE,ne_reduce4ODE,taxis] = ionAtLocation_4level(0.1e-10,0e-6,0,0);


function [Ion_all, ne_all, taxis] = ionAtLocation(zf,r, picture, peakIntensity)
% obtain ionization states of krypton for prescribed locations.
% Inputs:
% % zf: (m) longitudinal distance from focus position
% % r: (m) radial distance from center (m)
% % picture: = 1 for plotting; ~= 1 to turn off plotting
% % peakIntensity: (W/m^2) = 0 for default CO2 laser
% Outputs:
% % Ion_all: (%) Evolution of number densities for all ion 
% % ne_all: (per one neutral atom) Evolution of total number of electrons
% % taxis: (s) computational time steps

% constants
ep0 = 8.8542e-12; % (F/m) vacuume permitivoty
mu0 = pi*4e-7; % (H/m) vacuume permeability
c = 299792458.0;        % (m/s)


% Obtain electric field magnitude for entire laser duration at (zf,r)
if zf < eps     % adjust zf to avoid singularity
    zf = eps;
end
rsq = (r)^2; % (m^2) radial distance from center ^2

[taxis, Ex] = gaussian(zf, rsq, peakIntensity);
dt = taxis(2)-taxis(1);     % (s) timestep of computation

% Ionizated medium

% % Ar ionization energy
% % Znum = 18;
% % energyIon = [15.75962, 27.62967, 40.74, 59.58, 74.84, 91.29, 124.41, 143.457, 422.60, 479.76, 540.4, 619.0, 685.5, 755.13, 855.5, 918.375, 4120.666, 4426.224]; %(eV)
% % l_all = [0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0];
% 
% % Xe ionization energy
% % Znum = 56;
% % energyIon = [12.1298, 20.975, 31.05]; %(eV) higher levels omitted
% % l_all = [0, 1, 1];
% 
% Kr ionization energy
Znum = 36;
energyIon = [13.9996055, 24.35984, 35.838, 50.85, 64.69, 78.49, 109.13, 125.082, 233, 268, 308, 350, 391, 446, 492, 540, 591, 640]; % (eV) higher levels omitted
l_all = [0, 1, 1, 0, 1, 1, 0, 0, 0, 2, 3, 3, 2, 0, 2, 3, 3, 2];


% Initialization

% % compute ion parameters of each level

UH = 13.6;      % (eV)  ionization potential of Hydrogen
wa = 4.13e16;   % (s^-1) atomic unit for frequency
Ea = 510;       % (GV/m) atomic unit for electric field
e = exp(1.);

nstar = zeros(length(energyIon),1);
coeflm = zeros(length(energyIon),1);
coefnz = zeros(length(energyIon),1);
for Z = 1:length(energyIon)
    energy = energyIon(Z);              % ionization potential    
    l = l_all(Z);                       % orbital quantum number  
    m=0;                                % projection of l
    nstar(Z)  = Z * sqrt(UH/energy);    % effective principle quantum number
    coeflm(Z) = (2*l+1)*factorial(l+m)/( 2*pi*nstar(Z)*2^m*factorial(m)*factorial(l-m));
    coefnz(Z) = wa*sqrt(3/pi) * (Z/nstar(Z))^(-2.5+6*nstar(Z)-3*m) * (e/nstar(Z))^(2*nstar(Z)) * 2^(4*nstar(Z)-m-2);
end

% % arrays for single timestep
Ion = zeros(length(energyIon)+1,1);
Ion(1) = 1;     % initialize with neutral atom

newIon = zeros(length(energyIon)+1,1);
W = zeros(length(energyIon),1);

% % collection of all time sequence
W_all = [];
Ion_all = [];
ne_all = [];


% Solve time evolution of ionization
for n=1:length(taxis)
    El = Ex(n);
        
    if El-1e-15>0       % skip 0 electric field to avoid singularity      
        newIon(:) = 0;
        
        % compute ionization rate W for this electric field magnitude
        for Z = 1:length(energyIon)
            
            W(Z) = coeflm(Z)*coefnz(Z) * exp(-2*510/(3*El) * (Z/nstar(Z))^3) * (Ea/El)^(2*nstar(Z)-m-1.5);            
            W(Z) = W(Z)/ (3/pi * El/Ea * (UH/energy)^1.5)^0.5;      % if use W_{DC} instead of W_{AC}

        end        
        
        W_all = [W_all,W];
        
        %                         %%%% ionization solution
        %             % --- based on solution of full differential equation systems ---

        Wmatrix = ConstructMatrix(W);
        [Wvec, Wi] = eig(Wmatrix);  % Wvec : eigenvector V from the W matrix
        C = Wvec\Ion(:);            % solve for coefficient C from VC = N(t) 
        for Z = 1:(length(energyIon)+1)
            for i = 1:length(Wi)
                newIon(Z) = newIon(Z) + C(i)*Wvec(Z,i)*exp(Wi(i,i)*dt);
            end
        end
        
        % update ion density array
        Ion(:) = newIon;
        for ikr = 1:length(newIon)
            if Ion(ikr)<0
                Ion(ikr)=0;
            end
        end       
    end
        % compute total number of electrons at current ionization state
        ne = dot(Ion,(0:length(energyIon)));
        
        ne_all = [ne_all, ne]; 
        Ion_all = [Ion_all,Ion];

end

Ion_final = Ion;
ne_final = ne;

%% plot result

if picture == 1

    Ion = Ion_all;
      
    % number of electrons
    figure
    plot(taxis,ne_all,'LineWidth',2)
    title(['E0 = ', num2str(max(Ex))])
    xlabel('time (s)')
    ylabel('number of electrons')
    set(gca,'FontSize',20);

    figure
    %taxis = taxis(end:-1:1); % reverse t axis to show evolution of ionizationt
    taxis = taxis*1e12;
    plot(taxis,Ion(1,:), taxis,Ion(2,:), taxis,Ion(3,:), taxis,Ion(4,:), taxis,Ion(5,:),taxis,Ion(6,:), taxis,Ion(7,:),taxis,Ion(8,:), taxis,Ion(9,:), taxis,Ion(10,:), taxis,Ion(11,:), taxis,Ion(12,:), taxis,Ion(13,:), taxis,Ion(14,:), taxis,Ion(15,:), taxis,Ion(16,:) ,taxis,Ion(17,:), 'LineWidth',2)
    hold on
    plot(taxis, Ex./max(Ex),'LineStyle',':','Color','black')
    legend('n0','n1','n2','n3','n4','n5','n6','n7', 'n8', 'n9', 'n10', 'n11', 'n12', 'n13', 'n14', 'n15', 'n16', 'normalized |E|');    
    xlabel('time (ps)')
    ylabel('relative population of ionization levels')
    title(['r = ',num2str(r*1e6)])
    set(gca,'FontSize',20)
    
end

end
function [Ion_all, ne_all, taxis] = ionAtLocation_4level(zf,r, picture, peakIntensity)
% obtain ionization states for prescribed locations
% using the reduced 4 level ODE system
% Inputs:
% % zf: (m) longitudinal distance from focus position
% % r: (m) radial distance from center (m)
% % picture: = 1 for plotting; ~= 1 to turn off plotting
% % peakIntensity: (W/m^2) = 0 for default CO2 laser
% Outputs:
% % Ion_all: (%) Evolution of number densities for all ion 
% % ne_all: Evolution of total number of electrons
% % taxis: (s) computational time steps

% constants
ep0 = 8.8542e-12; % (F/m) vacuume permitivoty
mu0 = pi*4e-7; % (H/m) vacuume permeability
c = 299792458.0;        % (m/s)

% Obtain electric field magnitude for entire laser duration at (zf,r)
if zf < eps     % adjust zf to avoid singularity
    zf = eps;
end
rsq = (r)^2; % (m^2) radial distance from center ^2

[taxis, Ex] = gaussian(zf, rsq, peakIntensity);
dt = taxis(2)-taxis(1);

% Ar ionization energy
% Znum = 18;
% energyIon = [15.75962, 27.62967, 40.74, 59.58, 74.84, 91.29, 124.41, 143.457, 422.60, 479.76, 540.4, 619.0, 685.5, 755.13, 855.5, 918.375, 4120.666, 4426.224]; %(eV)
% l_all = [0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0];

% Xe ionization energy
% Znum = 56;
% energyIon = [12.1298, 20.975, 31.05]; %(eV) higher levels omitted
% l_all = [0, 1, 1];

% energyIon = [12.1298, 21.20979, 32.1230]; %(eV) for XenonZ=56 

% Kr ionization energy
Znum = 36;
energyIon = [13.9996055, 24.35984, 35.838, 50.85, 64.69, 78.49, 109.13, 125.082, 233, 268, 308, 350, 391, 446, 492, 540, 591, 640]; % (eV)
l_all = [0, 1, 1, 0, 1, 1, 0, 0, 0, 2, 3, 3, 2, 0, 2, 3, 3, 2];


% Initialization
% % arrays for single timestep
Ion = zeros(length(energyIon)+1,1);
Ion(1) = 1;     % initialize with neutral atom
p = 0;          % lowest ion level of the Reduced ODE, neutral 0
LowDel = 1e-5;  % lower threashold approximating fully ionized number density
newIon = zeros(4,1);
W = zeros(length(energyIon),1);

% % collection of all time sequence
W_all = [];
Ion_all = [];
ne_all = [];


% Compute ion parameters of each level

UH = 13.6;      % (eV)  ionization potential of Hydrogen
wa = 4.13e16;   % (s^-1) atomic unit for frequency
Ea = 510;       % (GV/m) atomic unit for electric field
e = exp(1.);

nstar = zeros(length(energyIon),1);
coeflm = zeros(length(energyIon),1);
coefnz = zeros(length(energyIon),1);
for Z = 1:length(energyIon)
    energy = energyIon(Z);          % ionization potential    
    l = l_all(Z);                   % orbital quantum number  
    m=0;                            % projection of l
    nstar(Z)  = Z * sqrt(UH/energy);    % effective principle quantum number
    coeflm(Z) = (2*l+1)*factorial(l+m)/( 2*pi*nstar(Z)*2^m*factorial(m)*factorial(l-m));
    coefnz(Z) = wa*sqrt(3/pi) * (Z/nstar(Z))^(-2.5+6*nstar(Z)-3*m) * (e/nstar(Z))^(2*nstar(Z)) * 2^(4*nstar(Z)-m-2);
end

% Solve time evolution of ionization
for n=1:length(taxis)
    El = Ex(n);
        
    if El-1e-15>0       % skip 0 electric field to avoid singularity      
        newIon(:) = 0;
        
        % compute ionization rate W for this electric field magnitude
        for Z = 1:length(energyIon)
           
            W(Z) = coeflm(Z)*coefnz(Z) * exp(-2*510/(3*El) * (Z/nstar(Z))^3) * (Ea/El)^(2*nstar(Z)-m-1.5);            
            W(Z) = W(Z)/ (3/pi * El/Ea * (UH/energy)^1.5)^0.5;      % if use W_{DC} instead of W_{AC}

        end        
        
        W_all = [W_all,W];
        
                     %%%% ionization solution
        % --- based on analytical solution of 4by4 ODE system ---
                
        newIon(1) = Ion(1+p) * exp(-W(1+p)*dt);
        % n^+, n^2+ , n^3+
        W12 = W(1+p)-W(2+p);
        W13 = -W(1+p)+W(3+p);
        W23 = W(2+p)-W(3+p);
        if abs(W12)>eps  % if ionization can happen for the lowest level
            newIon(2) = W(1+p)/W12 * (exp(-W(2+p)*dt)-exp(-W(1+p)*dt)) * Ion(1+p) + exp(-W(2+p)*dt)*Ion(2+p);
            newIon(3) = (exp(-W(3+p)*dt) * (W(3+p)^2*Ion(3+p)+W(1+p)*W(2+p)*Ion(1+p)+W(1+p)*W(2+p)*Ion(2+p)-W(2+p)*W(3+p)*Ion(2+p) + W(1+p)*W(2+p)*Ion(3+p)-W(1+p)*W(3+p)*Ion(3+p)-W(2+p)*W(3+p)*Ion(3+p)))/((W(1+p)-W(3+p))*(W(2+p)-W(3+p)))-(W(2+p)*exp(-W(2+p)*dt)*(W(1+p)*Ion(1+p)+W(1+p)*Ion(2+p)-W(2+p)*Ion(2+p)))/((W(1+p)- W(2+p))*(W(2+p) - W(3+p))) + (W(1+p)*W(2+p)*Ion(1+p)*exp(-W(1+p)*dt))/((W(1+p) - W(2+p))*(W(1+p) - W(3+p)));
            newIon(4) =  Ion(1+p) + Ion(2+p) + Ion(3+p) + Ion(4+p) - (exp(-W(3+p)*dt)*(W(3+p)^2*Ion(3+p) + W(1+p)*W(2+p)*Ion(1+p) + W(1+p)*W(2+p)*Ion(2+p) - W(2+p)*W(3+p)*Ion(2+p) + W(1+p)*W(2+p)*Ion(3+p) - W(1+p)*W(3+p)*Ion(3+p) - W(2+p)*W(3+p)*Ion(3+p)))/((W(1+p) - W(3+p))*(W(2+p) - W(3+p))) + (W(3+p)*exp(-W(2+p)*dt)*(W(1+p)*Ion(1+p) + W(1+p)*Ion(2+p) - W(2+p)*Ion(2+p)))/((W(1+p) - W(2+p))*(W(2+p) - W(3+p))) - (W(2+p)*W(3+p)*Ion(1+p)*exp(-W(1+p)*dt))/((W(1+p) - W(2+p))*(W(1+p) - W(3+p)));
        else
            newIon(2) = Ion(2+p);
            newIon(3) = Ion(3+p);
            newIon(4) = Ion(4+p);
         end
        Ion(1+p:4+p) = newIon;  % replace Ion list with the new population
        
        % shift 4 levels up when the current lowest level is considered to
        % be fully ionized.
        if ((newIon(1)<=LowDel) && (p+3<length(energyIon)))
            newIon(1) = newIon(2);
            newIon(2) = newIon(3);
            newIon(3) = newIon(4);
            newIon(4) = 0;
            p = p+1;
        end

    end
        % compute total number of electrons at current ionization state
        ne = dot(Ion,(0:length(energyIon)));
        
        ne_all = [ne_all, ne]; 
        Ion_all = [Ion_all,Ion];

end

Ion_final = Ion;
ne_final = ne_all(end);

%% plot result

if picture == 1

    Ion = Ion_all;
      
    % number of electrons
    figure
    plot(taxis,ne_all,'LineWidth',2)
    title(['E0 = ', num2str(max(Ex))])
    xlabel('time (s)')
    ylabel('number of electrons')
    set(gca,'FontSize',20);

    figure
    %taxis = taxis(end:-1:1); % reverse t axis to show evolution of ionizationt
    taxis = taxis*1e12;
    plot(taxis,Ion(1,:), taxis,Ion(2,:), taxis,Ion(3,:), taxis,Ion(4,:), taxis,Ion(5,:),taxis,Ion(6,:), taxis,Ion(7,:),taxis,Ion(8,:), taxis,Ion(9,:), taxis,Ion(10,:), taxis,Ion(11,:), taxis,Ion(12,:), taxis,Ion(13,:), taxis,Ion(14,:), taxis,Ion(15,:), taxis,Ion(16,:) ,taxis,Ion(17,:), 'LineWidth',2)
    hold on
    plot(taxis, Ex./max(Ex),'LineStyle',':','Color','black')
    legend('n0','n1','n2','n3','n4','n5','n6','n7', 'n8', 'n9', 'n10', 'n11', 'n12', 'n13', 'n14', 'n15', 'n16', 'normalized |E|');
    xlabel('time (ps)')
    ylabel('relative population of ionization levels')
    title(['r = ',num2str(r*1e6)])
    set(gca,'FontSize',20)
    
end

end

function RateMatrix = ConstructMatrix(W)
% construct lower bidiagonal matrix of ionization rate
% Input: W (s^-1) m*1 array of ionization rate for all levels
% Output: RateMatrix: (m+1)*(m+1) matrix
n = length(W)+1;
A = zeros(n,n);
W(W<eps) = 0;
% fill diagonal elements
A(1:n+1:end) = [-W;0];
% fill subdiagonal elements
A(2:n+1:end) = W;
RateMatrix = A;
end

function [t, Ex] = gaussian(zf, rsq, peakIntensity)
% generate electric field magnitude for linearly polarized gaussian laser pulse
% over a period of time
% Inputs:
%  % zf: (m) longitudinal distance from the focus position
%  % rsq: (m^2) square of radial distance from center
%  % peakIntensity: (W/m^2) peak intensity of laser
%  % % peakIntensity = 0 for default CO2 laser 
% Outputs:
%  % t: (s) time axis for Ex  
%  % Ex: (GV/m) electric field magnitude

% constants
ep0 = 8.8542e-12;       % (F/m) vacuume permitivoty
c = 299792458;          % (m/s) light speed

% CO2 laser, default laser parameter; customizable if stated otherwise
l0 = 9.2e-6;            % (m) wavelength
beam_waist = 40e-6/2;   % (m) beam waist
tau = 2e-12;            % (s) pulse duration (FWHM)
energy = 4;             % (J) laser energy

f0 = c/l0;              % (Hz) frequency
omega = 2*pi*f0;        % (rad) frequency
tau0 = 1/f0;            % (s) period
dt = tau0/50;           % (s) time resolution (<= tau0/20)

power = energy/tau;     % (w) average laser power
E0 = 2/beam_waist * sqrt(power/(pi*c*ep0));     %(V/m) max amplitude
if peakIntensity ~= 0
    E0 = sqrt(2*peakIntensity/(c*ep0));
end
E0 = E0*1e-9;           % (GV/m)

stand_dv = tau/(2*sqrt(2*log(2)));
tau = stand_dv*6;                   % (s) total laser pulse duration
t_num = ceil(tau/dt);
t = linspace(0,dt*t_num,t_num);     % time axis for total laser pulse

M2 = 1;
K = omega/c;                    % (m^-1) wave number
zR = pi*beam_waist^2*f0/c/M2;   % (m) rayleigh distance
R = zf.*(1+(zR./zf).^2);        % (m) radius of curvature at zf
phi = atan(-zf/zR);             % gouy phase
beam_waist_evolve = beam_waist*sqrt(1+(zf/zR).^2);  % spatial envelop of beam waist
envelop = 1./exp(18*((t-0.5*tau)/tau).^2);          % temporal envelop of laser intensity
c1 = E0*beam_waist./beam_waist_evolve./exp(rsq./beam_waist.^2);     % spatial profile of gaussian beam
c2 = cos(K*zf-0.5*rsq./R + phi + omega*(t-zf/c));   % modulation with laser frequency
Ex = abs(c1.*c2.*envelop);

end




