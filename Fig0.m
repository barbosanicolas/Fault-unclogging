clear all; close all; 

%% 1- Compute fluid relative displacement
% %--------------------- MATRIX PROPERTIES ----------------------------------%
% Berea Sandstone properties representative of Candela 2014
ks=37e9;                      % [Pa]   --> [Grain modulus]
rhos=2.65*1e3;                % [kg/m3]--> [Grain density]
kappa=10^(-14.553);           % [m^2]     --> [Permeability] in Candela 2014

S=3;                          % tortuosity
nj=8;                         % [-]     --> [Permeability-related parameter]

% Get porosity from relation with permeability according to plot of Bourbie
phi=0.05;
% Pride relation to get elastic moduli from porosity
mus=44e9; % grain shear modulus
c=4; % consolidation parameter
km=ks*(1-phi)/(1+c*phi);        % White's model
mum=mus*(1-phi)/(1+3*c*phi/2);  % White's model

% In the case of Candela 2014, F0 is in the order of 0.1-0.5MPa/50mm
F0=0.2E6/50E-3;                 % [Pa/m]

%--------------------- FLUID PROPERTIES ----------------------------------%
%Water
kf=2.25e9;                              % [Pa]  --> [Fluid modulus]       
rhof=1.09e3;                            % [kg/m3] --> [Fluid density]
eta=16e-4;                              % [Pa*s]    --> [Viscosity] in Candela 2014

%--------------------- BIOT PARAMETERS  ----------------------------------%
rob=(phi.*rhof+ (1-phi).*rhos);         % [kg/m3] --> [Bulk density] 
alfa=1-km./ks;
kav=1./((alfa-phi)./ks+phi./kf);
H_d=km+4*mum/3;                     % Dry P-wave modulus
H_u=H_d + alfa^2*kav;               % Saturated P-wave modulus
omega_biot=eta*phi/kappa/rhof/S;    % Biot freq
f_biot=omega_biot/2/pi;             % Biot freq in Hz
f_drained_undrained=4*kappa*km/eta/(50e-3)^2;                % According to Pimienta et al 2015.
N=kav*H_d/H_u;
Diff=kappa*N/eta;                   % Hydraulic diffusivity

%--------------------- WAVE PARAMETERS  ----------------------------------%
% STRAIN-FREQ DEPENDENT ANALYSIS
f_fixed=0.05;                   % Wave frequency in Hz used by Candela 2014
omega=2*pi*f_fixed;             % omega
frec_vec=f_fixed;               % linear freq range in Hz

% Strain range
s0=5e-7;        
s1=5e-6;
strain_per_decade=100;              % sampling
strain_expon=(log10(s0):1/strain_per_decade:log10(s1)); % strain exponents
strain_wave=10.^strain_expon;       % strain imposed by the seismic wave. Variable


% Solution of Biot's equations
% Dynamic permeability
k_dyn=kappa./(sqrt(1+1i*4*omega/nj./omega_biot)+1i*omega./omega_biot); 
g=1./omega.*imag(eta./k_dyn);
b=real(eta./k_dyn);

%% S-Wave number computation. NOT USED
% S_s=(1/mum*(rob-rhof^2./(g-1i*b./omega))).^(1/2);   % slowness
% k_s=omega.*S_s;                                     % wave number S wave
% v_s=1./real(S_s);                                   % vel S wave
% Q_s=-2*imag(S_s)./real(S_s);                        % attenuation S wave
% gamma_s=omega*rhof.*k_dyn./(1i*eta);                % gamma factor S wave
% 
% for istr=1:length(strain_wave)
%     for ifre=1:length(frec_vec)
% B_s=[];
% 
% A_s=2*strain_wave(istr)./(k_s(ifre)^2); 
% B_s=A_s*gamma_s(ifre);                                   % amplitude of S wave
% 
% ws(istr,ifre)=1i*k_s(ifre)*B_s;                          % Fluid relative displ S wave
% us(istr,ifre)=1i*k_s(ifre)*A_s;                          % Solid displ S wave
% vss(istr,ifre)=us(istr,ifre)*1i*omega(ifre);              % Particle velocity displ S wave
% 
%     end
% end
%% Fast and slow P-Wave number computation 
a1=kav*H_d;
b1=-rob*kav - H_u*g + H_u*1i*b./omega + 2*rhof*alfa*kav;
c1=rob*g - rob*1i*b./omega - rhof^2;
gamma=-b1./2./a1;

for iom=1:length(omega)
sol=roots([a1, 0, b1(iom), 0, c1(iom)]);
S_aux=sol(imag(sol)<0);
[ai bi]=max(S_aux);[ci di]=min(S_aux);
S_p2(iom)=S_aux(bi);
S_p1(iom)=S_aux(di);
end

k_p1=omega.*S_p1;                                   % fast P wavenumber
k_p2=omega.*S_p2;                                   % slow P wavenumber
gamma_p1=(omega.^2.*rob - H_u.*k_p1.^2)./(-omega.^2.*rhof + alfa*kav.*k_p1.^2);
gamma_p2=(omega.^2.*rob - H_u.*k_p2.^2)./(-omega.^2.*rhof + alfa*kav.*k_p2.^2);

for istr=1:length(strain_wave)
    for ifre=1:length(frec_vec)
B_p1=[];B_p2=[];
        
A_p1=strain_wave(istr)/(-1i*k_p1(ifre));
B_p1=A_p1*gamma_p1(ifre);
A_p2=strain_wave(istr)/(-1i*k_p2(ifre));
B_p2=A_p2*gamma_p2(ifre);

wp1(istr,ifre)=B_p1;                       % Fluid relative displ fast P wave
wp2(istr,ifre)=B_p2;                       % Fluid relative displ slow P wave
up1(istr,ifre)=A_p1;                       % Solid displ fast P wave
up2(istr,ifre)=A_p2;                       % Solid displ slow P wave
vsp1(istr,ifre)=A_p1*1i*omega(ifre);       % Particle velocity fast P wave
vsp2(istr,ifre)=A_p2*1i*omega(ifre);       % Particle velocity slow P wave


Pfp1(istr,ifre)=(alfa*kav*A_p1 + kav*B_p1)*1i*k_p1(ifre);                       % Fluid pressure fast P wave
Pfp2(istr,ifre)=(alfa*kav*A_p2 + kav*B_p2)*1i*k_p2(ifre);                      % Fluid pressure slow P wave
grad_pf_p2(istr,ifre)=Pfp2(istr,ifre)*(-1i*k_p2(ifre));       % pressure gradient due to slow P-wave
grad_pf_p1(istr,ifre)=Pfp1(istr,ifre)*(-1i*k_p1(ifre));       % pressure gradient due to slow P-wave
sigma11_p1(istr,ifre)=-(H_u*A_p1 + alfa*kav*B_p1)*1i*k_p1(ifre);                       % Stress fast P wave
sigma11_p2(istr,ifre)=-(H_u*A_p2 + alfa*kav*B_p2)*1i*k_p2(ifre);                       % Stress slow P wave
    end
end
v_p1=1./real(S_p1);
Q_p1=-2*imag(S_p1)./real(S_p1);

v_p2=1./real(S_p2);
Q_p2=-2*imag(S_p2)./real(S_p2);


%% Permeability change prediction for CANDELA ET AL 2014
data=load('Candela_2014_fig4_intact.txt','-ASCII');
Pf_candela_norm=data(:,1);
deltaperm_candela_norm=data(:,2);
y=log(deltaperm_candela_norm);x=[ones(length(y),1) log(Pf_candela_norm)];
coeff=x\y;
b_coeff=coeff(2);a_coeff=exp(coeff(1));

normal_perm_p1=a_coeff.*(abs(grad_pf_p1)./abs(F0)).^b_coeff;   % Candela criterion
normal_perm_p2=a_coeff.*(abs(grad_pf_p2)./abs(F0)).^b_coeff;   % Candela criterion


%% Plot Fig0
figure(1)
set(1,'Units','inches','Position',[10 10 10 10],'PaperPositionMode','auto');
semilogy(Pf_candela_norm, deltaperm_candela_norm,'ok','MarkerFaceColor','r','MarkerSize',15);grid on;hold on
plot(abs(grad_pf_p2)./abs(F0), normal_perm_p2,'--r','MarkerSize',15,'LineWidth',3);
plot(abs(grad_pf_p1)./abs(F0), normal_perm_p1,'--k','MarkerSize',15,'LineWidth',3);
legend('Candela et al. (2014)', 'Theoretical Slow P-wave', 'Theoretical Fast P-wave','Location','SouthEast') 
ylabel('Relative permeability change');axis square
xlabel('Normalized pressure amplitude');grid on;axis tight
xlim([0 0.7]);ylim([1e-20 1])
set(gca, 'fontsize', 22)
% print(sprintf('%s','Candela_normalizedpf_permchange_bothwaves_otherstrainrange'), '-dpng', '-r300'); %<-Save as PNG with 300 DPI


% %%
% figure
% 
% loglog(strain_wave, abs(grad_pf_p2)./abs(F0))
% ylim([0.2 0.7])