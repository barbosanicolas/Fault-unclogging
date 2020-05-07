clear all;
%% Physical properties 
% Make properties a function of porosity
phi0=0.05;
phib=phi0*it;                                   % porosity 

c=4;                                            % consolidation parameter 
B=0.0009;                                        % geometric factor
dgrain=150*1e-6;                                 % grain diameter
kappa_b=B*dgrain^2*phib^3/(1-phib)^2;           % Kozeny-Carman
mus=44;                                         % grain shear modulus

ks=[37 37 37];
ros=[2.7 2.7 2.7];
phi=[phib 0.25 phib];
km=[ks(1)*(1-phib)/(1+c*phib) 0.08*ks(1)*(1-phib)/(1+c*phib) ks(1)*(1-phib)/(1+c*phib)];            %SOFT
mum=[mus*(1-phib)/(1+3*c*phib/2) 0.15*mus*(1-phib)/(1+3*c*phib/2) mus*(1-phib)/(1+3*c*phib/2)];     %SOFT
% km=[ks(1)*(1-phib)/(1+c*phib) 0.66*ks(1)*(1-phib)/(1+c*phib) ks(1)*(1-phib)/(1+c*phib)];            %STIFF
% mum=[mus*(1-phib)/(1+3*c*phib/2) 0.66*mus*(1-phib)/(1+3*c*phib/2) mus*(1-phib)/(1+3*c*phib/2)];     %STIFF
% km=[ks(1)*(1-phib)/(1+c*phib) 1.891 ks(1)*(1-phib)/(1+c*phib)];
% mum=[mus*(1-phib)/(1+3*c*phib/2) 5 mus*(1-phib)/(1+3*c*phib/2)];
kappa=[kappa_b/0.987e-12 0.5 kappa_b/0.987e-12];
kf=[2.25 2.25 2.25];
rof=[1 1 1];
eta=[1e-2 1e-2 1e-2];
nj=[8 8 8];
m=[1.57 1.57 1.57];
%%%%%%%%%%%%%%%%%%%%%%
%% Additional Physical parameters
%%%%%%%%%%%%%%%%%%%%%%

% Unit conversion
ks=ks.*1e6;
ros=ros.*1e6;
km=km.*1e6;
mum=mum.*1e6;
kappa=kappa.*0.987e-12;
kf=kf.*1e6;
rof=rof.*1e6;
eta=eta.*1e-1;

% Biot Parameters
alfa=1-km./ks;
kav=1./((alfa-phi)./ks+phi./kf);
rob=phi.*rof+ (1-phi).*ros;

L=km+4*mum/3;                       % Dry modulus
lambda_dry_lame=L - 2*mum;
C=L + alfa.^2.*kav;                 % Undrained modulus
N=kav.*L./C;
E=9*km.*mum./(3*km + mum);          % Young Modulus
% E(1)/E(2)
Diff=kappa.*N./eta;                 % Diff=1 m2/s is considered a high diffusivity for a fault

eff=kappa./eta./sqrt(Diff);
Diff_ef=eff(1)^2/(eff(2)^2+eff(1)*eff(2))*Diff(2);

%% %%%%% DEFINE FAULT THICKNESS !%%%%%
H_fault=.1;
%% %%%%%%

omega_biot=eta.*phi./kappa./m./rof;
f_biot=omega_biot/2/pi*1e3;         % in Hz

vp_faults=[vp_faults,sqrt(C(2)/rob(2))];
vs_faults=[vs_faults,sqrt(mum(2)/rob(2))];
rob_faults=[rob_faults,rob(2)*1E-3];
factorE=[factorE,E(1)/E(2)];
normal_compliances=[normal_compliances, 0.1/(L(2))/1E3];  %[m/Pa]
shear_compliances=[shear_compliances, 0.1/(mum(2))/1E3];  %[m/Pa]

file_code=load(sprintf('%s%i%s','Output_files/slow_wave_w_',it,'.txt'));

%% Analisis of induced fluid relative displacement

omegas=file_code(:,1)*2*pi/1e3;                     % read frequencies from file
Arp2=complex(file_code(:,2),file_code(:,3));        % read slow p wave reflection coefficient
n=complex(file_code(:,4),file_code(:,5));           % read horizontal component of wavenumber
lp1=complex(file_code(:,6),file_code(:,7));         % read vertical component of fast wavenumber
lp2=complex(file_code(:,8),file_code(:,9));         % read vertical component of slow wavenumber
gamma_p2=complex(file_code(:,10),file_code(:,11));  % read gamma factor for slow wave in the background
Up2=complex(file_code(:,12),file_code(:,13));       
Dp2=complex(file_code(:,14),file_code(:,15));       
lp2_m=complex(file_code(:,16),file_code(:,17));      % read vertical component of slow wavenumber in the fault
gamma_p2_m=complex(file_code(:,18),file_code(:,19)); % read gamma factor for slow wave in the fault
Atp2=complex(file_code(:,20),file_code(:,21));       % read slow p wave transmitted coefficient


lp2_background_reference(it,:)=lp2;
lp2_fault_reference(it,:)=lp2_m;


linear_freq_Hz=omegas*1e3/2/pi;                     % Compute linear frequency in Hz.

result_veff=[];result_wdot=[];
for istrain=[10 50 100]
    
A_i=[];w1_dot=[];w2_dot=[];eps22_inc=[];eps11_inc=[];eps12_inc=[];Elast_energ=[];
factor_w2_max=[];factor_veff2_max=[];

strain_desired=istrain*1e-7;                        % Desired incident strain
A_i=-strain_desired./lp1.^2;                        % Corresponding incident potential amplitude given desired strain
w1_dot=(1i*omegas).*(-1i.*n).*Arp2.*gamma_p2.*A_i;      % Horizontal component of Fluid relative velocity
w2_dot=(1i*omegas).*1i.*lp2.*Arp2.*gamma_p2.*A_i;       % Vertical component of Fluid relative velocity
eps22_inc=-(lp1.^2).*A_i;                           % Incident strain vertical comp.
eps11_inc=-(n.^2).*A_i;                             % Incident strain horizontal comp.
eps12_inc=-n.*lp1.*A_i;                             % Incident strain shear comp.

w2_dot_m_up=(1i*omegas).*1i.*lp2_m.*gamma_p2_m.*A_i.*Up2; % Vertical component of Fluid relative velocity (fault)
w2_dot_m_down=-(1i*omegas).*1i.*lp2_m.*gamma_p2_m.*A_i.*Dp2; % Vertical component of Fluid relative velocity (fault)
w2_dot_m_transm=-(1i*omegas).*1i.*lp2_m.*gamma_p2.*A_i.*Atp2; % Vertical component of Fluid relative velocity (fault)

Elast_energ_dens=(rob(1).*(eps22_inc.*omegas./2./lp1).^2)*1e3;   % Elastic energy density with unit conversion to [J/m3]

% Compute effective pore velocity with change of units from m/ms to mm/s
factor_veff2_max=w2_dot*1e6/phi(1); %[mm/s]
factor_veff2_max_fault_d=w2_dot_m_down*1e6/phi(2); %[mm/s]
factor_veff2_max_fault_u=w2_dot_m_up*1e6/phi(2); %[mm/s]
factor_veff2_max_back_transm=w2_dot_m_transm*1e6/phi(1); %[mm/s]


% Which radius gets unclogged?
pf_max_fault=(w2_dot_m_down+w2_dot_m_up)/kappa(2)*eta(2);
radio_pore_fault=sqrt(0.1./abs(pf_max_fault)*eta(2)*8/1e6);     %[m]
pf_max_back=w2_dot/kappa(1)*eta(1);
radio_pore_back=sqrt(0.1./abs(pf_max_back)*eta(1)*8/1e6);     %[m]

pf_max_fault_u=w2_dot_m_up/kappa(2)*eta(2);
pf_max_fault_d=w2_dot_m_down/kappa(2)*eta(2);


%% Make output array
result_veff=[result_veff;linear_freq_Hz, factor_veff2_max, strain_desired*ones(length(omegas),1),factor_veff2_max_fault_d,...
    Elast_energ_dens, factor_veff2_max_fault_u, factor_veff2_max_back_transm];

result_wdot=[result_wdot;linear_freq_Hz, strain_desired*ones(length(omegas),1), w2_dot, pf_max_back, pf_max_fault];

end
% Put all results in one array
overall_result_veff.(sprintf('%s','case_',num2str(it)))=result_veff;
overall_result_wdot.(sprintf('%s','case_',num2str(it)))=result_wdot;
