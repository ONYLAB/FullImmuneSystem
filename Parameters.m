function Parameters(SimType,epitopes,HLA_DR)

% The protein-specific parameters in this current code were used for
% simulating in vivo immune response in human against FVII.

% parameters

% NA: Avogadro constant
NA = 6.0221367e23;

%% Therapeutic protein PK parameters

% F_Bio: bioavailability (protein-specific)
F_Bio=1;   % no unit

DosePerPatientkg = 83; %ug/kg
Patientkg = 75;

% Vp: plasma volume
Vpperkg = 39e-3; %L/kg 
Vp = Vpperkg*Patientkg;%2.75;% L This is the Volume of the sample

% AMT: amount of protein dosed (protein-specific)
AMT=Patientkg*DosePerPatientkg*1e-3;%40; % mg

% MW: molecular weight of the protein (protein-specific)
MW=45.18e3; %Daltons = g/mol

% Fitted concentration in plasma 
fitconc = 20.88e3;

% Dose: bioavailabe drug dose (protein-specific)
Dose = fitconc*Vp; %pmol = pmol/L * L (before=>AMT*1E9/MW*F_Bio;  % pmole)

EndotoxinContentoftheDrug = 0.1; %ng/mg drug
Endotoxin = EndotoxinContentoftheDrug*AMT; %ng = (ng/mg * mg)
VzEndo = 54e-3; %L/kg Volume of distribution - Eritoran
Endotoxin = SimType*Endotoxin/VzEndo*Vpperkg; %Scaling for available amount in plasma compartment ng

% Vec: volume of extra central compartment (protein-specific)
Vec= 0.36920; % L

% kaAg: absorption rate constant from injection site to plasma for therapeutic protein (protein-specific)
kaAg=1e4;%0.28; % day-1

%kel: elimation rate of therapeutic protein (protein-specific)
kel=6.43;%0.11370; % day-1

% k12: distribution rate constant from the plasma to the extra central compartment
k21= 0.0;%10000; % day-1

% k21: distribution rate constant from the extra central compartment  to the plasma (protein-specific)
k12= 0.0;%Vec/Vp*k21 ; % day-1

% k13: distribution rate constant form plasma to peripheral compartment (protein-specific)
k13=0.0;%0.39930 ; % day-1

% k31: distribution rate constant form peripheral to plasma (protein-specific)
k31=0.0;%0.43256; % day-1

%% T-epitope characteristics of therapeutic proteins

% N: the number of T-epitope (protein-specific)
% Number of HLA molecules
numHLADR = 1e5;
[kon,koff,N,ME0] = doNETMHCIIpan(epitopes,HLA_DR,SimType,NA,numHLADR);

% ME0:  initial amount of MHC-II molecule in a single mature dendritic cell
% ME0=[95.1E3/2; 95.1E3/2;  34E3/2; 34E3/2; 17.1E3/2; 17.1E3/2]/NA*1E12; % pmole

% % p1: DR allele frequency distribution (patient population-specific)
% p1=[0.089473684, 0.053383459,0.036090226,0.084962406,0.146616541,0.008270677,0.069172932,...
%   0.001503759, 0.043609023,0.00075188,0.008270677, 0.457894736];

% % randomly generate two DR alleles for a virture patient depending on p1
% [~,a] = histc(rand(1,2),[0;cumsum(p1(:))/sum(p1)]);
%
%   % T-epitope-MHC-II binding affinity table (Kd, unit nM): (protein-specific)
%   % For epitope #1:
%   Affinity1_DR=[123;78.52;180;124.73;57.44;75;306;112.43;317;53.7;148;4000];
%   Affinity1_DP=[4000 4000];
%   Affinity1_DQ=[4000 4000];
%
%    % For epitope #2:
%   Affinity2_DR=[85;147.85;38;104.16;101.5;77;292;4000;293;4000;4000;4000];
%   Affinity2_DP=[4000 4000];
%   Affinity2_DQ=[4000 4000];
%
%   %	koff:	off rate for for T-epitope-MHC-II binding
%    koff=8.64*1E-3*[Affinity1_DR(a(1)) Affinity1_DR(a(2)) Affinity1_DP Affinity1_DQ;
%                               Affinity2_DR(a(1)) Affinity2_DR(a(2)) Affinity2_DP Affinity2_DQ]*1E3; %  day-1
%
% % kon: on rate for for T-epitope-MHC-II binding
% kon=ones(N,6)*8.64*1E-3; %  pM-1day-1

%% Dendritic cells

% BetaMS: decay rate for the maturation signals (LPS)
BetaMS=0.3696;  % day-1

% BetaID:	death rate of immature dendritic cells.
BetaID=0.0924; % day-1

% DeltaID:	maximum activation rate of immature dendritic cells
DeltaID=1.66;%1.5; % day-1

% Km,MS: LPS concentration at which immature DC activation rate is 50% maximum.
KMS= 9.852E3; % ng/L

% BetaMD:	death rate of mature dendritic cells
BetaMD=0.2310 ; % day-1

%% Antigen presentation

% konN: on rate for competing peptide-MHC-II binding
konN=0.0;%ones(6,1)*8.64*1E-3; %  pM-1day-1

% koffN: off rate for competing peptide-MHC-II binding
koffN=8.64*1E-3*ones(6,1)*4000E3; %  day-1

% AlphaAgE:  Ag internalization rate constant by mature dendritic cells
AlphaAgE=14.4; %day-1

%	BetaAgE	: degradation rate for AgE in acidic vesicles
BetaAgE= 17.28; % day-1

%	Beta_p:	degradation rate for T-epitope peptide (same for all peptides)
Beta_p= 14.4; % day-1

% BetaM: degradation rate for MHC-II
BetaM= 1.663; % day-1.

% Beta_pM: degradation rate for pME
Beta_pM= 0.1663; % day-1

% kext: exocytosis rate for peptide-MHC-II complex from endosome to cell menbrane
kext = 28.8;  % day-1

% kin: internalization rate for peptide-MHC-II complex on DC membrane
kin = 14.4;  % day-1

% KpM_N: number of peptide-MHCII to achieve 50% activation rate of naïve helper T cells
KpM_N = 400;  % number of complex

%  KpM_M: number of peptide-MHCII to achieve 50% activation rate of memory helper T cells
KpM_M = 40;   % number of complex

% VD: volume of a dendritic cell
VD=2.5388e-012; % L

% VE: volume of endosomes in a dendritic cell
VE=4E-16; % L

%BetaNT: death rate of naive helper T cells
BetaNT=0.0029 ; % day-1

%DeltaNT: maximum activation rate of naive helper T cells
DeltaNT=1.5; % day-1

% RhoAT: maximum proliferation rate for activated helper T cells
RhoAT=0.5973; % day-1

%BetaAT: death rate of activated helper T cells
BetaAT=0.18; % day-1

%DeltaMT:	maximum activation rate of memory T cells
DeltaMT=1.5; % day-1

% BetaMT: death rate of memory helper T cells
BetaMT=2.7397e-004; % day-1

% BetaFT:	death rate of functional helper T cells
BetaFT=0.18; % day-1

% f1: differentiation percentage for activated helper T cells to become memory T cells
f1=0.5;
%% B cells

% J: the number of B cell subclones
J=17;

%   Ka(1:J): Association rate constant for Ag-BCR/Ag-Ab binding
for j = 1:17
    Ka(j,1)=1e-6*2^(j-(J+1)/2); % Ka is in the unit of pM-1(L/pmole)
end

% BRN: number of  B cell receptor on each B cell, including naïve B cells, activated B cells, as well
% as memory B cell
BRN=75000;  % number of BCR/cell

% DeltaNB:	maximum activation rate of naive B cells.
DeltaNB=3.0; % day-1

%	K_R: occupied BCR number to achieve 50% activation rate for naïve and memory B cells
K_R=1.0;  % number of BCR

% CC_N: the carrying capacity for 1 FT cell to stimulate the activation and proliferation of target naive B cells.
CC_N=10; % no unit

% CC_M: the carrying capacity for 1 FT cell to stimulate the activation and proliferation of target memory B cells.
CC_M=100; % no unit

%  RhoAB_N: maximum proliferation rate for activated B cells derived from naive B cells.
RhoAB_N=0.3333; % day-1

%  RhoAB_M: maximum proliferation rate for activated B cells derived from memory B cells.
RhoAB_M=0.7273; % day-1

%	BetaAB:	death rate of activated B cells
BetaAB=0.2518; % day-1

% DeltaMB: maximum activation rate of memory B cells
DeltaMB=3.0; % day-1

% BetaMB: death rate of memory B cells
BetaMB=7.8278e-005; % day-1

% g1: percentage for activated B cells to differentiate to memory B cells
g1=0.5; % no unit

% g2:	percentage for activated B cells to differentiate to short-lived plasma cells
g2=0.4; % no unit

%	 BetaSP: death rate of  short-lived plasma cells
BetaSP=0.2310; % day-1

% BetaLP: death rate of  long-lived plasma cells
BetaLP=0.0050; % day-1

% AlphaA: secretion rate of antibody by plasma cells
AlphaA=8.64E8; % day-1

%% Ab and immune complex

% BetaA: elimination rate of antibody
BetaA= 0.030130435; % day-1

% BetaC: elimation rate of immune complex
BetaC=10*kel; % day-1

%% Initial conditions for state variables in the differential equations:
% AgIS0: initial Ag in the injection site
AgIS0=Dose; % pmole

% Ag0: intitial therapeutic protein in the plasma compartment
Ag0=0; % pmole

% Agec0: initial therapeutic protein in the extra central compartment
Agec0=0; % pmole

% Agp0:initial Ag in the peripheral compartment
Agp0=0; % pmole

% MS0: Maturation signal (MS), particularly, endotoxin, LPS
MS0=Endotoxin;%5*70; % ng

% ID0: the initial number of immature dendritc cell
ID0=5e7; % cells

% MD0: the initial number of mature dendritic cells
MD0=0; % cells

% cp0: amount of endogenous competing protein in the plasma
cp0=0.0;%110E6*Vp; % pmole

% cpE0: initial amount of endogenous competing protein in endosome
cpE0=0;   % pmole

% cptE0: initial amount of endogenous competing peptide in endosome
cptE0=0;   % pmole

% cptME0: initial amount of endogenous peptide-MHC complex in endosome
cptME0=zeros(6,1); % pmole

% cptM0: initial amount of endogenous peptide-MHC complex on dendritic cell membrane
cptM0=zeros(6,1);  % pmole

% AgE0: initial amount of Ag in endosome
AgE0=0; % pmole

% pE0:	initial amount of T-epitope peptides from Ag digestion in endosome
pE0=ones(N,1)*0;  % pmole

% pME0:	initial amount of T-epitope-MHC-II complex in endosome
pME0=ones(6*N,1)*0;  % pmole

% pM0:	initial amount of T-epitope-MHC-II complex on dendritic cell membrane
pM0=ones(6*N,1)*0;  % pmole

% M0, free MHC-II molecule on dendritic cell membrane
M0=ones(6,1)*0; % pmole

% NT0: the initial number of naive T helper cell = total T helper cells * frequency
% of T-epitope-specific T helper cells (protein-specific)
% Human: Total naive T helper cells per kg body weight is: 876E6 cells/L*5L(blood
% volume)= 4.3800e+009 cells
NT0=ones(N,1)* 4.3800e+009*0.1/1E6; % cells changed for Factor VII

% AT_N0: initial number of activated T helper cell derived from naive T cells
AT_N0=ones(N,1)*0.0; % cells

% AT_M0: initial number of activated T helper cell derived from memory T cells
AT_M0=ones(N,1)*0.0; % cells

% MT0: initial number of memory T helper cell
MT0=ones(N,1)*0.0; % cells

% FT0: initial number of functional T helper cell
FT0=ones(N,1)*0.0; % cells

% NB0: initial number of naïve B cells
NB0=(normpdf(-3:6/16:3)./sum(normpdf(-3:6/16:3))*1.3e9*0.000004)'; % cells

% AB_N0: initial number of activated B cells derived from naive B cells
AB_N0=ones(J,1)*0.0; % cells

% AB_M0: initial number of activated B cells derived from memory B cells
AB_M0=ones(J,1)*0.0; % cells

% SP0: initial number of short-lived plasma (antibody secreting) B cells
SP0=ones(J,1)*0.0; % cells

%LP0: initial number of long-lived (antibody secreting) B cells
LP0=ones(J,1)*0.0; % cells

% BM0: initial number of memory B cells
MB0=ones(J,1)*0.0; % cells

% A0: initial total amount of antibody
A0=ones(J,1)*0.0; % pmole

%% Parameter vector
pars(1)=NA;
% Therapeutic protein PK parameters
pars(2)=F_Bio;
pars(3)=Dose;
pars(4)=kaAg;
pars(5)=kel;
pars(6)=k21;
pars(7)=k12;
pars(8)=k13;
pars(9)=k31;
pars(10)=Vp;
pars(11)=Vec;
% T-epitope characteristics of therapeutic proteins
pars(12)=N;
pars(13:(12+6*N))=reshape(kon,6*N,1);
pars((13+6*N):(12+12*N))= reshape(koff, 6*N,1);
% Dendritic cells
pars(13+12*N)=BetaMS;
pars(14+12*N)=BetaID;
pars(15+12*N)=DeltaID;
pars(16+12*N)=KMS;
pars(17+12*N)=BetaMD;
% Antigen presentation
pars((18+12*N):(23+12*N))=konN;
pars((24+12*N):(29+12*N))=koffN;
pars(30+12*N)=AlphaAgE;
pars(31+12*N)=BetaAgE;
pars(32+12*N)=Beta_p;
pars(33+12*N)=BetaM;
pars(34+12*N)=Beta_pM;
pars(35+12*N)=kext;
pars(36+12*N)=kin;
pars(37+12*N)=KpM_N;
pars(38+12*N)=KpM_M;
pars(39+12*N)=VD;
pars(40+12*N)=VE;
% T helper cells
pars(41+12*N)=BetaNT;
pars(42+12*N)=DeltaNT;
pars(43+12*N)=RhoAT;
pars(44+12*N)=BetaAT;
pars(45+12*N)=DeltaMT;
pars(46+12*N)=BetaMT;
pars(47+12*N)=BetaFT;
pars(48+12*N)=f1;
% B cells
pars(49+12*N)=J;
pars((50+12*N):(49+12*N+J))=Ka;
pars(50+12*N+J)=BRN;
pars(51+12*N+J)=DeltaNB;
pars(52+12*N+J)=K_R;
pars(53+12*N+J)=CC_N;
pars(54+12*N+J)=CC_M;
pars(55+12*N+J)=RhoAB_N;
pars(56+12*N+J)=RhoAB_M;
pars(57+12*N+J)=BetaAB;
pars(58+12*N+J)=DeltaMB;
pars(59+12*N+J)=BetaMB;
pars(60+12*N+J)=g1;
pars(61+12*N+J)=g2;
pars(62+12*N+J)=BetaSP;
pars(63+12*N+J)=BetaLP;
pars(64+12*N+J)=AlphaA;
% Ab and immune complex
pars(65+12*N+J)=BetaA;
pars(66+12*N+J)=BetaC;
% Initial conditions as parameters in the differential equations:
pars(67+12*N+J)=ID0;
pars(68+12*N+J)=cp0;
pars((69+12*N+J):(74+12*N+J))=ME0;
pars((75+12*N+J):(74+13*N+J))=NT0;

pars=pars';

save Parameters.mat;

% function [kon,koff,N,ME0] = doNETMHCIIpan(epitopes,HLA_DR,SimType,NA,numHLADR)
% %% Collect -on, -off rates, number of epitopes and amount of initial MHC
% % molecules
% EpitopePresent = 1;
% N = length(epitopes);
% if N<1
%     EpitopePresent = 0;
%     epitopes{1} = 'AAAAAAAAAAA'; %Dummy PolyA epitope
%     N=1;
% end
% 
% if strcmp(HLA_DR{1,1},HLA_DR{1,2}) %homozygot
%     disp('Homozygote');
%     ME0=[numHLADR; 0.0;  34E3/2; 34E3/2; 17.1E3/2; 17.1E3/2]/NA*1E12; % pmole
%     % kon: on rate for for T-epitope-MHC-II binding
%     kon=EpitopePresent*SimType*repmat([1 0 0 0 0 0],N,1)*8.64*1E-3; %  pM-1day-1
%     HLAtext = [HLA_DR{1,1} ',' HLA_DR{1,2}];
% else
%     ME0=[numHLADR/2; numHLADR/2;  34E3/2; 34E3/2; 17.1E3/2; 17.1E3/2]/NA*1E12; % pmole
%     % kon: on rate for for T-epitope-MHC-II binding
%     kon=EpitopePresent*SimType*repmat([1 1 0 0 0 0],N,1)*8.64*1E-3; %  pM-1day-1
%     HLAtext = [HLA_DR{1,1} ',' HLA_DR{1,2}];
% end
% 
% Affinity_DPQ=[4000 4000 4000 4000]; %place holder for other alleles (NOT USED)
% Affinity_DR = givekon(epitopes{1},HLAtext);
% %	koff:	off rate for for T-epitope-MHC-II binding
% koff=8.64*1E-3*[Affinity_DR Affinity_DPQ]*1E3; %  day-1
% 
% for i = 2:N
%     Affinity_DR = givekon(epitopes{i},HLAtext);
%     temp = 8.64*1E-3*[Affinity_DR Affinity_DPQ]*1E3;
%     koff = [koff;temp];
% end
% 
% 
% function Affinity_DR = givekon(epitopesequence,HLAtext)
% %% Run NetMHCIIPan if not done before, extract association rates
% % data(1).Sequence = 'epitopesequence'
% % data(1).Header = 'Seq'
% % fastawrite('example.fsa',data);
% if exist('out.dat')==0 %#ok<EXIST>
%     dlmwrite('example.fsa',epitopesequence,'')
%     command = ['netMHCIIpan -f example.fsa -inptype 1 -xls -xlsfile out.dat -a ' HLAtext];
%     [s,m] = unix(command);
% end
% table = readtable('out.dat');
% Affinity_DR(1) = table{1,'nM'};
% Affinity_DR(2) = table{1,'nM_1'};

function [kon,koff,N,ME0] = doNETMHCIIpan(epitopes,HLA_DR,SimType,NA,numHLADR)
%% Collect -on, -off rates, number of epitopes and amount of initial MHC
% molecules
EpitopePresent = 1;
N = length(epitopes);
if N<1
    EpitopePresent = 0;
    epitopes{1} = 'AAAAAAAAAAA'; %Dummy PolyA epitope
    N=1;
end

if strcmp(HLA_DR{1,1},HLA_DR{1,2}) %homozygot
    disp('Homozygote');
    nallele = 1;
    ME0=[numHLADR; 0.0;  34E3/2; 34E3/2; 17.1E3/2; 17.1E3/2]/NA*1E12; % pmole
    % kon: on rate for for T-epitope-MHC-II binding
    kon=EpitopePresent*SimType*repmat([1 0 0 0 0 0],N,1)*8.64*1E-3; %  pM-1day-1
    HLAtext = [HLA_DR{1,1} ',' HLA_DR{1,2}];
else
    nallele = 2;
    ME0=[numHLADR/2; numHLADR/2;  34E3/2; 34E3/2; 17.1E3/2; 17.1E3/2]/NA*1E12; % pmole
    % kon: on rate for for T-epitope-MHC-II binding
    kon=EpitopePresent*SimType*repmat([1 1 0 0 0 0],N,1)*8.64*1E-3; %  pM-1day-1
    HLAtext = [HLA_DR{1,1} ',' HLA_DR{1,2}];
end

Affinity_DPQ=[4000 4000 4000 4000]; %place holder for other alleles (NOT USED)
Affinity_DR = givekon(epitopes{1},nallele,HLA_DR);
%	koff:	off rate for for T-epitope-MHC-II binding
koff=8.64*1E-3*[Affinity_DR Affinity_DPQ]*1E3; %  day-1

for i = 2:N
    Affinity_DR = givekon(epitopes{i},nallele,HLA_DR);
    temp = 8.64*1E-3*[Affinity_DR Affinity_DPQ]*1E3;
    koff = [koff;temp];
end

function Affinity_DR = givekon(epitopesequence,nallele,HLAtext)
%% Run NetMHCIIPan if not done before, extract association rates
% data(1).Sequence = 'epitopesequence'
% data(1).Header = 'Seq'
% fastawrite('example.fsa',data);
if exist('IEDB_KD.dat')==0 %#ok<EXIST>
%     dlmwrite('example.fsa',epitopesequence,'')
%     command = ['netMHCIIpan -f example.fsa -inptype 1 -xls -xlsfile out.dat -a ' HLAtext];
%     [s,m] = unix(command);
    IEDBScraper(epitopesequence,nallele,HLAtext);
end
Affinity_DR = dlmread('IEDB_KD.dat');
Affinity_DR = Affinity_DR(:)';