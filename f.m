function dydt=f(t,y,pars)

global  solutionx
%% Parameter vector
NA=pars(1);
% Therapeutic protein PK parameters
F_Bio=pars(2);
Dose=pars(3);
kaAg=pars(4);
kel=pars(5);
k21=pars(6);
k12=pars(7);
k13=pars(8);
k31=pars(9);
Vp=pars(10);
Vec=pars(11);
% T-epitope characteristics of therapeutic proteins
N=pars(12);
kon=reshape(pars(13:(12+6*N)),N,6);
koff= reshape(pars((13+6*N):(12+12*N)), N,6);
% Dendritic cells
BetaMS=pars(13+12*N);
BetaID=pars(14+12*N);
DeltaID=pars(15+12*N);
KMS=pars(16+12*N);
BetaMD=pars(17+12*N);
% Antigen presentation
konN=pars((18+12*N):(23+12*N));
koffN=pars((24+12*N):(29+12*N));
AlphaAgE=pars(30+12*N);
BetaAgE=pars(31+12*N);
Beta_p=pars(32+12*N);
BetaM=pars(33+12*N);
Beta_pM=pars(34+12*N);
kext=pars(35+12*N);
kin=pars(36+12*N);
KpM_N=pars(37+12*N);
KpM_M=pars(38+12*N);
VD=pars(39+12*N);
VE=pars(40+12*N);
% T helper cells
BetaNT=pars(41+12*N);
DeltaNT=pars(42+12*N);
RhoAT=pars(43+12*N);
BetaAT=pars(44+12*N);
DeltaMT=pars(45+12*N);
BetaMT=pars(46+12*N);
BetaFT=pars(47+12*N);
f1=pars(48+12*N);
% B cells
J=pars(49+12*N);
Ka=pars((50+12*N):(49+12*N+J));
BRN=pars(50+12*N+J);
DeltaNB=pars(51+12*N+J);
K_R=pars(52+12*N+J);
CC_N=pars(53+12*N+J);
CC_M=pars(54+12*N+J);
RhoAB_N=pars(55+12*N+J);
RhoAB_M=pars(56+12*N+J);
BetaAB=pars(57+12*N+J);
DeltaMB=pars(58+12*N+J);
BetaMB=pars(59+12*N+J);
g1=pars(60+12*N+J);
g2=pars(61+12*N+J);
BetaSP=pars(62+12*N+J);
BetaLP=pars(63+12*N+J);
AlphaA=pars(64+12*N+J);
% Ab and immune complex
BetaA=pars(65+12*N+J);
BetaC=pars(66+12*N+J);
% Initial conditions as parameters in the differential equations:
ID0=pars(67+12*N+J);
cp0=pars(68+12*N+J);
ME0=pars((69+12*N+J):(74+12*N+J));
NT0=pars((75+12*N+J):(74+13*N+J));

%% State variables for the differential equations
% AgIS, antigenic protein in the injection site
AgIS=y(1);
% Ag, antigenic protein in the plasma compartment
Ag=y(2);
% Agec, antigenic protein in the extra central compartment
Agec=y(3);
% Agp, antigenic protein in the peripheral compartment
Agp=y(4);
% MS, maturation signal for Immature dendritic cells
MS=y(5);
% ID,	Immature dendritic cells
ID=y(6);
% MD,	mature dendritic cells
MD=y(7);
% cpE, endogenous competing protein in the endosomes
cpE=y(8);
% cptE, endogenous competing peptide in the endosomes
cptE=y(9);
% cptME, competing peptide-MHC complex in the endosomes
cptME=y(10:15);
% cptM, competing peptide-MHC complex on dendritic cell membrane
cptM=y(16:21);
%AgE, Ag in the endosomes
AgE=y(22);
% pE,	free epitope peptides from Ag digestion
pE=y((22+1):(22+N));
% ME(1:6),	free MHC II molecule in endosome
ME=y((23+N):(28+N));
% pME(N,6),	epitope peptide-MHC II complex in endosomes
% pME is in a matrix form, N epitopes against 6 possible MHC alleles
pME=reshape(y((29+N):(28+7*N)),N,6);
% pM(N,6),	epitope peptide-MHC II complex on dendritic cell membrane
pM=reshape(y((29+7*N):(28+13*N)),N,6);
% M, free MHC II molecule on dendritic cell menbrane
M=y((29+13*N):(34+13*N));
% NT, naïve helper T cells
NT=y((35+13*N):(34+14*N));
% AT_N	activated helper T cells
AT_N=y((35+14*N):(34+15*N));
% AT_M	activated helper T cells
AT_M=y((35+15*N):(34+16*N));
% MT, memory helper T cells
MT=y((35+16*N):(34+17*N));
% FT, functional helper T cells
FT=y((35+17*N):(34+18*N));
% NB,	naïve B cells
NB=y((35+18*N):(34+18*N+J));
% AB_N,	activated B cells
AB_N=y((35+18*N+J):(34+18*N+2*J));
% AB_M,	activated B cells
AB_M=y((35+18*N+2*J):(34+18*N+3*J));
% SP, short-lived plasma (antibody secreting) B cells
SP=y((35+18*N+3*J):(34+18*N+4*J));
% LP, long-lived plasma (antibody secreting) B cells
LP=y((35+18*N+4*J):(34+18*N+5*J));
% MB, memory B cells
MB=y((35+18*N+5*J):(34+18*N+6*J));
% A: antibody concentration (including free and bound antibody)
A= y((35+18*N+6*J):(34+18*N+7*J));

%% Calculate functions for helper T cells activation, proliferation, or differentiation

% Adding up all the pM molecules from one epitope against 6 different MHC alleles

pM_NUMBER=pM*ones(6,1)*1E-12*NA;

% The activation function D for helper T cells
D_N=(MD/(MD+sum(NT)+sum(AT_N)+sum(AT_M)+sum(MT)))*(pM_NUMBER./(pM_NUMBER+KpM_N));   % for naive T cells
D_M=(MD/(MD+sum(NT)+sum(AT_N)+sum(AT_M)+sum(MT)))*(pM_NUMBER./(pM_NUMBER+KpM_M));  % for memory T cells

% The proliferation/differentiation function E for helper T cells
E_N=(MD/(MD+sum(NT)+sum(AT_N)+sum(AT_M)+sum(MT)))*(pM_NUMBER-KpM_N)./(pM_NUMBER+KpM_N);  % for naive T cells
E_M=(MD/(MD+sum(NT)+sum(AT_N)+sum(AT_M)+sum(MT)))*(pM_NUMBER-KpM_M)./(pM_NUMBER+KpM_M); % for memory T cells

%% Calculate functions for B cells activation, proliferation, or differentiation

% Calculate the number of occupied BCR for B cell activation
BCR=(NB+AB_N+AB_M+MB)/NA*1E12*BRN; % BCR is the amount of BCR, pmole
Eq1 = @(x)(Ag/Vp-x*(1+2*sum(Ka.*A/Vp./(1+Ka*x))+ sum(Ka.*BCR/Vp./(1+Ka*x))));
Options = optimset('TolX',5e-16);
if t==0
    Agfree_con = fzero(Eq1,Ag/Vp,Options);
else
    Agfree_con = fzero(Eq1,solutionx,Options);
end
solutionx=Agfree_con;

% Calculate the receptor occupancy "ro"
ro = Ka*Agfree_con./ (1 + Ka*Agfree_con);
% Calculate occupied BCR number "R"
R=ro.*BRN;

% Calculate the functions based on BCR occupancy for B cell activation, proliferation or differentiation
F=R./(K_R+R);
G=(1-ro).*F;
H=(R-K_R)./(R+K_R);

% Functions based on help from functional T cells, for B cell activation, proliferation or differentiation
% function for naive B cells
P_N=CC_N*sum(FT)./(CC_N*sum(FT)+sum(NB)+sum(AB_N)+sum(AB_M)+sum(MB));

% function for memory B cells
P_M=CC_M*sum(FT)./(CC_M*sum(FT)+sum(NB)+sum(AB_N)+sum(AB_M)+sum(MB));

%% Differential equations
% AgIS, y(1), total amount of antigenic protein in the injection site, pmole
dydt(1,1)=-kaAg*AgIS;

% Ag, y(2), total amount of antigenic protein in the plasma, including free and bound, pmole
dydt(2,1)=kaAg*AgIS -(kel+k13+k12)*Agfree_con*Vp+k21*Agec+k31*Agp-kin*sum(ro(1:J).*BCR(1:J))...
    -BetaC*sum(ro(1:J)*2.*A(1:J));

% Agec,  y(3), antigenic protein in the extra central compartment, pmole
dydt(3,1)=k12*Agfree_con*Vp-k21*Agec;

% Agp, y(4), Amount of Ag in the peripheral compartment, pmole
dydt(4,1)=k13*Agfree_con*Vp-k31*Agp;

% MS, y(5), maturation signal, particularly, LPS, for immature dendritic cells, ng
dydt(5,1)=-BetaMS*MS;

% ID, y(6), immature dendritic cells, cells
dydt(6,1)=BetaID*(ID0-ID)-DeltaID*ID*(MS/Vp)/((MS/Vp)+KMS);

% MD,	y(7), mature dendritic cells , cells
dydt(7,1)=DeltaID*ID*(MS/Vp)/((MS/Vp)+KMS)-BetaMD*MD;

% cpE, y(8), endogenous competing protein in endosome, pmole
dydt(8,1)=AlphaAgE*VD*cp0/Vp -BetaAgE*cpE;

% cptE, y(9), endogenous competing peptide in endosome, pmole
dydt(9,1)=BetaAgE*cpE -Beta_p*cptE -(konN.*(cptE*(ME/VE)))'*ones(6,1)+(koffN.*cptME)'*ones(6,1);

% cptME, y(10:15), endogenous competing peptide-MHC complex in endosome, pmole
dydt((10:15),1)=konN*cptE.*(ME/VE)-koffN.*cptME -Beta_pM*cptME -kext*cptME;

% cptM, y(16:21), endogenous peptide-MHC complex on dendritic cell membrane, pmole
dydt((16:21),1)=kext*cptME -koffN.*cptM;

%AgE, y(22), Ag in the endosome, pmole
dydt(22,1)=AlphaAgE*VD*(Ag/Vp) -BetaAgE*AgE;

% pE, y(23:(22+N)), T-epitope peptides from Ag digestion, pmole
dydt((23:(22+N)),1)=BetaAgE*AgE  -Beta_p*pE...
    -kon.*(pE*(ME'/VE))*ones(6,1)+koff.*pME*ones(6,1);

% ME, y((23+N):(28+N)),	free MHC-II molecule in endosome, pmole
dydt((23+N):(28+N),1)= BetaM*(ME0-ME) -(ones(1,N)*(kon.*(pE*(ME'/VE))))'+(ones(1,N)*(koff.*pME))'...
    -konN*cptE.*(ME/VE)+koffN.*cptME +kin*M;

% pME, y((29+N):(28+7*N)), T-epitope-MHC-II complex in endosome, pmole
dydt(((29+N):(28+7*N)),1)= reshape(kon.*(pE*(ME'/VE))-koff.*pME  -Beta_pM*pME -kext*pME,6*N,1);

% pM,	y((29+7*N):(28+13*N)), T-epitope-MHC-II complex on dendritic cell membrane, pmole
dydt((29+7*N):(28+13*N),1)=reshape(kext*pME -koff.*pM,6*N,1);

% M, y((29+13*N):(34+13*N)), free MHC-II molecule on dendritic cell menbrane, pmole
dydt((29+13*N):(34+13*N),1)=-kin*M +(ones(1,N)*(koff.*pM))'+koffN.*cptM;

% NT, y((35+13*N):(34+14*N)), naïve helper T cells specific for T-epitope 1 to N,  cells
dydt((35+13*N):(34+14*N),1)=BetaNT*(NT0-NT) -DeltaNT*D_N.*NT;

%AT_N, y((35+14*N):(34+15*N)),	activated helper T cells derived from NT,  cells
dydt((35+14*N):(34+15*N),1)=DeltaNT*D_N.*NT +RhoAT*E_N.*AT_N -BetaAT*AT_N;

% AT_M, y((35+15*N):(34+16*N)),	activated helper T cells derived from MT, cells
dydt((35+15*N):(34+16*N),1)=DeltaMT*D_M.*MT +RhoAT*E_M.*AT_M -BetaAT*AT_M;

% MT, y((35+16*N):(34+17*N)), memory T helper cell, cells
dydt((35+16*N):(34+17*N),1)=RhoAT*(1-E_N)*f1.*AT_N +RhoAT*(1-E_M)*f1.*AT_M...
    -DeltaMT*D_M.*MT -BetaMT*MT;

%FT, y((35+17*N):(34+18*N)), functional helper T cell, cells
dydt((35+17*N):(34+18*N),1)=RhoAT*(1-E_N)*(1-f1).*AT_N +RhoAT*(1-E_M)*(1-f1).*AT_M -BetaFT*FT;

% NB, y((35+18*N):(34+18*N+J)), naïve B cells,   cells
dydt(((35+18*N):(34+18*N+J)),1)=  -DeltaNB*F*P_N.*NB;

% AB_N,	y((35+18*N+J):(34+18*N+2*J)), activated B cells derived from naive B cells,  cells
dydt(((35+18*N+J):(34+18*N+2*J)),1)=DeltaNB*G*P_N.*NB +  RhoAB_N*H*P_N.*AB_N  -BetaAB*AB_N;

% AB_M, y((35+18*N+2*J):(34+18*N+3*J)), activated B cells derived from memory B cells, cells
dydt(((35+18*N+2*J):(34+18*N+3*J)),1)=DeltaMB*G*P_M.*MB +RhoAB_M*H*P_M.*AB_M -BetaAB*AB_M;

% SP, y((35+18*N+3*J):(34+18*N+4*J)), short-lived plasma (antibody secreting) B cells,  cells
dydt(((35+18*N+3*J):(34+18*N+4*J)),1)=RhoAB_N*(1-H)*P_N*g2.*AB_N +RhoAB_M*(1-H)*P_M*g2.*AB_M -BetaSP*SP;

% LP,	y((35+18*N+4*J):(34+18*N+5*J)), long-lived (antibody secreting) B cells , cells
dydt(((35+18*N+4*J):(34+18*N+5*J)),1)=RhoAB_N*(1-H)*P_N*(1-g1-g2).*AB_N...
    +RhoAB_M*(1-H)*P_M*(1-g1-g2).*AB_M -BetaLP*LP;

% MB, y((35+18*N+5*J):(34+18*N+6*J)), memory B cells, cells
dydt(((35+18*N+5*J):(34+18*N+6*J)),1)=RhoAB_N*(1-H)*P_N*g1.*AB_N +RhoAB_M*(1-H)*P_M*g1.*AB_M...
    -DeltaMB*F*P_M.*MB -BetaMB*MB ;

% A, y((35+18*N+6*J):(34+18*N+7*J)), total antibody, pmole
dydt(((35+18*N+6*J):(34+18*N+7*J)),1)=AlphaA*(SP/NA*1E12) + AlphaA*(LP/NA*1E12)...
    -BetaA*(1-2*ro+ro.^2).*A -BetaC*(2*ro-ro.^2).*A ;
