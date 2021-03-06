
function Main_human(SimType,epitopes,HLA_DR,regimen)

% clc
close all

%% Main simulation

% Set a global variable to store the previous solution for free Ag concentration
global solutionx

numdoses = length(regimen)-1;

% Load the parameters
Parameters(SimType,epitopes,HLA_DR); %SimType=1 if with Sample, 0 if without
load Parameters.mat; %#ok<LOAD>

% Run the ODEs
options = odeset('RelTol',1e-10, 'AbsTol',1e-10);

% First dose
% Initial condition vector
yic1 = [AgIS0;Ag0;Agec0;Agp0; MS0; ID0; MD0; cpE0; cptE0;cptME0;cptM0;AgE0;pE0; ME0;pME0; pM0; M0;NT0; AT_N0; AT_M0;MT0; FT0;NB0; AB_N0; AB_M0; SP0; LP0;MB0;A0];

% First dosing interval, protein specific
tspan1 =linspace(regimen(1),regimen(2),100);

% Call ODE
[T1,Y1]=ode15s(@f, tspan1, yic1, options,pars);

% Record the results into new matrix
t_record=T1;
y_record=Y1;


% Subsequent doses
for num_dose=2:numdoses  % number of doses, protein specific
    
    % Initial conditions
    AgIS1=y_record(numel(t_record), 1)+Dose;
    Ag1=y_record(numel(t_record),2);
    Agce1=y_record(numel(t_record),3);
    Agp1=y_record(numel(t_record),4);
    MS1=y_record(numel(t_record),5)+MS0;
    ID1=y_record(numel(t_record),6);
    MD1=y_record(numel(t_record),7);
    cpE1=y_record(numel(t_record),8);
    cptE1=y_record(numel(t_record),9);
    cptME1=y_record(numel(t_record),(10:15));
    cptM1=y_record(numel(t_record),(16:21));
    AgE1=y_record(numel(t_record),22);
    pE1=y_record(numel(t_record),(22+1):(22+N));
    ME1=y_record(numel(t_record),(23+N):(28+N));
    pME1=y_record(numel(t_record),(29+N):(28+7*N));
    pM1=y_record(numel(t_record),(29+7*N):(28+13*N));
    M1=y_record(numel(t_record),(29+13*N):(34+13*N));
    NT1=y_record(numel(t_record),(35+13*N):(34+14*N));
    AT_N1=y_record(numel(t_record),(35+14*N):(34+15*N));
    AT_M1=y_record(numel(t_record),(35+15*N):(34+16*N));
    MT1=y_record(numel(t_record),(35+16*N):(34+17*N));
    FT1=y_record(numel(t_record),(35+17*N):(34+18*N));
    NB1=y_record(numel(t_record),(35+18*N):(34+18*N+J));
    AB_N1=y_record(numel(t_record), (35+18*N+J):(34+18*N+2*J));
    AB_M1=y_record(numel(t_record), (35+18*N+2*J):(34+18*N+3*J));
    SP1=y_record(numel(t_record), (35+18*N+3*J):(34+18*N+4*J));
    LP1=y_record(numel(t_record), (35+18*N+4*J):(34+18*N+5*J));
    MB1=y_record(numel(t_record), (35+18*N+5*J):(34+18*N+6*J));
    A1= y_record(numel(t_record), (35+18*N+6*J):(34+18*N+7*J));
    
    
    % Initial condition vector
    yic2 = [AgIS1,Ag1, Agce1,Agp1,MS1, ID1, MD1, cpE1, cptE1,cptME1,cptM1,AgE1,pE1,  ME1,pME1,  pM1, M1,NT1, AT_N1,AT_M1,MT1,FT1,NB1, AB_N1, AB_M1, SP1, LP1,MB1,A1]';
    
    % dosing interval
    t_start=regimen(num_dose);
    t_end=regimen(num_dose+1);
    tspan2 = linspace(t_start, t_end, 100);
    
    % Call ODE
    [T2,Y2]=ode15s(@f, tspan2, yic2, options,pars);
    
    % Record the result into new matrix
    t_record=[t_record;T2];
    y_record=[y_record; Y2];
    
end
%% Output all state variables

%AgIS, Therapeutic protein in the injection site, pmole
AgIS_vector=y_record(:, 1);
% Ag, Ag in the plasma compartment, pmole
Ag_vector=y_record(:,2);
% Agec, antigenic protein in the extra central compartment, pmole
Agec_vector=y_record(:,3);
% Agp, Ag in the peripheral compartment, pmole
Agp_vector=y_record(:,4);
% MS, Maturation signal for activation of immature dendritic cells
% (Here is LPS, ng)
MS_vector=y_record(:,5);
% ID,	Immature dendritic cells, number of cells
ID_vector=y_record(:,6);
% MD,	Mature dendritic cells, number of cells
MD_vector=y_record(:,7);
% cpE, endogenous competing protein in endosome, pmole
cpE_vector=y_record(:,8);
% cptE, endogenous competing peptide in endosome, pmole
cptE_vector=y_record(:,9);
% cptME, endogenous peptide-MHC complex in endosome, pmole
cptME_vector=y_record(:,(10:15));
% cptM, endogenous peptide-MHC complex on dendritic cell membrane, pmole
cptM_vector=y_record(:,(16:21));
%AgE, Ag in the endosomes, pmole
AgE_vector=y_record(:,22);
% pE,	free epitope peptides from Ag digestion , pmole
pE_vector=y_record(:,(22+1):(22+N));
% ME,	free MHC II molecule in endosome , pmole
ME_vector=y_record(:,(23+N):(28+N));
% pME,	T-epitope-MHC-II complex in endosome, pmole
pME_vector=y_record(:,(29+N):(28+7*N));
% pM,	T-epitope-MHC-II on dendritic cell membrane, pmole
pM_vector=y_record(:,(29+7*N):(28+13*N));
% M, free MHC II molecule on dendritic cell menbrane, pmole
M_vector=y_record(:,(29+13*N):(34+13*N));
% NT,	Na�ve helper T cells, number of cells
NT_vector=y_record(:,(35+13*N):(34+14*N));
% AT_N,	activated helper T cells derived from naive T cells, number of cells
AT_N_vector=y_record(:,(35+14*N):(34+15*N));
% AT_M,	activated helper T cells derived from memory T cells, number of cells
AT_M_vector=y_record(:,(35+15*N):(34+16*N));
% MT, memory helper T cells, number of cells
MT_vector=y_record(:,(35+16*N):(34+17*N));
% FT, functional helper T cells, number of cells
FT_vector=y_record(:,(35+17*N):(34+18*N));
% NB,	Na�ve B cells, number of cells
NB_vector(:,1:J)=y_record(:,(35+18*N):(34+18*N+J));
% AB_N,	Activated B cells derived form naive B cells, number of cells
AB_N_vector(:,1:J)=y_record(:, (35+18*N+J):(34+18*N+2*J));
% AB_M, Activated B cells derived from memory B cells, number of cells
AB_M_vector(:,1:J)=y_record(:, (35+18*N+2*J):(34+18*N+3*J));
% SP,	Short-lived plasma (antibody secreting) B cells, number of cells
SP_vector(:,1:J)=y_record(:, (35+18*N+3*J):(34+18*N+4*J));
% LP,	Long-lived plasma (antibody secreting) B cells, number of cells
LP_vector(:,1:J)=y_record(:, (35+18*N+4*J):(34+18*N+5*J));
% MB, Memory B cells, number of cells
MB_vector(:,1:J)=y_record(:, (35+18*N+5*J):(34+18*N+6*J));
% A: Total antibody, pmole
A_vector(:,1:J)= y_record(:, (35+18*N+6*J):(34+18*N+7*J));

%  A_total: The total amount of 17 clones of antibody, pmole
A_total=A_vector(:,1:J)*ones(J,1);

%% Post-processing calculation

% Average binding affinity of Ab during the experiment time
Ka_average=10.^(A_vector*log10(Ka)./A_total);


% Antigen presentation processes
pM_NUMBER_M1=pM_vector(:,(1:N))*1E-12*NA; % T-epitope-MHC II complex number associated with MHCII molecule 1
pM_NUMBER_M2=pM_vector(:,(N+1):(2*N))*1E-12*NA; % T-epitope-MHC II complex number associated with MHCII molecule 2
pM_NUMBER_M3=pM_vector(:,(2*N+1):(3*N))*1E-12*NA; % T-epitope-MHC II complex number associated with MHCII molecule 3
pM_NUMBER_M4=pM_vector(:,(3*N+1):(4*N))*1E-12*NA; % T-epitope-MHC II complex number associated with MHCII molecule 4
pM_NUMBER_M5=pM_vector(:,(4*N+1):(5*N))*1E-12*NA; % T-epitope-MHC II complex number associated with MHCII molecule 5
pM_NUMBER_M6=pM_vector(:,(5*N+1):(6*N))*1E-12*NA; % T-epitope-MHC II complex number associated with MHCII molecule 6

Total_pM(:,1:N)=pM_NUMBER_M1+pM_NUMBER_M2+pM_NUMBER_M3+pM_NUMBER_M4+pM_NUMBER_M5...
    +pM_NUMBER_M6;   % Total number of T-epitope-MHC II complex


%  Ag, Ab, and B cell receptor (BCR) binding equilibrium
Rows=numel(t_record);
% Asign intial empty vectors for the storage of varible values
R_vector=zeros(Rows,J);
ro_vector=zeros(Rows,J);
Agfree_vector= zeros(Rows,1);
Afree_vector=zeros(Rows,J);
Complex_vector=zeros(Rows,J);
BCR_vector=zeros(Rows,J);

for i=1:Rows
    Ag=Ag_vector(i);
    NB=NB_vector(i, :);
    AB_N=AB_N_vector(i,:);
    AB_M=AB_M_vector(i,:);
    MB=MB_vector(i,:);
    A=A_vector(i,:);
    
    % Calculate the free Ag concentration
    BCR(1:J)=(NB+AB_N+AB_M+MB)/NA*1E12*BRN; % BCR is the amount of BCR, pmole
    Eq1 = @(x)(Ag/Vp-x*(1+2*sum(Ka'.*A/Vp./(1+Ka'*x))+ sum(Ka'.*BCR/Vp./(1+Ka'*x))));
    Options = optimset('TolX',5e-16);
    if i==1
        Agfree_con = fzero(Eq1,Ag/Vp,Options); % the unit of Agfree is in concentration, pM.
    else
        Agfree_con = fzero(Eq1,xsolution2,Options);
    end
    xsolution2=Agfree_con;
    
    
    ro(1:J) = Ka*Agfree_con./ (1 + Ka*Agfree_con);% Calculate the B cell receptor occupancy "ro"
    R(1:J)=ro(1:J).*BRN; % Occupied BCR number
    R_vector(i,:)=R(1:J);
    ro_vector(i,:)=ro(1:J);
    
    Complex_con=ro(1:J).*(2*A/Vp);  % Antibody-Antigen complex concentration, pM
    BCR_con=ro(1:J).*(BCR/Vp);   % BCR concentration, pM
    Afree_con=(2*(A/Vp)-Complex_con)/2; % The scaler "2" converts the antibody binding site concentration to antibody concentration
    
    Agfree_vector(i,:)= Agfree_con*Vp;% convert concentration to amount, pmole
    Afree_vector(i,:)=Afree_con*Vp;  % convert concentration to amount, pmole
    Complex_vector(i,:)=Complex_con*Vp; % convert  concentration to amount, pmole
    BCR_vector(i,:)=BCR_con*Vp;% convert concentration to amount, pmole
end

%% Save the results
savefile='results.mat';

save(savefile, 'koff', 't_record','AgIS_vector','Ag_vector','Agec_vector','Agp_vector', 'MS_vector','ID_vector','MD_vector',...
    'cpE_vector', 'cptE_vector', 'cptME_vector', 'cptM_vector', 'AgE_vector','pE_vector', ...
    'ME_vector','pME_vector', 'pM_vector', 'M_vector', 'NT_vector', 'AT_N_vector', 'AT_M_vector', 'MT_vector',...
    'FT_vector', 'NB_vector','AB_N_vector','AB_M_vector','SP_vector','LP_vector','MB_vector','A_vector','A_total',...
    'Ka_average','Total_pM', 'Agfree_vector', 'Afree_vector', 'Complex_vector',...
    'BCR_vector','R_vector','ro_vector');
