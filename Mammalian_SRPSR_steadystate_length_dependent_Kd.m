

%% system parameters

n=150; % lenth of polypeptide chain
ss1=1; % start of signal sequence
ss2=1; % end of signal sequence
buried=35; % residues accomodated by ribosome


srp=0.5*10^-6; % [SRP]
sr=0.5*10^-6;     % [SR]


%% initializing state probability
rnc = zeros([4 n]);
rnc(1,1) = 1;

ksrpon = zeros([n 1]);
ksrpoff = zeros([n 1]);
ksron = zeros([n 1]);
ksroff = zeros([n 1]);
kd = zeros([n 1]);
%% SRP Kd length dependence
kd_base = 4 * 10^-9;


%% defining rate constants
        T = zeros([4 4]);
        A = zeros([4 4]);


for i = 1:n
    
        ksrpon(i)=1*10^7;
        ktrans=6;
        ktar=0.2;
    
        residue=min(ss2-ss1,max(0,i-buried-ss1));% kon koff khydro are all dependent on exposed hydrophobic sequence
        kd(i) = kd_base*(1+(14-residue)*6);      % Kd based on SRPRNC titration, assuming NAC bound(M)
        %kd(i) = kd_base;
        ksrpoff(i) = ksrpon(i) * kd(i);% parameters are from empirical fitting to experimental data
        ksron(i) = (7330*(1/(exp(residue*0.14-0.85755)+1))+1.28*10^6*(1-1/(exp(residue*0.14-0.85755)+1)));    % kon based on SRP-SRabdTM
        ksroff(i) = 100*10^-9*ksron(i);
        khydro = 0; % disabled for now

        
        D = [ktrans+ksrpon(i)*srp ktrans+ksrpoff(i)+ksron(i)*sr ktrans+ksroff(i)+ktar ktrans];
        
        T(:,:,i)=[-ksrpon(i)*srp ksrpoff(i) 0 0;
                  ksrpon(i)*srp -ksrpoff(i)-ksron(i)*sr ksroff(i) 0;
                  0 ksron(i)*sr -ksroff(i)-ktar 0;
                  0 0 ktar 0];
              
        translation = [ktrans ktrans ktrans ktrans];      
        A(:,:,i)=diag(translation);
        
       
      
    
end


%% solving steady state probability

for i = 1:n-1
    
    rnc(:,i+1) = (A(:,:,i)-T(:,:,i))\A(:,:,i)*rnc(:,i);   
    
    
    
    
    
    
end
