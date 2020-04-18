
function [target] = Mammalian_SRPSR_2R_NAC_steadystate_function(SRP,SR,elongation,ktarget,KdSRP,konSRP,KdSR)
%% system parameters

n=150; % lenth of polypeptide chain
ss1=1; % start of signal sequence
ss2=1; % end of signal sequence (1 for ssmt, 15 for ss)
buried=35; % residues accomodated by ribosome


srp=SRP;      % [SRP]
sr=SR;     % [SR]


%% initializing state probability
rnc = zeros([4 n]);
rnc(1,1) = 1;

ksrpon = zeros([n 1]);
ksrpoff = zeros([n 1]);
ksron = zeros([n 1]);
ksroff = zeros([n 1]);



%% defining rate constants
        T = zeros([4 4]);
        A = zeros([4 4]);


for i = 1:n
    
        ksrpon(i)=konSRP;
        ktrans=elongation;
        ktar=ktarget;
    
        residue=min(ss2-ss1,max(0,i-buried-ss1));% kon koff khydro are all dependent on exposed hydrophobic sequence
        kd = KdSRP;      % Kd based on SRPRNC titration, assuming NAC bound(M)
        ksrpoff(i) = ksrpon(i) * kd;% parameters are from empirical fitting to experimental data
        ksron(i) = 13400;    % kon based on SRP-SRabdTM
        ksroff(i) = KdSR*ksron(i);
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
target = rnc(4,:);

end
