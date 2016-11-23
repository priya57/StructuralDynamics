
% Moment of Inertia for Design (A)
E = 210 * 10^9;  % N/m^2
rho = 7850;  % kg/m^3
I  = 9.4751 * 10^(7) * 10^(-12); % m^4

% Moment of Inertia for Design (B)
I_b(1) = 73671200 * 10^(-12); % m^4
I_b(2) = 4.3778 * 10^7* 10^(-12) ;
I_b(3) =  2.3990e+07 * 10^(-12) ;

% Original Area of the beams
A = 6014; % mm^2

% Disceretize areas
A1 = 5400 ;
A2 = 4420;
A3 = 3520;

% Mass of the I Beams in Kg/m
ma = 47.21 ; % ^kg/m

% Discretize masses
m(1) = 42.39;
m(2) = 34.697;
m(3) = 27.632;
% Mass of the Concrete Slab + Fittings

Mcs = 500 ; % ^kg/m

% Mass of  the People
Mp = 593.75 ; % ^ kg/m

% Mass of the Design (a)

M_a = ma+Mcs+Mp;

% Mass of the Design B

for i = 1:3
    M_b(i) = m(i)+(Mcs+Mp);
end

% Mass Matrix 

l= 1.6;
L = 4.8;
M_Facta = M_a*l/420;
K_Facta = (E*I)/l^3;

E = 210 * 10^9;
for i = 1:3
    M_Factb(i) =  M_b(i)*l/420;
    K_Factb(i) = (E*I_b(i))/l^3;
end



% Global Mass Matrix for Design (a)
prompt='Enter No of nodes from 1 to 4';
n=input(prompt)

M=zeros(4,4);

M     =     [156          22*l        54          -13*l  ;  ...,
             22*l         4*l*l       13*l        -3*l*l ;  ...,
             54           13*l        156         -22*l  ;  ...,       
            -13*l        -3*l*l      -22*l         4*l*l ]  ;
         
K    =      [12            6*l       -12         6*l    ;  ...,
             6*l           4*l*l     -6*l        2*l*l  ;  ...,
            -12           -6*l        12        -6*l    ;  ...,
             6*l           2*l*l     -6*l        4*l*l  ];
         
   Ma=zeros(8,8);
   Mb=zeros(8,8);
   
   Ka=zeros(8,8);
   Kb=zeros(8,8);
   
   M1=zeros(8,8);
   M2=zeros(8,8); 
   M3=zeros(8,8);
   
   K1=zeros(8,8);
   K2=zeros(8,8); 
   K3=zeros(8,8); 
   
   
   M1(1:4,1:4)=M;
   M2(3:6,3:6)=M;
   M3(5:8,5:8)=M;
      
   K1(1:4,1:4)=K;
   K2(3:6,3:6)=K;
   K3(5:8,5:8)=K;
       
         i=n;
            if(i==2) 
               
                Ma  = M_Facta*M1+Ma;
                Mb  = M_Factb(i-1)*M1+Mb;
                
                Ka  = K_Facta*K1+Ka;
                Kb  = K_Factb(i-1)*K1+Kb;
                
                MGa = Ma(1:4,1:4);
                MGb = Mb(1:4,1:4);
                KGa = Ka(1:4,1:4);
                KGb = Kb(1:4,1:4);

            elseif(i==3)
                
                Ma  = Ma+M_Facta*M1+M_Facta*M2;
                Mb  = Mb+M_Factb(i-2)*M1+M_Factb(i-1)*M2;
                
                Ka  = Ka+K_Facta*K1+K_Facta*K2;
                Kb  = Kb+K_Factb(i-2)*K1+K_Factb(i-1)*K2;
                
                
                MGa = Ma(3:6,3:6);
                MGb = Mb(3:6,3:6);
                
                    
                KGa = Ka(3:6,3:6);
                KGb = Kb(3:6,3:6);
                
            elseif(i==4)
                
                 Ma  = Ma+M_Facta*M1+M_Facta*M2+M_Facta*M3;
                 Mb  = Mb+M_Factb(i-3)*M1+M_Factb(i-2)*M2+M_Factb(i-1)*M3;
                 
                 Ka  = Ka+K_Facta*K1+K_Facta*K2+K_Facta*K3;
                 Kb  = Kb+K_Factb(i-3)*K1+K_Factb(i-2)*K2+K_Factb(i-1)*K3;
                 
                 MGa = Ma(3:8,3:8);
                 MGb = Mb(3:8,3:8);
                 
                 KGa = Ka(3:8,3:8);
                 KGb = Kb(3:8,3:8);
                 

            end
             % Calculating Eigen values and Eigen Vectors
 
%Eigenvectors and eigenvalues:
[Xa,lambdaA]=eig(MGa\KGa,'vector');
[Xb,lambdaB]=eig(MGb\KGb,'vector');


%Natural circular frequencies:
omegaA=sqrt(lambdaA);
omegaB=sqrt(lambdaB);

%  Natural frequencies:
fnA=omegaA/(2*pi);
fnB=omegaB/(2*pi);

FnA = real(fnA); 
FnB = real(fnB); 

% Normalizing the Eigen Vectors for both design (A) and (B)

XAnew = Xa;
XBnew = Xb;

for i = 1:6
    XAnew(:,i) = Xa(:,i)/min(Xa(:,i));
    
    XBnew(:,i) = Xb(:,i)/min(Xb(:,i));
end                 
 
% Length of the discretised beam

X = [0 1.6 3.2 4.8];

% Time Steps by using NEWMARK Time Discretization Method
for design = 1:2
    
for k = 1:6
ti = 0. ;
if design == 1
 tf = 4*2*pi/omegaA(k) ;
else 
   tf = 4*2*pi/omegaB(k) ; 
end
t = linspace(ti,tf) ;
dt = t(2)-t(1);

nstep = length(t) ;
n = length(MGa) ;

U = zeros(n,nstep) ;
V = zeros(n,nstep) ;
A = zeros(n,nstep) ;

% Initial Conditions
A(:,1) = zeros ;
V(:,1) = zeros ;

% Which Eigen vector should be used ? 
% Normalized eigen vector or normal eigen vector 

if design == 1
    U(:,1) = Xa(:,k);
else 
    U(:,1) = Xb(:,k);
end


for i = 1:nstep-1
    
    if design == 1
        
        Numerator   = 4/dt^2 * MGa * U(:,i) + 4/dt * MGa * V(:,i) + MGa * A(:,i);
        Denominator = 4/dt^2 * MGa + KGa;
    else 
        Numerator   = 4/dt^2 * MGb * U(:,i) + 4/dt * MGb * V(:,i) + MGb * A(:,i);
        Denominator = 4/dt^2 * MGa + KGb;
    end
    
    U(:,i+1) = (Denominator)\(Numerator) ;
    
    V(:,i+1) = 2/dt*(U(:,i+1)-U(:,i))-V(:,i);
    
    A(:,i+1) = 4/dt^2*(U(:,i+1)-U(:,i)) - 4/dt*V(:,i) - A(:,i);
    
end

u1=zeros(1,nstep);
Disp =[u1; U(1,:); U(3,:); U(5,:)];

if design == 1
    figure(k)
    title(['Mode ' k])
    plot(t,Disp)  
    hold on
    xlabel('Time in sec') ;ylabel('U') ;
    legend('Node1','Node2','Node3','Node4')
else 
    figure(k+9)
    title(['Mode for deisgn B ' k])
    plot(t,Disp)  
    hold on
    xlabel('Time in sec') ;ylabel('U') ;
    legend('Node1','Node2','Node3','Node4')
end

end
end






        
          

