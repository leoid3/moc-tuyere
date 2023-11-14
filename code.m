clc
cla

%Hypothèses: 
%-Isentropique
%-Tuyère adaptée et amorcée au niveau de la mer
%-Gaz coloriquement parfait
%-Pas de choc en sortie de tuyère

% Données d'entrée
%Pc = input('Pression de la chambre (en Pa) : ');
%Tc = input('Température de la chambre (en K) : ');
%F = input('Poussée désirée (en N) : ');
%alt = input('Altitude (en m) : ');
Pc=100000;
Tc=2000;
F=2000;
alt=0;
gamma = 1.4;
R = 8.314;
T0=288;
P0=101325;
g=9.81;

if alt<=0
    Pe= P0;
elseif (alt>0) && (alt<=11000)
    Tz= T0 + (-0.0065)*alt;
    Pe= P0 *((1+(Tz/T0)*alt))^(-5.2558);
elseif (alt>11000) && (alt<=20000)
  return
end

%Calculs preliminaires
Pe=101325;
RP=Pe/Pc;
Me=sqrt((2/(gamma-1))*((RP)^((gamma-1)/(-gamma))-1)); %Mach de sortie
RS=(1/Me)*((2/(gamma+1))*(1+((gamma-1)/2)*(Me)^2))^((gamma+1)/(2*(gamma-1)));%Rapport de section
phi1=1.268;
phi2=RP*RS*(1+gamma*(Me)^2); %Dynalpie a la sortie
Ac=F/((phi2-phi1)*Pc); %Section au col
Rc=sqrt(Ac/(2*pi)); %Rayon au col
Ae=Ac*RS; %Section a la sortie
Re=sqrt(Ae/(2*pi)); %Rayon a la sortie
Me=2.4;
theta_max=rad2deg(sqrt((gamma+1)/(gamma-1))*atan(sqrt((gamma-1)/(gamma+1)*(Me^2-1)))-atan(sqrt(Me^2-1)))/2; %Angle max de Prandlt-Meyer

%Methode des caracteristiques
nc=4; %Nombre de characteristique
theta_0=theta_max/nc; %angle à l'origine
dtheta=(theta_max-theta_0)/(nc-1); %Pas d'angle
node=0.5*nc*(4+nc-1); %Nombre de noeuds

theta=zeros(1, node);
v=zeros(1, node);
Km=zeros(1, node); %caracteristique descendente (C_)
Kp=zeros(1, node); %caracteristique montante (C+)

%Calcule les C_,C+ et les angles de Prandtl-Meyer pour les nc premiers noeuds (ceux sur la première caractéristique)
for i=1:nc
    theta(i)=theta_0 + (i-1)*dtheta;
    v(i)=theta(i);
    Km(i)=theta(i)-v(i);
    Kp(i)=theta(i)+v(i);
end 

%Calcul du premier noeud sur la paroi
j=nc+1; %Variable pour les noeuds à la paroi
theta(j)=theta(nc);
v(j)=v(nc);
Kp(j)=Kp(nc);
Km(j)=Km(nc);
q=nc+2; %Variable pour les noeuds sur l'axe de symetrie/de revolution de la tuyere
m=q+1;%Variable pour les noeuds entre la paroi et l'axe de symetrie/de revolution de la tuyere
j=j+(nc);
k=0;

%Calcule les C_,C+ et les angles de Prandtl-Meyer pour les noeuds restants
%une caracteristique à la fois
for i=q:node
    if q<node %noeuds sur l'axe de révolution
        theta(q)=0;
        Kp(q)=Kp(q-(nc-k));
        v(q)=Kp(q)-theta(q);
        Km(q)=theta(q)-v(q);      
        for l=m: j -1 %noeuds entre la paroi et l'axe de symetrie/de revolution de la tuyere    
            Kp(l)=Kp(l-(nc-k));
            Km(l)=Km(l-1);
            theta(l)=(Km(l)+Kp(l))/2;
            v(l)=(Kp(l)-Km(l))/2;
            h=l;
        end
    end
    if j<=node %noeuds sur la paroi
        theta(j)=theta(h);
        v(j)=v(h);
        Kp(j)=Kp(h);
        Km(j)=Km(h);
    else
        break
    end
    q=q+(nc-k);
    m=q+1;
    if m==node
        h=q;
    end
    j=j+(nc-k-1);
    k=k+1;
end

%Calcul des coordonées (x,y) des noeuds
