%% In Plane Stress Calculation in Composite Laminates for Classical Lamination Theory
% Orthotropic relations hold through the layers since 
% Classical lamination Theory is based upon that.
% Plane stress is assumed since Classical lamination theory is based upon
% that.

%Read the data from txt file and assign them into arrays
file='CompositeProject.txt';
X=readtable(file);
M=table2array(X);

 for i=1:numel(M(:,1))
l(i,1)=M(i,1);  % layer number from bottom to top,
T(i,1)=M(i,2);  % layer orientations from bottom to top
t(i,1)=M(i,3);  % layer thicknesses from bottom to top
C(i,1)=M(i,4);  % elastic coefficients of the material in order:E1,E2,E3,G13,G23,G12,v12,v23,v21
R(i,1)=M(i,5);  % Resultant Force-Moment Vector in order:Nx,Ny,Nxy,Mx,My,Mxy
 end

% To clear the NaN terms from the arrays and getting the arrays 
% in correct size is important for matrix operations 
 l= l';
 l = l(~isnan(l))';
 T = T';
 T = T(~isnan(T))';
 t = t';
 t = t(~isnan(t))';
 C = C';
 C = C(~isnan(C))';
 R = R';
 R = R(~isnan(R))';
%totalt is laminate thickness and nl is Number of Ply
%totalt2 is half of the Laminate thickness
nl=max(l);
totalt=0;
for i=1:nl
totalt=t(i,1)+totalt;
end
totalt2=totalt/2;
h(1)=-totalt2;
% By this formulation, code is able to solve for different ply thicknesses 
for i=1:nl
h(i+1)=h(i)+t(i);
end
% Request an input in meters! to user to specify the z coordinate of interest 
prompt = 'Where you want to find the in plane stresses? (in meters) ';
z = input(prompt);

E1=C(1);E2=C(2);E3=C(3);G13=C(4);G23=C(5);
G12=C(6);v12=C(7);v23=C(8);v21=C(9);
%% Q Matrix is obtained 
% 3 indice accunts for 6 indice of the contracted notation

Q(1,1)=E1/(1-v12*v21);
Q(2,2)=E2/(1-v12*v21);
Q(1,2)=v12*E2/(1-v12*v21);Q(2,1)=Q(1,2);
Q(3,3)=G12;
Q(1,3)=0;Q(3,1)=0;
Q(2,3)=0;Q(3,2)=0;
% O(4,4) and O(5,5) are simply Q44 and Q55 elements.
O(4,4)=G23;
O(5,5)=G13;

for i=1:3
    for j=1:3
        A(i,j)=0.0000000;
        B(i,j)=0.0000000;
        D(i,j)=0.0000000;
    end
end

%Calculate A,B and D matrices
%For each layer transformations are being applied with help of T matrix
for i=1:3
    for j=1:3
     for k=1:nl
         m=cosd(T(k));n=sind(T(k));
Qbar(1,1)=m^4*Q(1,1)+2*(Q(1,2)+2*Q(3,3))*m^2*n^2+Q(2,2)*n^4;
Qbar(2,2)=n^4*Q(1,1)+2*(Q(1,2)+2*Q(3,3))*m^2*n^2+Q(2,2)*m^4;
Qbar(1,2)=m^4*Q(1,2)+Q(1,2)*n^2+(Q(1,1)+Q(2,2)-4*Q(3,3))*m^2*n^2;
Qbar(3,3)=Q(3,3)*(m^2-n^2)^2+(Q(1,1)+Q(2,2)-2*Q(1,2))*m^2*n^2;
Qbar(1,3)=-m*n^3*Q(2,2)+m^3*n*Q(1,1)-m*n*(m^2-n^2)*(Q(1,2)+2*Q(3,3));
Qbar(2,3)=-n*m^3*Q(2,2)+n^3*m*Q(1,1)+m*n*(m^2-n^2)*(Q(1,2)+2*Q(3,3));
Qbar(3,2)=Qbar(2,3);Qbar(3,1)=Qbar(1,3);Qbar(2,1)=Qbar(1,2);
                  
sumA=Qbar(i,j)*(h(k+1)-h(k));
sumB=Qbar(i,j)*(h(k+1)^2-h(k)^2)/2;
sumD=Qbar(i,j)*(h(k+1)^3-h(k)^3)/3;

A(i,j)=A(i,j)+sumA;
B(i,j)=B(i,j)+sumB;
D(i,j)=D(i,j)+sumD;
  
     end
    end
end

% Calculation of Transverse shear stifness elements A44, A45, A55
% In the code, they are symbolized with P(4,4),P(4,5) and P(5,5)
P(4,5)=0;P(4,4)=0;P(5,5)=0;
 for k=1:nl
      m=cosd(T(k));n=sind(T(k));
Obar(4,4)=O(4,4)*m^2+O(5,5)*n^2;
Obar(4,5)=(O(5,5)-O(4,4))*m*n
Obar(5,5)=O(5,5)*m^2+O(4,4)*n^2;
     sum45=1.25*Obar(4,5)*((h(k+1)-h(k)-(4/3)*(h(k+1)^3-h(k)^3)/(h(5)-h(1))^2));
    sum44=1.25*Obar(4,4)*((h(k+1)-h(k)-(4/3)*(h(k+1)^3-h(k)^3)/(h(5)-h(1))^2)); 
 sum55=1.25*Obar(5,5)*((h(k+1)-h(k)-(4/3)*(h(k+1)^3-h(k)^3)/(h(5)-h(1))^2));
 
 P(4,5)=P(4,5)+sum45;
 P(4,4)=P(4,4)+sum44;
 P(5,5)=P(5,5)+sum55;
 end


%% Create K matrix (Laminate stiffness Matrix) consists of A,B and D matrices

for i=1:3
    for j=1:3
        K(i,j)=A(i,j);
    end
end
for i=1:3
    for j=4:6
        K(i,j)=B(i,j-3);
    end
end
for i=4:6
    for j=1:3
        K(i,j)=B(i-3,j);
    end
end
for i=4:6
    for j=4:6
        K(i,j)=D(i-3,j-3);
    end
end


%% Calculation of the in plane layer stresses

MPdef=inv(K)*R;

%Define mid plane strains and curvatures
for i=1:3
Eps(i)=MPdef(i);
    Kap(i)=MPdef(i+3);
end


for k=1:nl
    if(h(k)<z<h(k+1))
    m=cosd(T(k));n=sind(T(k));
Qbar(1,1)=m^4*Q(1,1)+2*(Q(1,2)+2*Q(3,3))*m^2*n^2+Q(2,2)*n^4;
Qbar(2,2)=n^4*Q(1,1)+2*(Q(1,2)+2*Q(3,3))*m^2*n^2+Q(2,2)*m^4;
Qbar(1,2)=m^4*Q(1,2)+Q(1,2)*n^2+(Q(1,1)+Q(2,2)-4*Q(3,3))*m^2*n^2;
Qbar(3,3)=Q(3,3)*(m^2-n^2)^2+(Q(1,1)+Q(2,2)-2*Q(1,2))*m^2*n^2;
Qbar(1,3)=-m*n^3*Q(2,2)+m^3*n*Q(1,1)-m*n*(m^2-n^2)*(Q(1,2)+2*Q(3,3));
Qbar(2,3)=-n*m^3*Q(2,2)+n^3*m*Q(1,1)+m*n*(m^2-n^2)*(Q(1,2)+2*Q(3,3));
Qbar(3,2)=Qbar(2,3);Qbar(3,1)=Qbar(1,3);Qbar(2,1)=Qbar(1,2);

    end
end

sig=Qbar*Eps.'+z*Qbar*Kap.';
%R;
%disp(inv(K));
%'K';,disp(K(:,:));
%'KAP',disp(Kap(:,:))
%'EPS',disp(Eps(:,:))

'[A]=',disp(A(:,:))
'[B]=',disp(B(:,:))
'[D]=',disp(D(:,:))
%'Transverse Shear Stiffness Coefficients'
%'A44=',disp(P(4,4))
%'A45=',disp(P(4,5))
%'A55=',disp(P(5,5))

%'In plane stress vector is',disp(sig(:,:))
% disp(Qbar(:,:))
% disp(Q(:,:))

