% Performance Characteristics of a NACA 4 digit airfoil
% Coded by: Mir Owais Ali
% Concordia University
% Can be used for educational purposes only

clear
clc
close all

%% Inputs
TypeNACA = '1410';                                                         %Specify NACA 4-digit airfoil
cr = 1;                                                                    %Root chord (m)
nodes = 40;                                                                %Number of points along airfoil
alpha_min = -5;                                                            %Specify the minimum alpha for distribution (no more than 0)
alpha_max = 16;                                                            %Specify the maximum alpha for distribution (no less than 16)
n = 15;                                                                    %Number of points along the span
b = 8;                                                                     %Wingspan (m)
TR = 0.7;                                                                  %Taper Ratio

%% Initialization
z_c = zeros(nodes,1);
dz_dx = zeros(nodes,1);
theta = zeros(nodes,1);
t = zeros(nodes,1);
x_u = zeros(nodes,1);
z_u = zeros(nodes,1);
x_l = zeros(nodes,1);
z_l = zeros(nodes,1);
X = zeros(2*nodes-2,1);
Z = zeros(2*nodes-2,1);
S = zeros(2*nodes-2,1);
theta_n = zeros(2*nodes-2,1);
RHS = zeros(2*nodes-2,1);
CN1 = zeros(2*nodes-2,2*nodes-2);
CN2 = zeros(2*nodes-2,2*nodes-2);
CT1 = zeros(2*nodes-2,2*nodes-2);
CT2 = zeros(2*nodes-2,2*nodes-2);
AN = zeros(2*nodes-1,2*nodes-1);
AT = zeros(2*nodes-1,2*nodes-1);
V = zeros(1,2*nodes-2);
CP = zeros(1,2*nodes-2);
Cl_v = zeros(1,alpha_max-alpha_min+1);
Cm_v = zeros(1,alpha_max-alpha_min+1);
alpha_d = zeros(1,alpha_max-alpha_min+1);
Cl_th = zeros(1,alpha_max-alpha_min+1);
Cm_th = zeros(1,alpha_max-alpha_min+1);
CL_fw = zeros(1,alpha_max-alpha_min+1);
Cm_fw = zeros(1,alpha_max-alpha_min+1);
x_cp = zeros(1,alpha_max-alpha_min+1);
Coeff_A = zeros(n,2*n-1);
CD_i = zeros(1,alpha_max-alpha_min+1);
LHS = zeros(n,1);
Circ = zeros(n,1);

%% Airfoil Coordinates
c = cr;
m_c = str2double(TypeNACA(1))/100;                                         %Maximum Camber
p_c = str2double(TypeNACA(2))/10;                                          %Position of Maximum Camber
t_max = str2double(TypeNACA(3:4))/100;                                     %Maximum Thickness

psi = linspace(0,pi,nodes)';                                               %Distributing more nodes at the L.E. and T.E. 
x(:,1) = (0.5*(1-cos(psi)));

for i = 1:1:nodes
    if(x(i,1) >= 0 && x(i,1) < p_c)
        z_c(i,1) = (m_c/p_c^2)*((2*p_c*x(i,1))-x(i,1)^2);
        dz_dx(i,1) = (2*m_c/(p_c)^2)*(p_c-x(i,1));
    elseif(x(i,1) >= p_c && x(i,1) <= 1)
        z_c(i,1) = (m_c/(1-p_c)^2)*(1-(2*p_c)+(2*p_c*x(i,1))-(x(i,1)^2));
        dz_dx(i,1) = (2*m_c/(1-p_c)^2)*(p_c-x(i,1));
    end
theta(i,1) = atan(dz_dx(i,1));                                             %Airfoil Angle

t(i,1) = 5*t_max*(0.2969*sqrt(x(i,1))-0.126*x(i,1)-0.3516*x(i,1)^2....     %Thickness Distribution
                +0.2843*x(i,1)^3-0.1036*x(i,1)^4);

x_u(i,1) = x(i,1)-t(i,1)*sin(theta(i,1));                                  %Upper Surface
z_u(i,1) = z_c(i,1)+t(i,1)*cos(theta(i,1));

x_l(i,1) = x(i,1)+t(i,1)*sin(theta(i,1));                                  %Lower Surface
z_l(i,1) = z_c(i,1)-t(i,1)*cos(theta(i,1));

end
z_c = z_c*c;                                                               %Multilplying by chord to show the actual profile
x = x*c;
x_u = x_u*c;                                         
z_u = z_u*c;
x_l = x_l*c;                                          
z_l = z_l*c;

XB = [flip(x_l(2:nodes,1));x_u];                                           %Reassigning airfoil profile
ZB = [flip(z_l(2:nodes,1));z_u];

f1 = figure(1);
str = sprintf('NACA %s Airfoil Profile for chord =%d m',TypeNACA,c);
title(str);
axis equal
hold on
grid on
grid minor;
plot(x,z_c,'r-');                                                          %Mean Camber Line
plot([0,c],[0,0],'b-');                                                    %Chord Line
plot(XB,ZB,'k-o','LineWidth',1.5,'MarkerSize',2);                          %Airfoil Profile
hold off

%% Vortex Panel Method
XB = XB./c;                                                                %Reassigning airfoil profile
ZB = ZB./c;
M = length(XB);
z = 2;
k = 1;
for alpha = alpha_min:1:alpha_max                                                                           
    alpha = degtorad(alpha);                                               %Assigning A.O.A

                                                                           %Calculating Control Points
    for i=1:M-1                                                            % and the panel length
        X(i,1) = 0.5*(XB(i)+XB(i+1)); 
        Z(i,1) = 0.5*(ZB(i)+ZB(i+1));
        S(i,1) = sqrt((XB(i+1)-XB(i))^2+(ZB(i+1)-ZB(i))^2);
        theta_n(i,1) = atan2(ZB(i+1)-ZB(i),XB(i+1)-XB(i));
        RHS(i,1) = sin(theta_n(i)-alpha);
    end
                                                                           %Calculating coefficients
    for i=1:M-1
        for j=1:M-1
            if (i == j)
                CN1(i,j) = -1;
                CN2(i,j) = 1;
                CT1(i,j) = (1/2)*pi;
                CT2(i,j) = (1/2)*pi;
            else
                A = - (X(i) - XB(j))*(cos(theta_n(j))) - (Z(i) - ZB(j))*(sin(theta_n(j)));
                B = (X(i) - XB(j))^2 + (Z(i) - ZB(j))^2;
                C = sin(theta_n(i) - theta_n(j));
                D = cos(theta_n(i) - theta_n(j));
                E = (X(i) - XB(j))*sin(theta_n(j)) - (Z(i) - ZB(j))*cos(theta_n(j));
                F = log(1 + ((S(j))^2 + (2*A*S(j))) / B);
                G = atan2((E*S(j)) , (B + A*S(j)));
                P = ((X(i) - XB(j)) * sin(theta_n(i) - 2*theta_n(j))) + ((Z(i) - ZB(j)) * cos(theta_n(i) - 2*theta_n(j)));
                Q = ((X(i) - XB(j)) * cos(theta_n(i) - 2*theta_n(j))) - ((Z(i) - ZB(j)) * sin(theta_n(i) - 2*theta_n(j)));

                CN2(i,j) = D + ((0.5*Q*F)/S(j)) - ((A*C + D*E)*(G/S(j)));
                CN1(i,j) = 0.5*D*F + C*G - CN2(i,j);
                CT2(i,j) = C + ((0.5*P*F)/S(j)) + ((A*D - C*E)*(G/S(j)));
                CT1(i,j) = 0.5*C*F - D*G - CT2(i,j);
            end
        end
    end
                                                                           %Computation of Influence Coefficients
    for i = 1:M-1
        AN(i,1) = CN1(i,1);
        AN(i,M) = CN2(1,M-1);
        AT(i,1) = CT1(i,1);
        AT(i,M) = CT2(i,M-1);
        for j = 2:M-1
            AN(i,j) = CN1(i,j) + CN2(i,j-1);
            AT(i,j) = CT1(i,j) + CT2(i,j-1);
        end
    end
    AN(M,1) = 1;
    AN(M,M) = 1;
    for j = 2:M-1
        AN(M,j) = 0;
    end
    RHS(M) = 0;
                                                                           %Solving for a system of linear equations
    Gama = AN\RHS;                                                             

    for i = 1:M-1
        V(i) = cos(theta_n(i)-alpha);
        for j = 1:M
            V(i) = V(i) + AT(i,j)*Gama(j);
            CP(i) = 1 - (V(i))^2;
        end
    end
                                                                           %Plotting Cp vs x/c
    if (radtodeg(alpha) == 0 ||radtodeg(alpha) == 8 || radtodeg(alpha) == 16) 
        figure(z)
        plot(X,CP);
        hold on  
    end
    %Calculation of Lift Coefficient
    CPl = CP(1:((M-1)/2));
    CPl = flip(CPl);
    CPu = CP((((M-1)/2)+1):end);

    dCP = (CPl - CPu)';
    dx = X((((M-1)/2)+1):end);

    Cl_v(k) = trapz(dx,dCP);
    Cm_v(k) = trapz(dx,-dCP.*dx);
    alpha_d(k) = alpha;
    k = k+1;
end
    str = sprintf('C_p distribution (NACA %s)',TypeNACA);
    title(str);
    set(gca,'Ydir','reverse');
    xlabel('x/c');
    ylabel('Pressure Coefficient');
    legend('alpha = 0^0','alpha = 8^0','alpha = 16^0','Location','northwest');
    grid on;
    grid minor;
    hold off
z= z+1; 
k = 1;

%% Thin Airfoil Theory
e = acos(1-2*p_c);

th = e;
I1_u = (m_c/p_c^2)*(sin(th)+(2*p_c-1)*th);
I3_u = (m_c/(2*p_c^2))*((cos(th)+4*p_c-2)*sin(th)+th);
I5_u = (m_c/(3*p_c^2))*(cos(th)*(2*cos(th)+6*p_c-3)+1)*sin(th);

th = 0;
I1_l = (m_c/p_c^2)*(sin(th)+(2*p_c-1)*th);
I3_l = (m_c/(2*p_c^2))*((cos(th)+4*p_c-2)*sin(th)+th);
I5_l = (m_c/(3*p_c^2))*(cos(th)*(2*cos(th)+6*p_c-3)+1)*sin(th);

I1 = I1_u-I1_l;
I3 = I3_u-I3_l;
I5 = I5_u-I5_l;

th = pi;
I2_u = (m_c/(1-p_c)^2)*(sin(th)+(2*p_c-1)*th);
I4_u = (m_c/(2*(1-p_c)^2))*((cos(th)+4*p_c-2)*sin(th)+th);
I6_u = (m_c/(3*(1-p_c)^2))*(cos(th)*(2*cos(th)+6*p_c-3)+1)*sin(th);

th = e;
I2_l = (m_c/(1-p_c)^2)*(sin(th)+(2*p_c-1)*th);
I4_l = (m_c/(2*(1-p_c)^2))*((cos(th)+4*p_c-2)*sin(th)+th);
I6_l = (m_c/(3*(1-p_c)^2))*(cos(th)*(2*cos(th)+6*p_c-3)+1)*sin(th);

I2 = I2_u-I2_l;
I4 = I4_u-I4_l;
I6 = I6_u-I6_l;
for alpha = alpha_min:1:alpha_max
    alpha = degtorad(alpha);
    
    A0 = alpha - (1/pi)*(I1+I2);
    A1 = (2/pi)*(I3+I4);
    A2 = (2/pi)*(I5+I6);

    Cl_th(k) = pi*(2*A0+A1);
    Cm_th(k) = -(pi/2)*(A0+A1-A2/2);
    x_cp(k) = (c/4)*(1+(A1-A2)*pi/Cl_th(k));

    k = k+1;
end
alpha_L0 = alpha-Cl_th(k-1)/(2*pi);                                        %Zero-Lift Angle of Attack

%% Finite Wing Method
if (TR == 1)
    cr = c;
    ct = c;
elseif (TR < 1)
    cr;                                                                    %Specify root chord
    ct = TR*cr;
end
AR = (2*b)/(cr*(1+TR));                                                    %Aspect Ratio
S = b^2/AR;                                                                %Wing Area (m^2)
phi = linspace(1E-10,pi/2,n);                                              %Span angle distribution
c_d = ct+((cr-ct)/(pi/2))*phi;                                             %Chord Distribution
mue = c_d*2*pi/(4*b);
k = 1;
for alpha = alpha_min:1:alpha_max
    alpha = degtorad(alpha);
    for i = 1:1:n
        for j = 1:2:(2*n-1)
            Coeff_A(i,j) = sin(j*phi(i))*(j*mue(i)+sin(phi(i)));
        end
        LHS(i,1) = mue(i)*(alpha-alpha_L0)*(sin(phi(i)));
    end 
    A_fw_d = Coeff_A\LHS;
    delta = 0;
  
    for j = 3:2:(2*n-1)
        del = j*A_fw_d(j)^2/A_fw_d(1)^2;
        delta = delta+del;
        del = delta;
    end
    for i = 1:1:n
        for j = 1:2:(2*n-1) 
                Cir = A_fw_d(j)*sin(j*phi(i));
                Circ(i) = Circ(i)+Cir;
                Cir = Circ(i);
        end
    end
    CL_fw(k) = A_fw_d(1)*pi*AR;
    Cm_fw(k) = -CL_fw(k)*x_cp(k);
    CD_i(k) = (CL_fw(k)^2/(pi*AR))*(1+delta);
    k = k+1;
end
figure(z)                                                                  %Plotting Cd_i vs alpha
plot(radtodeg(alpha_d),CD_i);
str = sprintf('C_d_i vs alpha (NACA %s)',TypeNACA);
title(str);
xlabel('Angle of Attack (Degrees)');
ylabel('Induced Drag Coefficient');
grid on;
grid minor;
z = z+1;

figure(z)                                                                  %Plotting Circulation along the span
plot(0.5*(1-cos(phi)),Circ,'-b');
hold on
str = sprintf('Non-dimensional circulation distribution (NACA %s)',TypeNACA);
plot(0.5*(1-cos(pi+phi)),Circ,'-b')
title(str);
xlabel('y/b');
ylabel('Circulation coefficient');
grid on;
grid minor;
z = z+1;

%% Comparison
figure(z)                                                                  %Plotting Cl vs alpha
plot(radtodeg(alpha_d),Cl_v);
hold on
plot(radtodeg(alpha_d),Cl_th);
plot(radtodeg(alpha_d),CL_fw);
str = sprintf('C_l vs alpha (NACA %s)',TypeNACA);
title(str);
xlabel('Angle of Attack (Degrees)');
ylabel('Lift Coefficient');
legend('Vortex-Panel','Thin-Airfoil','Finite Wing','Location','northwest');
grid on;
grid minor;
z = z+1;

figure(z)                                                                  %Plotting Cm vs alpha
plot(radtodeg(alpha_d),Cm_v);
set(gca,'Ydir','reverse');
hold on
plot(radtodeg(alpha_d),Cm_th);
plot(radtodeg(alpha_d),Cm_fw);
str = sprintf('C_m vs alpha (NACA %s)',TypeNACA);
title(str);
xlabel('Angle of Attack (Degrees)');
ylabel('Moment Coefficient at L.E.');
legend('Vortex-Panel','Thin-Airfoil','Finite Wing');
grid on;
grid minor;
