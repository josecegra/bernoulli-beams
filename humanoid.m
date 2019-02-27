clear
%number of nodes
n = 20;
%length of the beams
L = 1;
%lattice spacing
h = L/(n-1);
%lattice
xh = linspace(0,L,n)';

%mass density
mu = 1;
%vertical shear force
Qv = 0.5;
%shear momentum
Mv = 0.5;
%horizontal shear force
Qh = -0;
%beam in which the force is applied
bf = 5;

%final time
T = 20;
%number of time steps
nt = 200;
%time step
ht = T/nt;

%geometry of the network
%number of junctions
N_coupl = 11;
%number of beams
N = 12;
%angle of the first beam wrt x axis
alpha0 = degtorad(30);
%coordinates of the nodes
x_pos = zeros(2,12);
x_pos(1,1) = 0;
x_pos(2,1) = 0;
x_pos(1,2) = L*cos(alpha0);
x_pos(2,2) = L*sin(alpha0);
x_pos(1,3) = L*cos(alpha0);
x_pos(2,3) = L*sin(alpha0);
x_pos(1,4) = x_pos(1,3);
x_pos(2,4) = x_pos(2,3)+L;
x_pos(1,5) = x_pos(1,4)+L;
x_pos(2,5) = x_pos(2,4);
x_pos(1,6) = x_pos(1,4);
x_pos(2,6) = x_pos(2,4);
x_pos(1,7) = x_pos(1,6)-L;
x_pos(2,7) = x_pos(2,6);
x_pos(1,8) = x_pos(1,4);
x_pos(2,8) = x_pos(2,4);
x_pos(1,9) = x_pos(1,8);
x_pos(2,9) = x_pos(2,8)+L;
x_pos(1,10) = x_pos(1,9);
x_pos(2,10) = x_pos(2,9);
x_pos(1,11) = x_pos(1,9)-L*cos(degtorad(45));
x_pos(2,11) = x_pos(2,9)+L*sin(degtorad(45));
x_pos(1,12) = x_pos(1,10)+L*cos(degtorad(45));
x_pos(2,12) = x_pos(2,10)+L*sin(degtorad(45));

%angles
alpha = zeros(12,1);
alpha(1) = alpha0;
alpha(2) = -alpha0;
alpha(3) = degtorad(90);
alpha(4) = 0;
alpha(5) = alpha0;
alpha(6) = pi;
alpha(7) = pi-alpha0;
alpha(8) = degtorad(90);
alpha(9) = degtorad(135);
alpha(10) = pi/4;
alpha(11) = pi/4;
alpha(12) = 3*pi/4;

%junctions (i,j) for beams i and j
coupl = zeros(11,2);
coupl(1,1) = 1;  coupl(1,2) = 3; 
coupl(2,1) = 1;  coupl(2,2) = 2; 
coupl(3,1) = 3;  coupl(3,2) = 4; 
coupl(4,1) = 4;  coupl(4,2) = 5; 
coupl(5,1) = 3;  coupl(5,2) = 6; 
coupl(6,1) = 6;  coupl(6,2) = 7;
coupl(7,1) = 3;  coupl(7,2) = 8; 
coupl(8,1) = 8;  coupl(8,2) = 9; 
coupl(9,1) = 8;  coupl(9,2) = 10; 
coupl(10,1) = 9;  coupl(10,2) = 11; 
coupl(11,1) = 10;  coupl(11,2) = 12; 


%declaration of the 3rd degree polynomials
syms x
syms phi1(x)
syms phi2(x)
syms phi3(x)
syms phi4(x)
phi1(x) = 1-3*x^2+2*x^3;
phi2(x) = x*(x-1)^2;
phi3(x) = 3*x^2-2*x^3;
phi4(x) = (x^2)*(x-1);
%Assembling of the basis functions; 
%they are stored in a (2*n,2) matrix. 
%
%The row index i represent the
%index of the basis function (i = 1,...,2*n)
%
%The column index j represent the domain of definition
%of each function: 1 [xi-1,xi]
%                  2 [xi,xi+1]
phi=zeros(2*n,2);
phi=sym(phi);
phi(1,1) = 0;
phi(2,1) = 0;
phi(1,2) = phi1(x/h);
phi(2,2) = h*phi2(x/h);
for i = 2:(n-1)        
    phi(2*i-1,1) = phi3((x-xh(i-1))/h);
    phi(2*i-1,2) = phi1((x-xh(i))/h);
    phi(2*i,1) = h*phi4((x-xh(i-1))/h);
    phi(2*i,2) = h*phi2((x-xh(i))/h);    
end
phi(2*n-1,1) = phi3((x-xh(n-1))/h);
phi(2*n,1) = h*phi4((x-xh(n-1))/h);
phi(2*n-1,2) = 0;
phi(2*n,2) = 0;
phi = simplify(phi);
%Second derivative of the basis functions
ddphi = simplify(diff(phi,x,2));

%Stiffness matrix
%Only the non zero elements will be computed, being those corresponding
%to the overlap between the basis functions.
%only the upper triangular part will be computed
S = zeros(2*n,2*n);
S(1,1) = vpa(int(ddphi(1,2)*ddphi(1,2),x,0,h));
S(2,2) = vpa(int(ddphi(2,2)*ddphi(2,2),x,0,h));
S(1,2) = vpa(int(ddphi(1,2)*ddphi(2,2),x,0,h));
S(1,3) = vpa(int(ddphi(1,2)*ddphi(3,1),x,0,h));
S(1,4) = vpa(int(ddphi(1,2)*ddphi(4,1),x,0,h));
S(2,3) = vpa(int(ddphi(2,2)*ddphi(3,1),x,0,h));
S(2,4) = vpa(int(ddphi(2,2)*ddphi(4,1),x,0,h));
for i = 2:n-1      
    S(2*i-1,2*i-1) = vpa(int(ddphi(2*i-1,1)*ddphi(2*i-1,1),x,xh(i-1),xh(i)))+...
        vpa(int(ddphi(2*i-1,2)*ddphi(2*i-1,2),x,xh(i),xh(i+1)));    
    S(2*i-1,2*i) = vpa(int(ddphi(2*i-1,1)*ddphi(2*i,1),x,xh(i-1),xh(i)))+...
        vpa(int(ddphi(2*i-1,2)*ddphi(2*i,2),x,xh(i),xh(i+1)));   
    S(2*i-1,2*i+1) = vpa(int(ddphi(2*i-1,2)*ddphi(2*i+1,1),x,xh(i),xh(i+1)));    
    S(2*i-1,2*i+2) = vpa(int(ddphi(2*i-1,2)*ddphi(2*i+2,1),x,xh(i),xh(i+1)));    
    S(2*i,2*i) = vpa(int(ddphi(2*i,1)*ddphi(2*i,1),x,xh(i-1),xh(i)))+...
        vpa(int(ddphi(2*i,2)*ddphi(2*i,2),x,xh(i),xh(i+1)));    
    S(2*i,2*i+1) = vpa(int(ddphi(2*i,2)*ddphi(2*i+1,1),x,xh(i),xh(i+1)));    
    S(2*i,2*i+2) = vpa(int(ddphi(2*i,2)*ddphi(2*i+2,1),x,xh(i),xh(i+1)));
end
S(2*n-1,2*n-1) = vpa(int(ddphi(2*n-1,1)*ddphi(2*n-1,1),x,xh(n-1),L));
S(2*n-1,2*n) = vpa(int(ddphi(2*n-1,1)*ddphi(2*n,1),x,xh(n-1),L));
S(2*n,2*n) = vpa(int(ddphi(2*n,1)*ddphi(2*n,1),x,xh(n-1),L));
%building the stiffness matrix from the upper triangular part
S = S + S';
for i = 1:size(S)
    S(i,i) = S(i,i)/2;
end
Sw1 = S;

%mass matrix
%only the upper triangular part will be computed
M = zeros(2*n,2*n);
M(1,1) = vpa(int(phi(1,2)*phi(1,2),x,0,h));
M(2,2) = vpa(int(phi(2,2)*phi(2,2),x,0,h));
M(1,2) = vpa(int(phi(1,2)*phi(2,2),x,0,h));
M(1,3) = vpa(int(phi(1,2)*phi(3,1),x,0,h));
M(1,4) = vpa(int(phi(1,2)*phi(4,1),x,0,h));
M(2,3) = vpa(int(phi(2,2)*phi(3,1),x,0,h));
M(2,4) = vpa(int(phi(2,2)*phi(4,1),x,0,h));
for i = 2:n-1      
    M(2*i-1,2*i-1) = vpa(int(phi(2*i-1,1)*phi(2*i-1,1),x,xh(i-1),xh(i)))+...
        vpa(int(phi(2*i-1,2)*phi(2*i-1,2),x,xh(i),xh(i+1)));    
    M(2*i-1,2*i) = vpa(int(phi(2*i-1,1)*phi(2*i,1),x,xh(i-1),xh(i)))+...
        vpa(int(phi(2*i-1,2)*phi(2*i,2),x,xh(i),xh(i+1)));   
    M(2*i-1,2*i+1) = vpa(int(phi(2*i-1,2)*phi(2*i+1,1),x,xh(i),xh(i+1)));    
    M(2*i-1,2*i+2) = vpa(int(phi(2*i-1,2)*phi(2*i+2,1),x,xh(i),xh(i+1)));    
    M(2*i,2*i) = vpa(int(phi(2*i,1)*phi(2*i,1),x,xh(i-1),xh(i)))+...
        vpa(int(phi(2*i,2)*phi(2*i,2),x,xh(i),xh(i+1)));    
    M(2*i,2*i+1) = vpa(int(phi(2*i,2)*phi(2*i+1,1),x,xh(i),xh(i+1)));    
    M(2*i,2*i+2) = vpa(int(phi(2*i,2)*phi(2*i+2,1),x,xh(i),xh(i+1)));
end
M(2*n-1,2*n-1) = vpa(int(phi(2*n-1,1)*phi(2*n-1,1),x,xh(n-1),L));
M(2*n-1,2*n) = vpa(int(phi(2*n-1,1)*phi(2*n,1),x,xh(n-1),L));
M(2*n,2*n) = vpa(int(phi(2*n,1)*phi(2*n,1),x,xh(n-1),L));
%building the mass matrix from the upper triangular part
M = M + M';
for i = 1:size(M)
    M(i,i) = M(i,i)/2;
end
%multiplying by the mass density (it should be included inside the integrals)
M = mu*M;
Mw1 = M;
%boundary conditions
eL = zeros(2*n,1);
eL = sym(eL);
eL(2*n-3) = phi(2*n-3,2); eL(2*n-2) = phi(2*n-2,2);
eL(2*n-1) = phi(2*n-1,1); eL(2*n) = phi(2*n,1);
dL = diff(eL,x);
dLf = matlabFunction(dL);
dL = dLf(L);
eLf = matlabFunction(eL);
eL = eLf(L);
%right hand side
fw1 = zeros(2*n,1);

%Horizontal direction
%beam parameters
syms x
syms E 
syms A
%Young's module
E = 1;
%area
A = 1;
%basis functions
phi = zeros(n,2);
phi = sym(phi);
phi(1,1) = 0; phi(1,2) = (h-x)/h;
for i = 2:n-1
    phi(i,1) = (x-xh(i-1))/h;
    phi(i,2) = (xh(i+1)-x)/h;
end
phi(n,1) = (x-(L-h))/h; phi(n,2) = (L+h-x)/h;
dphi = diff(phi,x,1);

%stiffness matrix
S = zeros(n,n);
S(1,1) = vpa(int(E*A*dphi(1,2)*dphi(1,2),x,xh(1),xh(2)));
S(1,2) = vpa(int(E*A*dphi(1,2)*dphi(2,1),x,xh(1),xh(2)));
for i = 2:n-1 
    S(i,i) = vpa(int(E*A*dphi(i,1)*dphi(i,1),x,xh(i-1),xh(i)))...
    +vpa(int(E*A*dphi(i,2)*dphi(i,2),x,xh(i),xh(i+1)));

    S(i,i+1) = vpa(int(E*A*dphi(i,2)*dphi(i+1,1),x,xh(i),xh(i+1)));
end
S(n,n) = vpa(int(E*A*dphi(n,1)*dphi(n,1),x,xh(n-1),xh(n)));
S = S + S';
for i = 1:size(S)
    S(i,i) = S(i,i)/2;
end
Sv1 = S;

%mass matrix
M = zeros(n,n);
M(1,1) = vpa(int(phi(1,2)*phi(1,2),x,xh(1),xh(2)));
M(1,2) = vpa(int(phi(1,2)*phi(2,1),x,xh(1),xh(2)));
for i = 2:n-1 
    M(i,i) = vpa(int(phi(i,1)*phi(i,1),x,xh(i-1),xh(i)))...
    +vpa(int(phi(i,2)*phi(i,2),x,xh(i),xh(i+1)));

    M(i,i+1) = vpa(int(phi(i,2)*phi(i+1,1),x,xh(i),xh(i+1)));
end
M(n,n) = vpa(int(phi(n,1)*phi(n,1),x,xh(n-1),xh(n)));
M = M + M';
for i = 1:size(M)
    M(i,i) = M(i,i)/2;
end
M = mu*M;
Mv1 = M;

%inhomogeneity
syms f
f = 0*exp(x);
F = zeros(n,1);
F(1) = vpa(int(f*phi(1,2),x,xh(1),xh(2)));
for i = 2:n-1    
    F(i) = vpa(int(f*phi(i,1),x,xh(i-1),xh(i)))...
    +vpa(int(f*phi(i,2),x,xh(i),xh(i+1)));
end
F(n) = vpa(int(f*phi(n,1),x,xh(n-1),xh(n)));
%right hand side
fv1 = F;


%Assembling into a block matrix
S0 = blkdiag(Sv1,Sw1);
M0 = blkdiag(Mv1,Mw1);
S = S0;
M = M0;
F0 = [fw1;fv1];
F = F0;
for i = 1:N-1
    S = blkdiag(S,S0);
    M = blkdiag(M,M0);
    F = [F;F0]; 
end
%index of the beginning of the ith beam. 
%S(i)+1 is the first node of beam i 
ind = zeros(N+1,1);
for i = 1:N+1
    ind(i) = (i-1)*3*n;
end
%right hand side
F = [F;zeros(N_coupl*4+1*3,1)];
%coupling and boundary conditions
C = zeros(size(S,1),N_coupl*4+1*3);

%assembling of the coupling conditions
%junctions
for k = 1:N_coupl
    i = coupl(k,1);
    j = coupl(k,2);
    
    dk = (k-1)*4;
    
    C(ind(i+1),dk+1) = 1;
    C(ind(j)+n+2,dk+1) = -1;
    
    C(ind(j)+n,dk+2) = cos(alpha(j));
    C(ind(j+1)-1,dk+2) = -sin(alpha(j));
    
    C(ind(i)+n,dk+3) = cos(alpha(i));
    C(ind(i+1)-1,dk+3) = -sin(alpha(i));
    C(ind(j)+1,dk+3) = -cos(alpha(j));
    C(ind(j)+n+1,dk+3) = sin(alpha(j));
    
    C(ind(i)+n,dk+4) = sin(alpha(i));
    C(ind(i+1)-1,dk+4) = cos(alpha(i));
    C(ind(j)+1,dk+4) = -sin(alpha(j));
    C(ind(j)+n+1,dk+4) = -cos(alpha(j));
end   
%fixing beam 1
C(ind(1)+1,N_coupl*4+1) = 1;
C(ind(1)+n+1,N_coupl*4+2) = 1;
C(ind(1)+n+2,N_coupl*4+3) = 1;
%fixing beam 1
F(ind(N+1)+(N_coupl-1)*4+1) = x_pos(1,1);

%adding Neumann boundary conditions
%vertical
F(ind(bf+1)-2) = Qv;
F(ind(bf+1)-1) = Mv;
%horizontal
E = 1;
F(ind(bf)+n) = A*E*Qh;

%extended stiffness and mass matrix
S_e = [[S,C];
    [C',zeros(size(C,2),size(C,2))]];

M_e = blkdiag(M,zeros(size(C,2),size(C,2)));
%solution
z = S_e\F;
%solution without boundary conditions
beam_wb = z(1:end-size(C,2));

%plotting static case
figure(1)
for i = 1:N   
    beam = beam_wb(ind(i)+1:ind(i)+3*n);
    %horizontal solution
    v = beam(1:n);
    %vertical solution
    w = beam(n+1:end);
    w = w(1:2:end);
    %global coordinates
    xg = x_pos(1,i) + cos(alpha(i))*(xh+v)-sin(alpha(i))*w;
    yg = x_pos(2,i) + sin(alpha(i))*(xh+v)+cos(alpha(i))*w;    
    figure(1)
    plot(xg,yg,'b')
    hold on
end

%time evolution
u = z;
%inhomogeneity
p = zeros(size(S_e,1),1);

figure(2)
%Newmark method
beta = 0.25; %[0,0.5]
gamma = 0.5; %[0,1]
%first derivative
ud = zeros(size(S_e,1),1);
%second derivative
udd = zeros(size(S_e,1),1);
for j = 1:nt-1
    %j-dependent vectors
    v = u + ud*ht + (0.5-beta)*udd*ht^2;
    vd = ud + (1-gamma)*ht*udd;
    %left hand side of the system
    At = M_e + beta*S_e*ht^2;
    %right hand side of the system
    Bt = p - S_e*v;    
    udd = At\Bt;
    %update of u
    u = v + beta*udd*ht^2;
    %update of the derivative of u
    ud = vd + gamma*udd*ht;
    
    %plotting
    %solution without boundary conditions
    beam_wb = u(1:end-size(C,2));
    for i = 1:N
        beam = beam_wb(ind(i)+1:ind(i)+3*n);
        %horizontal solution
        v = beam(1:n);
        %vertical solution
        w = beam(n+1:end);
        w = w(1:2:end);
        %global coordinates
        xg = x_pos(1,i) + cos(alpha(i))*(xh+v)-sin(alpha(i))*w;
        yg = x_pos(2,i) + sin(alpha(i))*(xh+v)+cos(alpha(i))*w;
        figure(2)
        plot(xg,yg,'b')
        xlim([-2*L,3*L])
        ylim([-2*L,7*L])
        hold on        
    end
    hold off
    drawnow
     
end
    










    
    












