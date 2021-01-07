
function [Iter,Vx, Vy, Tole] = newton(data,step) 

maxIters = 12;                                                             % Max Iteration bound.................

%% Initial seed / germ for starting TENR Algorithm.........................

testcase = data;
mpc = loadcase(testcase);
mpc = ext2int(mpc);

%% Initializing ...........................................................

tnr = tnr_init(mpc);                                                       % Initial data format .......
[pv, pq, npv, npq]  = deal(tnr.pv, tnr.pq, tnr.npv, tnr.npq);              % PV, PQ buses and ..........

Vx = real(tnr.V0);                                                         % real component of complex voltage phasor..............
Vy = imag(tnr.V0);                                                         % imag component of complex voltage phasor..............
l = 1.0;                                                                   % Loadability parameter.................................

%% Newtpn Iterations ................................

Tole = 1;  
Iter = 1;
counter = 0; 

while (Tole > 1e-4)  

F = tnr.F( Vx,Vy,l,tnr );                                                  % Power mismatch vector ........

[J] = tnr.J(Vx,Vy,tnr);                                                    % Power flow Jacobian .......

dv = - J \ F;                                                               % Newton correction step ....

%% Updating state variables ..........

Vx([pv;pq])= Vx([pv;pq]) + step*dv(1:npv+npq);
Vy([pv;pq])= Vy([pv;pq]) + step*dv(npv+npq+1:end);

%% Checking Tolerance.......

Iter  = Iter + 1;
 Tole = max(abs(F));
 counter = counter + 1;
 if counter ==30
     break;
 end
end

if Iter > maxIters % check for non-convergence
    disp('Newton did not converge, change step size'); 
 else
     disp ('Newton converged within 12 iterations')
end
end
