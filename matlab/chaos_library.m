function [t,y]=chaos_library(S,t,IC,p)
% function chaos_library(S,t,IC,p)
% inputs  - S, string name of chaotic attractor
%         - t, time, either of the form [to,tf] or [to:f:tf]
%         - IC, initial conditions
%         - p, coefficients of the system of differential equations used.
% outputs - t, time
%         - y, column oriented time series
% Remarks
% - This MATLAB m file uses ODE solvers to calculate various chaotic 
%   attractors. The user is able to change all components of the chaotic 
%   attractor, including the time scale, the initial conditions and the 
%   internal parameters of the attractor equations. The m file will
%   also sketch out the plot bfore displaying a final version of it so the
%   user will be able to see the path a particle would take, and observe
%   patterns.
% - The first function handle contains the code that governs the input
%   recognition from the user. The rest of the handles contain more
%   information on the chaotic attractors, and includes the vector
%   components and typical parameters and initial conditions.
% - The discrete maps, Hennon and Logistic, take no step size. For those t
%   specifies the number of data points, either as a count or as a vector of
%   sample indices with one element per iterate. The first row of y is the
%   initial condition, so a request for n points returns n rows and n-1
%   iterations. The rest of the function inputs can be used normally.
% - y is column oriented for every system: n-by-3 for the flows, n-by-2 for
%   Hennon ([x y]) and n-by-1 for Logistic. t has one element per row of y.
% Future Work
% - More systems could be added.
% Jun 2016 - Created by Christopher Cunningham
% Jul 2016 - Modified by Ben Senderling, quarc@unomaha.edu
%          - Reformated comments section.
%          - Hennon, Logistic, Aizawa attractors added.
% Jul 2021 - Modified by Ben Senderling, quarc@unomaha.edu
%          - Removed plotting.
% Aug 2026 - Modified by Aaron Likens, alikens@unomaha.edu
%          - Fixed the two discrete maps. Logistic never assigned y and so
%            errored on every call; Hennon returned only its second
%            coordinate and rejected a vector t. Both now return a column
%            oriented y and a matching t. An unknown system name errors.
% Copyright (c) 2021-2026 Quantitative Analysis Research Core,
% Center for Human Movement Variability, University of Nebraska at Omaha.
% MIT licence. See LICENSE.txt.
%% Begin code

switch S
    case 'Rossler'
        [t,y] = ode45(@Rossler,t,IC,[],p);
    case 'Lorenz'
        [t,y] = ode45(@Lorenz,t,IC,[],p);
    case 'Hennon'
        % The map state is held in x and z, not in the output variable y.
        % Writing the second coordinate straight into y returned that
        % coordinate alone, as a row, with x discarded and the returned t
        % still the caller's iteration count.
        n = map_length(t);
        x = zeros(n,1);
        z = zeros(n,1);
        x(1) = IC(1);
        z(1) = IC(2);
        for c = 1:n-1
            x(c+1) = 1 - p(1)*x(c)^2 + z(c);
            z(c+1) = p(2)*x(c);
        end
        y = [x z];
        t = map_time(t,n);
    case 'Logistic'
        % The iterate was accumulated in x and never copied to y, so the
        % function threw MATLAB:unassignedOutputs on every call.
        n = map_length(t);
        x = zeros(n,1);
        x(1) = IC(1);
        for c = 1:n-1
            x(c+1) = p(1)*x(c)*(1-x(c));
        end
        y = x;
        t = map_time(t,n);
    case 'Aizawa'
        [t,y] = ode23(@Aizawa,t,IC,[],p);
    case 'DequanLi'
        
        [t,y] = ode23(@DequanLi,t,IC,[],p);
    case 'NoseHoover'
        [t,y] = ode45(@NoseHoover,t,IC,[],p);
    case 'QiChen'
        [t,y] = ode45(@QiChen,t,IC,[],p);
    case'YuWang'
        [t,y] = ode45(@YuWang,t,IC,[],p);
    case 'LuChen'
        [t,y] = ode45(@LuChen,t,IC,[],p);
    case 'Arneodo'
        [t,y] = ode45(@Arneodo,t,IC,[],p);
    case 'TSUCSI'
        [t,y] = ode45(@TSUCSI,t,IC,[],p);
    otherwise
        % Without this an unrecognised name fell out of the switch with both
        % outputs unassigned, which reads as a fault in the solver rather
        % than a misspelled system.
        error('chaos_library:unknownSystem', ...
            'Unknown system ''%s''.', S);
end
end

%% Support for the discrete maps

function n = map_length(t)
% Number of data points to generate for a discrete map. The maps take no
% step size, so t is either a count or a vector of sample indices; each
% element of a vector is one iterate. A vector reached the colon expression
% directly before this and raised MATLAB:colon:operandsNotRealScalar, even
% though the header documents t as being allowed in vector form.
if isscalar(t)
    n = round(t);
else
    n = numel(t);
end
if n < 1
    error('chaos_library:badLength', ...
        't must ask for at least one data point.');
end
end

function t = map_time(t,n)
% Sample indices for a discrete map, one per row of y, so the two outputs
% can be used together. A vector t is returned as given, in column form.
if isscalar(t)
    t = (0:n-1)';
else
    t = t(:);
end
end

function [yp] = Rossler(t,y,p)
% function [yp] = Rossler(t,y,p)
% inputs  - 
% outputs - 
% Remarks
% - Typical IC's are: [-9 0 0]
% - Typical Parameters are: [0.2 0.2 5.7]
% - Source: https://en.wikipedia.org/wiki/R%C3%B6ssler_attractor
% - Equations
%       x'= -y -z
%       y'= x + ay
%       z'= b+ z(x-c)
%% Begin code

yp = zeros(3,1);
yp(1) = -y(2)-y(3);
yp(2) = y(1) + p(1)*y(2);
yp(3) = p(2) + y(3)*(y(1)-p(3));

end

%LORENZ================================================
function [ yp ] = Lorenz(t,y,p)
%   Typical IC's are: 0 -0.01 9
%   Typical Params are: 10 28 8/3
%   Source: https://en.wikipedia.org/wiki/Lorenz_system
%   EQUATION
%       x' = o(y-x)
%       y' = x(p-z)-y
%       z' = xy - Bz
%% Begin code

yp = zeros(3,1);
yp(1) = p(1)*( y(2) - y(1) );
yp(2) = -y(1).*y(3) + p(2).*y(1) - y(2);
yp(3) = (y(1).*y(2)) - p(3)*y(3);

end

%DEQUAN-LI=============================================
function [ yp ] = DequanLi(t,y,p)
%   Typical IC's are: Unknown, start with a point other than [0 0 0]
%   Typical Params are: 40 0.16 55 20 1.833 0.65
%   Source: https://www.researchgate.net/publication/223710675_A_three-scroll_chaotic_attractor
%   EQUATION
%       x' =  a(y-x)+dxz
%       y' = px + Ly - xz
%       z = Bz + xy = ex^3
%% Begin code

yp = zeros(3,1);
yp(1) = p(1)*(y(2) - y(1)) + p(2)*(y(1)*y(3));
yp(2) = p(3)*y(1) + p(4)*y(2) - y(1)*y(3);
yp(3) = p(5)*y(3) + y(1)*y(2) - p(6)*(y(1)*y(1));
end

%NOSE-HOOVER===========================================
function [ yp ] = NoseHoover(t,y,p)
%   Typical IC's are: Unknown, start with a point other than [0 0 0]
%   Typical Params are: 1.5
%   Source: http://williamhoover.info/Scans1980s/1986-4.pdf
%   EQUATION:
%       x' = y
%       y' = -x +yz
%       z' = a - y^2
%% Begin code

yp = zeros(3,1);
yp(1) = y(2);
yp(2) = -y(1) + y(2)*y(3);
yp(3) = p(1) - (y(2)*y(2));
end

%QI-CHEN===============================================
function [ yp ] = QiChen(t,y,p)
%   Typical IC's are: Unknown, start with a point other than [0 0 0]
%   Typical Params are: 38 8/3 80
%   Source: https://www.emis.de/journals/HOA/MPE/Volume2012/438328.pdf
%   EQUATION
%       x' = a(y-x) + yz
%       y' = cx + y -xz
%       z' = xy - Bz
%% Begin code

yp = zeros(3,1);
yp(1) = p(1)*(y(2)-y(1)) + y(2)*y(3);
yp(2) = p(3)*y(1) + y(2) - y(1)*y(3);
yp(3) = y(1)*y(2) - p(2)*y(3);
end

%YU-WANG===============================================
function [ yp ] = YuWang(t,y,p)
%   Typical IC's are: Unknown, start with a point other than [0 0 0]
%   Typical Params are: 10 40 2 2.5
%   EQUATION:
%       x'= a(y-x)
%       y' = Bx - cxz
%       z' = e^(xy) - dz
%% Begin code

yp = zeros(3,1);
yp(1) = p(1)*(y(2)-y(1));
yp(2) = p(2)*y(1) - p(3)*y(1)*y(3);
yp(3) = exp(y(1)*y(2)) - p(4)*y(3);
end

%LU-CHEN===============================================
function [ yp ] = LuChen(t,y,p)
%   Typical IC's are: Unknown, start with a point other than [0 0 0]
%   Typical Params are: -10 -4 18.1
%   Source: https://en.wikipedia.org/wiki/Multiscroll_attractor
%   EQUATION
%       x' = -((aBx)/(a+B))-yz + c
%       y' = ay + xz
%       z' = Bz + xy
%% Begin code

yp = zeros(3,1);
yp(1) = -1*(p(1)*p(2)*y(1)/(p(1)+p(2))) - y(2)*y(3) + p(3);
yp(2) = p(1)*y(2) + y(1)*y(3);
yp(3) = p(2)*y(3) + y(1)*y(2);
end

%ARNEODO===============================================
function [ yp ] = Arneodo(t,y,p)
%   Typical IC's are: Unknown, start with a point other than [0 0 0]
%   Typical Params are: -5.5 3.5 -1
%   EQUATION
%       x' = y
%       y' = z
%       z' = =ax -By - z + dx^3
%% Begin code

yp = zeros(3,1);
yp(1)= y(2);
yp(2) = y(3);
yp(3) = -p(1)*y(1) - p(2)*y(2)-y(3) + p(3)*y(1)*y(1)*y(1);
end

%TSUCSI================================================
function [yp] = TSUCSI(t,y,p)
%   Typical IC's are: Unknown, start with a point other than [0 0 0]
%   Typical Params are: 40 0.833 0.5 0.65 20
%   Source: https://www.researchgate.net/publication/265811017_A_New_Three-Scroll_Unified_Chaotic_System_Coined
%   EQUATION
%       x' = a(y-x) +dxz
%       y' = Ly - xz
%       z' = Bz + xy - ex^2
%% Begin code

yp = zeros(3,1);
yp(1) = p(1)*(y(2)-y(1)) + p(3)*y(1)*y(3);
yp(2) = p(5)*y(2) - y(1)*y(3);
yp(3) = p(2)*y(3) + y(1)*y(2) - p(4)*y(1)*y(1);
end

%AIZAWA================================================
function [yp] =Aizawa(t,y,p)
%   Typical IC's are: Unknown, start with a point other than [0 0 0]
%   Typical Params are: 0.95 0.9 0.6 3.5 0.25 0.1 
%   EQUATION
%       x' = (z-b)*x - d*y
%       y' = d*y + (z-b)*y 
%       z' = c + a*z - (z^3/3)+f*z*x^3
%% Begin code

yp=zeros(3,1);
yp(1)=(y(3)-p(2))*y(1)-p(4)*y(2);
yp(2)=p(4)*y(1)+(y(3)-p(2))*y(2);
yp(3)=p(3)+p(1)*y(3)-(y(3))^3/3-(y(1))^2+p(6)*y(3)*(y(1))^3;

end

%HENNON================================================
%   Typical IC's are: Unknown, start with a point other than [0 0]
%   Typical Params are: 1.4 0.3
%   Source: https://en.wikipedia.org/wiki/H%C3%A9non_map
%   EQUATION
%       x(n+1) = 1 - ax(n)^2 + y(n)
%       y(n+1) = bx(n)

%LOGISTIC==============================================
%   Typical IC is: Use a value between 0 and 1
%   Typical Params are: Any positive value, usually between 0 and 4
%   Source: https://en.wikipedia.org/wiki/Logistic_map
%   EQUATION
%       x(n+1) = r*x(n)*(1-x(n))




