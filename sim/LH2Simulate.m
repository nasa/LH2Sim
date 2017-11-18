function data = LH2Simulate(tspan,name)
% data = LH2Simulate(tspan,name)
%	Simulates LH2 model
%	tspan is optional argument specifying time array. If not given,
%		or given as [], then defaults to [0:1:LH2Model.tFinal].
%	name is an optional string argument specifying scenario name.

% obtain model parameters structure
try
	P = evalin('base','LH2Model');
catch ME
	if strcmp(ME.identifier,'MATLAB:UndefinedFunction')
		evalin('base','LH2ModelParams');
		P = evalin('base','LH2Model');
	else
		error(ME.message);
	end
end

% set default tspan
if nargin<1 || isempty(tspan)
	tspan = 0:1:P.tFinal;
end

% set default name
if nargin<2
	name = 'LH2';
end

P.tEnd = tspan(end);
P.waitbar = waitbar(0,['Simulating ' name '...']);

% set up initial state
mL0 = P.mL0;
TLB0 = P.TL0*ones(P.k,1);
%TLL0 = mean([P.TL0 P.Tw0])*ones(P.k-2,1);
TLL0 = P.TL0*ones(P.k-2,1);
TLI0 = P.TL0*ones(P.nL,1);
mv0 = P.mv0;
mg0 = P.mg0;
TvB0 = P.Tv0*ones(P.n,1);
%TvL0 = mean([P.Tv0 P.Tw0])*ones(P.n-2,1);
TvL0 = P.Tv0*ones(P.n-2,1);
TvI0 = P.Tv0*ones(P.nV,1);
Ts0 = 20; %P.T_c*(P.p0/P.p_c)^(1/P.lambda);
x0 = [	mL0;
		TLB0;
		TLL0;
		TLI0;
		mv0;
		mg0;
		TvB0;
		TvL0;
		TvI0;
		Ts0;
		P.TwL0;
		P.Twv0;
		P.Tv0;
	];

% simulate equations
P.relTol = 1e-4; %doesnt work well with 1e-3
options = odeset('MaxStep',1,'RelTol',P.relTol);
rhs = @(t,x) LH2dxdt(P,t,x);
[t,x] = ode15s(rhs,tspan,x0,options);

% close waitbar
close(P.waitbar);

% configure data struct
data.name = name;
data.t = t;

% extract state variables
data.mL = x(:,1);
data.TLB = x(:,2:1+P.k);
data.TLL = x(:,2+P.k:1+2*P.k-2);
data.TLI = x(:,2+2*P.k-2:1+2*P.k-2+P.nL);
data.mv = x(:,2+2*P.k-2+P.nL);
data.mg = x(:,3+2*P.k-2+P.nL);
data.TvB = x(:,4+2*P.k-2+P.nL:3+2*P.k-2+P.nL+P.n);
data.TvL = x(:,4+2*P.k-2+P.nL+P.n:3+2*P.k-2+P.nL+2*P.n-2);
data.TvI = x(:,4+2*P.k-2+P.nL+2*P.n-2:3+2*P.k-2+P.nL+2*P.n-2+P.nV);
data.Ts = x(:,4+2*P.k-2+P.nL+2*P.n-2+P.nV);
data.TwL = x(:,5+2*P.k-2+P.nL+2*P.n-2+P.nV);
data.Twv = x(:,6+2*P.k-2+P.nL+2*P.n-2+P.nV);
data.TvG = x(:,7+2*P.k-2+P.nL+2*P.n-2+P.nV);

TLL = [data.TLB(:,1) data.TLL data.TLB(:,P.k)];
TvL = [data.TvB(:,1) data.TvL data.TvB(:,P.n)];

% tank calculations for liquid
VL = data.mL/P.rhoL;
hL = VL/P.A;
hLi = hL/P.k;
pL = P.rhoL*P.g*hL;
nuL = P.muL/P.rhoL;			% note, this is not Nusselt number (Nu)
PrL = P.muL*P.c_L/P.kappaL;
PsiL = (1+(0.492/PrL)^(9/16))^(-16/9);

% tank calculations for vapor
VUllage = P.VTotal-VL;
hvi = (P.H-hL)/P.n;
Prv = P.mu*P.c_p/P.kappa_v;
rhov = data.mv./VUllage;
rhog = data.mg./VUllage;
rho = rhov+rhog;
nuv = P.mu./rho;				% note, this is not Nusselt number (Nu)
Psiv = (1+(0.492/Prv)^(9/16))^(-16/9);

% compute the deltaL terms
for i=1:P.k
	xL(:,i) = (data.TwL>=data.TLB(:,i)).*hLi*i + (data.TwL<data.TLB(:,i)).*(hL-hLi*i);
	GrL(:,i) = P.g*P.betaL*(data.TwL-data.TLB(:,i)).*xL(:,i).^3./nuL^2;
	RaL(:,i) = abs(GrL(:,i))*PrL;
	deltaL(:,i) = (RaL(:,i)<10^9).*xL(:,i).*3.93.*((0.952+PrL)./abs(GrL(:,i))./PrL.^2).^(1/4) + ...
		(RaL(:,i)>=10^9).*xL(:,i).*0.565.*((1+0.494*PrL.^(2/3))./abs(GrL(:,i))).^(1/10)./PrL.^(8/15);
end

% compute the deltav terms
for i=1:P.n
	xv(:,i) = (data.Twv>=data.TvB(:,i)).*hvi*i + (data.Twv<data.TvB(:,i)).*((P.H-hL)-hvi*i);
	Grv(:,i) = P.g*1./TvL(:,i).*(data.Twv-data.TvB(:,i)).*xv(:,i).^3./nuv.^2;
	Rav(:,i) = abs(Grv(:,i))*Prv;
	deltav(:,i) = (Rav(:,i)<10^9).*xv(:,i).*3.93.*((0.952+Prv)./abs(Grv(:,i))./Prv.^2).^(1/4) + ...
		(Rav(:,i)>=10^9).*xv(:,i).*0.565.*((1+0.494*Prv.^(2/3))./abs(Grv(:,i))).^(1/10)./Prv.^(8/15);
end

% compute liquid volumes for bulk and boundary layers (first is 0 for
% convenience, there is no boundary layer 1)
VLB(:,1) = VL./P.k;
for i=2:P.k-1
	VLB(:,i) = pi*(P.R-deltaL(:,i)).^2.*hLi;
	VLL(:,i) = (pi*P.R^2-pi*(P.R-deltaL(:,i)).^2).*hLi;
end
VLB(:,P.k) = VL./P.k;

% compute vapor/gas volumes for bulk and boundary layers (first is 0 for
% convenience, there is no boundary layer 1)
VvB(:,1) = VUllage./P.n;
for i=2:P.n-1
	VvB(:,i) = pi*(P.R-deltav(:,i)).^2.*hvi;
	VvL(:,i) = (pi*P.R^2-pi*(P.R-deltav(:,i)).^2).*hvi;
end
VvB(:,P.n) = VUllage./P.n;
VvL(:,P.n) = 0; % doing this only to make vectors same length for vector operations below, will ignore this value)

sumv = 0;
for i=1:size(VvB,2)
	sumv = sumv + VvB(:,i)./data.TvB(:,i) + VvL(:,i)./TvL(:,i);
end
data.pv = data.mv.*P.R_v./sumv;
data.pg = data.mg.*P.R_g./sumv;
data.p = data.pv + data.pg;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% differential equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dxdt = LH2dxdt(P,t,x)

%fprintf('Time = %g\n',t);
waitbar(t/P.tEnd,P.waitbar);

% obtain state variables
mL = x(1);
TLB = x(2:1+P.k);
TLL = x(2+P.k:1+2*P.k-2);
TLI = x(2+2*P.k-2:1+2*P.k-2+P.nL);
mv = x(2+2*P.k-2+P.nL);
mg = x(3+2*P.k-2+P.nL);
TvB = x(4+2*P.k-2+P.nL:3+2*P.k-2+P.nL+P.n);
TvL = x(4+2*P.k-2+P.nL+P.n:3+2*P.k-2+P.nL+2*P.n-2);
TvI = x(4+2*P.k-2+P.nL+2*P.n-2:3+2*P.k-2+P.nL+2*P.n-2+P.nV);
Ts = x(4+2*P.k-2+P.nL+2*P.n-2+P.nV);
TwL = x(5+2*P.k-2+P.nL+2*P.n-2+P.nV);
Twv = x(6+2*P.k-2+P.nL+2*P.n-2+P.nV);
TvG = x(7+2*P.k-2+P.nL+2*P.n-2+P.nV);

if mL<=1e-5
	error('Out of LH2 at t=%g!',t);
end

% setup temperatures of boundary layers (setting first and last to bulks
% for convenience)
TLL = [0; TLL; 0];
TvL = [0; TvL; 0];
% TLL = (Tw+TLB)/2;
% TvL = (Tw+TvB)/2;
TLL(1) = TLB(1);
TLL(P.k) = TLB(P.k);
TvL(1) = TvB(1);
TvL(P.n) = TvB(P.n);

% tank calculations for liquid
VL = mL/P.rhoL;
hL = VL/P.A;
hLi = hL/P.k;
pL = P.rhoL*P.g*hL;
nuL = P.muL/P.rhoL;			% note, this is not Nusselt number (Nu)
PrL = P.muL*P.c_L/P.kappaL;
PsiL = (1+(0.492/PrL)^(9/16))^(-16/9);

% tank calculations for vapor
VUllage = P.VTotal-VL;
hvi = (P.H-hL)/P.n;
Prv = P.mu*P.c_p/P.kappa_v;
rhov = mv/VUllage;
rhog = mg/VUllage;
rho = rhov+rhog;
nuv = P.mu/rho;				% note, this is not Nusselt number (Nu)
Psiv = (1+(0.492/Prv)^(9/16))^(-16/9);

% compute the deltaL terms
for i=1:P.k
	if TwL>=TLB(i)
		xL(i) = hLi*i;
	else
		xL(i) = hL-hLi*i;
	end
	GrL(i) = P.g*P.betaL*(TwL-TLB(i))*xL(i)^3/nuL^2;
	RaL(i) = abs(GrL(i))*PrL;
	if RaL(i)<10^9
		% laminar
		deltaL(i) = xL(i)*3.93*((0.952+PrL)/abs(GrL(i))/PrL^2)^(1/4);
	else
		% turbulent
		deltaL(i) = xL(i)*0.565*((1+0.494*PrL^(2/3))/abs(GrL(i)))^(1/10)/PrL^(8/15);
	end
end

% compute the deltav terms
for i=1:P.n
	if Twv>=TvB(i)
		xv(i) = hvi*i;
	else
		xv(i) = (P.H-hL)-hvi*i;
	end
	Grv(i) = P.g*1/TvL(i)*(Twv-TvB(i))*xv(i)^3/nuv^2;
	Rav(i) = abs(Grv(i))*Prv;
	if Rav(i)<10^9
		% laminar
		deltav(i) = xv(i)*3.93*((0.952+Prv)/abs(Grv(i))/Prv^2)^(1/4);
	else
		% turbulent
		deltav(i) = xv(i)*0.565*((1+0.494*Prv^(2/3))/abs(Grv(i)))^(1/10)/Prv^(8/15);
	end
end

% compute liquid volumes for bulk and boundary layers (first is 0 for
% convenience, there is no boundary layer 1)
VLB(1) = VL/P.k;
for i=2:P.k-1
	% bulk layer
	VLB(i) = pi*(P.R-deltaL(i))^2*hLi;
	% boundary layer
	VLL(i) = (pi*P.R^2-pi*(P.R-deltaL(i))^2)*hLi;
end
VLB(P.k) = VL/P.k;

% compute vapor/gas volumes for bulk and boundary layers (first is 0 for
% convenience, there is no boundary layer 1)
VvB(1) = VUllage/P.n;
for i=2:P.n-1
	% bulk layer
	VvB(i) = pi*(P.R-deltav(i))^2*hvi;
	% boundary layer
	VvL(i) = (pi*P.R^2-pi*(P.R-deltav(i))^2)*hvi;
end
VvB(P.n) = VUllage/P.n;
VvL(P.n) = 0; % doing this only to make vectors same length for vector operations below, will ignore this value)

% compute liquid masses
mLB = VLB*P.rhoL;
mLL = VLL*P.rhoL;

% compute common pressures (p = m/(V1/(RT1)+V2/(RT2)+...))
pv = mv*P.R_v/(sum(VvB./TvB')+sum(VvL./TvL'));
pg = mg*P.R_g/(sum(VvB./TvB')+sum(VvL./TvL'));
p = pv+pg;

% compute local masses and densities
mvB = pv*VvB/P.R_v./TvB';
mgB = pg*VvB/P.R_g./TvB';
mvL = pv*VvL/P.R_v./TvL';
mgL = pg*VvL/P.R_g./TvL';
rhovB = mvB./VvB;
rhogB = mgB./VvB;
rhovL = mvL./VvL;
rhogL = mgL./VvL;
rhovL(1) = rhovB(1);
rhogL(1) = rhogB(1);

% flow out tank bottom
JLOut = P.JLOut; 

% determine new vent valve state
P.VentState = evalin('base','LH2Model.VentState;');
P.VentState = getVentState(P,p);

% compute vapor/gas flows for ET vent
%Jve = rhov*P.SVent*P.VentState*sqrt(P.gamma*(p-P.pAtm))/P.Gamma/sqrt(P.KVent*p);
%Jge = rhog*P.SVent*P.VentState*sqrt(P.gamma*(p-P.pAtm))/P.Gamma/sqrt(P.KVent*p);
gamma = (mv*P.c_p+mg*P.c_gp)/(mv*P.c_v+mg*P.c_gv);
Jve = P.VentState*gasFlow(P.SVent,gamma,rhov,p,P.pAtm);
Jge = P.VentState*gasFlow(P.SVent,gamma,rhog,p,P.pAtm) - P.Pressurization*gasFlow(P.SVent,gamma,rhog,3*P.pAtm,p);

% surface temperature
Ts0 = P.T_c*(pv/P.p_c)^(1/P.lambda);	% note this is function of pv, not p
dTsdt = (Ts0-Ts)/P.tmin;

% optimal grid for liquid (count from interface to bottom)
lmin = sqrt(P.kappaL*P.tmin/P.c_L/P.rhoL);
l_L(1) = lmin/(1+exp(pi/2/sqrt(P.nL+1)));			% l_0
l12_L(1) = lmin;									% l_1/2
for i=2:P.nL+1
	l12_L(i) = l12_L(i-1)*exp(pi/sqrt(P.nL+1));		% l_i+1/2
	l_L(i) = sqrt(l12_L(i-1)*l12_L(i));				% l_i
end

% optimal grid for vapor/gas (count from interface to top)
lmin = sqrt(P.kappa_v*P.tmin/P.c_v/rho);
l_V(1) = lmin/(1+exp(pi/2/sqrt(P.nV+1)));			% l_0
l12_V(1) = lmin;									% l_1/2
for i=2:P.nV+1
	l12_V(i) = l12_V(i-1)*exp(pi/sqrt(P.nV+1));		% l_i+1/2
	l_V(i) = sqrt(l12_V(i-1)*l12_V(i));				% l_i
end

% heat transfer liquid->surface
alphaL_cond = P.kappaL/l12_L(1);
QdotLS_cond = alphaL_cond*P.A*(TLI(1)-Ts) - l_L(1)*P.c_L*P.rhoL*dTsdt;
%GrLs = P.g*P.betaL*(TLB(P.k)-Ts)*P.R^3/nuL^2;
GrLs = P.g*P.betaL*(TLI(1)-Ts)*P.R^3/nuL^2;
RaLs = abs(GrLs)*PrL;
if RaLs>1e5 && RaLs<2e9
	alphaL_conv = P.kappaL/P.R*(0.54*(RaLs)^(1/4));
else %if RaLs>=1e9 && RaLs<3e10
	alphaL_conv = P.kappaL/P.R*(0.15*(RaLs)^(1/3));
%else
%	error('Bad RaLs at time %g',t);
end
QdotLS_conv = pi*P.R^2*alphaL_conv*(TLB(P.k)-Ts);
if TLI(1)>Ts %TLB(P.k)>Ts
	QdotLS = max(QdotLS_conv,QdotLS_cond);
else
	QdotLS = QdotLS_cond; % Q_dotLS2_conv is 0 here
end

% heat transfer vapor->surface
alphav_cond = P.kappa_v/l12_V(1);
QdotVS_cond = alphav_cond*P.A*(TvI(1)-Ts) - l_V(1)*P.c_v*rho*dTsdt;
%Grvs = P.g*1/TvB(1)*(Ts-TvB(1))*P.R^3/nuv^2;
Grvs = P.g*1/TvB(1)*(Ts-TvI(1))*P.R^3/nuv^2;
Ravs = abs(Grvs)*Prv;
if Ravs>1e5 && Ravs<2e9
	alphav_conv = P.kappa_v/P.R*(0.54*(Ravs)^(1/4));
else %if RaLs>=1e9 && RaLs<3e10
	alphav_conv = P.kappa_v/P.R*(0.15*(Ravs)^(1/3));
%else
%	error('Bad RaLs at time %g',t);
end
QdotVS_conv = pi*P.R^2*alphav_conv*(TvB(1)-Ts);
if Ts>TvI(1) %Ts>TvB(1)
	QdotVS = min(QdotVS_conv,QdotVS_cond);
else
	QdotVS = QdotVS_cond; % Q_dotVS2_conv is 0 here
end

% condensation flow
qh = P.qh0*sqrt((P.T_c-Ts)/(P.T_c-P.TL0));
Jcd = -(QdotLS+QdotVS)/qh;

% heat transfer wall->liquid
if RaL(1)<2*10^9 % && RaL(1)>10^5 % seems I can get below the lower limit
	alphaL = P.kappaL/P.R*0.59*RaL(1)^(0.24);
elseif RaL(1)>=2*10^9 %&& RaL(1)<2*10^14
	alphaL = P.kappaL/P.R*0.14*RaL(1)^(0.33);
else
	error('Bad RaL=%g at time %g',RaL(1),t);
end
QdotWL1 = pi*P.R^2*alphaL*(TwL-TLB(1));
for i=2:P.k-1
	Si = 2*pi*P.R*hLi;
	if RaL(i)>10^4 && RaL(i)<10^9
		Nui = 0.68+0.503*(RaL(i)*PsiL)^(1/4);
	else %if RaL(i)>=10^9 && RaL(i)<10^14
		Nui = 0.15*RaL(i)^(1/3);
	%else
	%	error('Bad Rav at %g',t);
	end
	alphai = P.kappaL/xL(i)*Nui;
	QdotWL(i) = alphai*Si*(TwL-TLB(i));
end
QdotWLTotal = sum(QdotWL)+QdotWL1;

% heat transfer wall->vapor/gas
if Rav(P.n)>1e5 && Rav(P.n)<2e9
	alphav = P.kappa_v/P.R*(0.54*(Rav(P.n))^(1/4));
else
	alphav = P.kappa_v/P.R*(0.15*(Rav(P.n))^(1/3));
end
QdotWVn = pi*P.R^2*alphav*(Twv-TvB(P.n));
for i=2:P.n-1
	Si = 2*pi*P.R*hvi;
	if Rav(i)>10^4 && Rav(i)<10^9
		Nui = 0.68+0.503*(Rav(i)*Psiv)^(1/4);
	else %if Rav(i)>=10^9 && Rav(i)<10^14
		Nui = 0.15*Rav(i)^(1/3);
	%else
	%	error('Bad Rav at %g',t);
	end
	alphai = P.kappa_v/xv(i)*Nui;
	QdotWV(i) = alphai*Si*(Twv-TvB(i));
end
QdotWVTotal = sum(QdotWV)+QdotWVn;

% heat transfer environment->wall
QdotEW = P.QdotEW;

% total mass flows
JLTotal = Jcd - JLOut;
dVLdt = JLTotal/P.rhoL;
dVUllagedt = -dVLdt;
mdotv = -Jcd - Jve;
mdotg = -Jge;

% compute the JLi terms
for i=1:P.k-1
	Vi = sign(GrL(i))*1.185*nuL/xL(i)*sqrt(abs(GrL(i))/(1+.494*PrL^(2/3)));
	if RaL(i)<10^9
		% laminar
		JLL(i) = 2*pi*P.R*P.rhoL*Vi*deltaL(i)*0.0833;
	else
		% turbulent
		JLL(i) = 2*pi*P.R*P.rhoL*Vi*deltaL(i)*0.1436;
	end
end

% compute JLBL mass flows (first 0 for convenience)
for i=2:P.k-1
	JLBL(i) = VLL(i)/VL*JLTotal - JLL(i-1) + JLL(i);
end

% compute JLB mass flows (first 0)
JLB(2) = VLB(1)/VL*JLTotal + JLL(1) + JLOut;
for i=3:P.k
	JLB(i) = VLB(i-1)/VL*JLTotal + JLBL(i-1) + JLB(i-1);
end

% compute mdot for bottom liquid layer
mdotLB(1) = JLB(2) - JLL(1) - JLOut;
% compute mdot for middle bulk and boundary/lateral liquid layers
for i=2:P.k-1
	% bulk
	mdotLB(i) = JLB(i+1) - JLB(i) - JLBL(i);
	% boundary
	mdotLL(i) = JLL(i-1) + JLBL(i) - JLL(i);
end
% compute mdot for top liquid layer
mdotLB(P.k) = JLL(P.k-1) - JLB(P.k) + Jcd;

% compute the Jvi,Jgi terms
for i=1:P.n-1
	Vi = sign(Grv(i))*1.185*nuv/xv(i)*sqrt(abs(Grv(i))/(1+.494*Prv^(2/3)));
	if Rav(i)<10^9
		% laminar
		JvL(i) = 2*pi*P.R*(rhovL(i))*Vi*deltav(i)*0.0833;
		JgL(i) = 2*pi*P.R*(rhogL(i))*Vi*deltav(i)*0.0833;
	else
		% turbulent
		JvL(i) = 2*pi*P.R*(rhovL(i))*Vi*deltav(i)*0.1436;
		JgL(i) = 2*pi*P.R*(rhogL(i))*Vi*deltav(i)*0.1436;
	end
end

% compute JvBL,JgBL mass flows (first 0 for convenience)
for i=2:P.n-1
	JvBL(i) = VvL(i)/VUllage*mdotv - JvL(i-1) + JvL(i);
	JgBL(i) = VvL(i)/VUllage*mdotg - JgL(i-1) + JgL(i);
end

% compute JvB,JgB mass flows (first 0)
JvB(2) = VvB(1)/VUllage*mdotv + JvL(1) + Jcd;
JgB(2) = VvB(1)/VUllage*mdotg + JgL(1);
for i=3:P.n
	JvB(i) = VvB(i-1)/VUllage*mdotv + JvBL(i-1) + JvB(i-1);
	JgB(i) = VvB(i-1)/VUllage*mdotg + JgBL(i-1) + JgB(i-1);
end

% compute mdot for bottom vapor/gas layer
mdotvB(1) = JvB(2) - JvL(1) - Jcd;
% compute mdot for middle bulk and boundary/lateral vapor/gas layers
for i=2:P.n-1
	% bulk
	mdotvB(i) = JvB(i+1) - JvB(i) - JvBL(i);
	mdotgB(i) = JgB(i+1) - JgB(i) - JgBL(i);
	% boundary
	mdotvL(i) = JvL(i-1) + JvBL(i) - JvL(i);
	mdotgL(i) = JgL(i-1) + JgBL(i) - JgL(i);
end
% compute mdot for top vapor/gas layer
mdotvB(P.n) = JvL(P.n-1) - JvB(P.n) - Jve;
mdotgB(P.n) = JgL(P.n-1) - JgB(P.n) - Jge;

% compute enthalpies for mass flows (liquid)
for i=2:P.k
	hLB(i) = P.c_L*( (JLB(i)>0)*TLB(i) + (JLB(i)<0)*TLB(i-1) );
end
for i=2:P.k-1
	hLBL(i) = P.c_L*( (JLBL(i)>0)*TLB(i) + (JLBL(i)<0)*TLL(i) );
end
hLL(1) = P.c_L*( (JLL(1)>0)*TLB(1) + (JLL(1)<0)*TLL(2) );
for i=2:P.k-1
	hLL(i) = P.c_L*( (JLL(i)>0)*TLL(i) + (JLL(i)<0)*TLL(i+1) );
end

% compute enthalpies for mass flows (vapor/gas)
for i=2:P.n
	hvB(i) = P.c_p*( (JvB(i)>0)*TvB(i) + (JvB(i)<0)*TvB(i-1) );
	hgB(i) = P.c_gp*( (JgB(i)>0)*TvB(i) + (JgB(i)<0)*TvB(i-1) );
end
for i=2:P.n-1
	hvBL(i) = P.c_p*( (JvBL(i)>0)*TvB(i) + (JvBL(i)<0)*TvL(i) );
	hgBL(i) = P.c_gp*( (JgBL(i)>0)*TvB(i) + (JgBL(i)<0)*TvL(i) );
end
hvL(1) = P.c_p*( (JvL(1)>0)*TvB(1) + (JvL(1)<0)*TvL(2) );
hgL(1) = P.c_gp*( (JgL(1)>0)*TvB(1) + (JgL(1)<0)*TvL(2) );
for i=2:P.n-1
	hvL(i) = P.c_p*( (JvL(i)>0)*TvL(i) + (JvL(i)<0)*TvL(i+1) );
	hgL(i) = P.c_gp*( (JgL(i)>0)*TvL(i) + (JgL(i)<0)*TvL(i+1) );
end

% compute Tdot for bottom liquid layer
dTLBdt(1) = 1/mLB(1)/P.c_L*( ...
	QdotWL1 ...
	+ JLB(2)*hLB(2) ...
	- JLL(1)*hLL(1) ...
	- JLOut*P.c_L*TLB(1) ...
	- p*dVLdt*VLB(1)/VL ...
	- mdotLB(1)*P.c_L*TLB(1) );
% compute Tdot for middle bulk liquid layers
for i=2:P.k-1
	dTLBdt(i) = 1/mLB(i)/P.c_L*( ...
		JLB(i+1)*hLB(i+1) ...
		- JLB(i)*hLB(i) ...
		- JLBL(i)*hLBL(i) ...
		- p*dVLdt*VLB(i)/VL ...
		- mdotLB(i)*P.c_L*TLB(i));
end
% compute Tdot for middle boundary liquid layers
for i=2:P.k-1
	dTLLdt(i-1) = 1/mLL(i)/P.c_L*( ...
		QdotWL(i) ...
		+ JLL(i-1)*hLL(i-1) ...
		- JLL(i)*hLL(i) ...
		+ JLBL(i)*hLBL(i) ...
		- p*dVLdt*VLL(i)/VL ...
		- mdotLL(i)*P.c_L*TLL(i));
end
% compute Tdot for top liquid layer
dTLBdt(P.k) = 1/mLB(P.k)/P.c_L*( ...
	- QdotLS ...
	+ JLL(P.k-1)*hLL(P.k-1) ...
	- JLB(P.k)*hLB(P.k) ...
	+ Jcd*P.c_p*Ts ...
	- p*dVLdt*VLB(P.k)/VL ...
	- mdotLB(P.k)*P.c_L*TLB(P.k));

% compute Tdot for bottom vapor/gas layer
dTvBdt(1) = 1/(mvB(1)*P.c_v+mgB(1)*P.c_gv)*( ...
	- QdotVS ...
	+ (JvB(2)*hvB(2)+JgB(2)*hgB(2)) ...
	- (JvL(1)*hvL(1)+JgL(1)*hgL(1)) ...
	- Jcd*P.c_p*Ts ...
	- p*dVUllagedt*VvB(1)/VUllage ...
	- (mdotvB(1)*P.c_v*TvB(1)+mdotgB(1)*P.c_gv*TvB(1)) );
% compute Tdot for middle bulk vapor/gas layers
for i=2:P.n-1
	dTvBdt(i) = 1/(mvB(i)*P.c_v+mgB(i)*P.c_gv)*( ...
		(JvB(i+1)*hvB(i+1)+JgB(i+1)*hgB(i+1)) ...
		- (JvB(i)*hvB(i)+JgB(i)*hgB(i)) ...
		- (JvBL(i)*hvBL(i)+JgBL(i)*hgBL(i)) ...
		- p*dVUllagedt*VvB(i)/VUllage ...
		- (mdotvB(i)*P.c_v*TvB(i) + mdotgB(i)*P.c_gv*TvB(i)) );
end
% compute Tdot for middle boundary vapor/gas layers
for i=2:P.n-1
	dTvLdt(i-1) = 1/(mvL(i)*P.c_v+mgL(i)*P.c_gv)*( ...
		QdotWV(i) ...
		+ (JvL(i-1)*hvL(i-1)+JgL(i-1)*hgL(i-1)) ...
		- (JvL(i)*hvL(i)+JgL(i)*hgL(i)) ...
		+ (JvBL(i)*hvBL(i)+JgBL(i)*hgBL(i)) ...
		- p*dVUllagedt*VvL(i)/VUllage ...
		- (mdotvL(i)*P.c_v*TvL(i) + mdotgL(i)*P.c_gv*TvL(i)) );
end
% compute Tdot for top vapor/gas layer
dTvBdt(P.n) = 1/(mvB(P.n)*P.c_v+mgB(P.n)*P.c_gv)*( ...
	QdotWVn ...
	+ (JvL(P.n-1)*hvL(P.n-1)+JgL(P.n-1)*hgL(P.n-1)) ...
	- (JvB(P.n)*hvB(P.n)+JgB(P.n)*hgB(P.n)) ...
	- (Jve*P.c_v*TvB(P.n)+Jge*P.c_gv*TvB(P.n)) ...
	- p*dVUllagedt*VvB(P.n)/VUllage ...
	- (mdotvB(P.n)*P.c_v*TvB(P.n) + mdotgB(P.n)*P.c_gv*TvB(P.n)) );

% wall masses
mwL = P.rhow*2*pi*((P.R+P.tw)^2-P.R^2)*hL + P.A*P.tw;
mwv = P.rhow*2*pi*((P.R+P.tw)^2-P.R^2)*(P.H-hL) + P.A*P.tw;
% height derivatives
dhLdt = JLTotal/P.A/P.rhoL;
dhvdt = -JLTotal/P.A/P.rhoL;
% wall mass derivatives
dmwLdt = P.rhow*2*pi*((P.R+P.tw)^2-P.R^2)*dhLdt;
dmwvdt = P.rhow*2*pi*((P.R+P.tw)^2-P.R^2)*dhvdt; % this is not reallt dmwvdt but rather -dmwLdt
% "enthalpy" terms
hmw = P.cw*Twv*dmwLdt*(dmwLdt>0) + P.cw*TwL*dmwLdt*(dmwLdt<0);
hvw = P.cw*TwL*dmwvdt*(dmwvdt>0) + P.cw*Twv*dmwvdt*(dmwvdt<0);
% compute QdotEW terms
QdotEWL = QdotEW*mwL/(mwL+mwv);
QdotEWv = QdotEW*mwv/(mwv+mwL);
% compute Tdot for wall
dTwLdt = 1/mwL*(1/P.cw*(QdotEWL-QdotWLTotal+hmw) - TwL*dmwLdt);
dTwvdt = 1/mwv*(1/P.cw*(QdotEWv-QdotWVTotal+hvw) - Twv*dmwvdt);

% liquid optimal grid temperatures (layers contained within TLB(P.k))
for i=1:P.nL
	if i==1
		TLIim1 = Ts;
	else
		TLIim1 = TLI(i-1);
	end
	if i==P.nL
		TLIip1 = TLB(P.k);
	else
		TLIip1 = TLI(i+1);
	end
	dTLIdt(i) = P.kappaL/l_L(i)/P.c_L/P.rhoL*((TLIip1-TLI(i))/l12_L(i+1)-(TLI(i)-TLIim1)/l12_L(i));
end

% vapor/gas optimal grid temperatures (layers contained within TvB(1))
for i=1:P.nV
	if i==1
		TvIim1 = Ts;
	else
		TvIim1 = TvI(i-1);
	end
	if i==P.nV
		TvIip1 =  TvB(1);
	else
		TvIip1 = TvI(i+1);
	end
	rhoi = pv/P.R_v/TvI(i) + pg/P.R_g/TvI(i);
	dTvIdt(i) = P.kappa_v/l_V(i)/P.c_v/rhoi*((TvIip1-TvI(i))/l12_V(i+1)-(TvI(i)-TvIim1)/l12_V(i));
end

GrvG = P.g*1/TvG*(Twv-TvG)*P.R^3/nuv^2;
RavG = abs(GrvG)*Prv;
alphavG = P.kappa_v/P.R*(0.15*(RavG)^(1/3));
QdotWVG = (pi*P.R^2+2*pi*P.R*(P.H-hL))*alphavG*(Twv-TvG);
dTvGdt = 1/(mv*P.c_v+mg*P.c_gv)*( ...
	- QdotVS ...
	+ QdotWVG ...
	- Jcd*P.c_p*Ts ...
	- p*dVUllagedt ...
	- (Jve*P.c_v*TvG+Jge*P.c_gv*TvG) ...
	- (mdotv*P.c_v*TvG + mdotg*P.c_gv*TvG));

% state derivatives
dxdt(1) = JLTotal;
dxdt(2:1+P.k) = dTLBdt;
dxdt(2+P.k:1+2*P.k-2) = dTLLdt;
dxdt(2+2*P.k-2:1+2*P.k-2+P.nL) = dTLIdt;
dxdt(2+2*P.k-2+P.nL) = mdotv;
dxdt(3+2*P.k-2+P.nL) = mdotg;
dxdt(4+2*P.k-2+P.nL:3+2*P.k-2+P.nL+P.n) = dTvBdt;
dxdt(4+2*P.k-2+P.nL+P.n:3+2*P.k-2+P.nL+2*P.n-2) = dTvLdt;
dxdt(4+2*P.k-2+P.nL+2*P.n-2:3+2*P.k-2+P.nL+2*P.n-2+P.nV) = dTvIdt;
dxdt(4+2*P.k-2+P.nL+2*P.n-2+P.nV) = dTsdt;
dxdt(5+2*P.k-2+P.nL+2*P.n-2+P.nV) = dTwLdt;
dxdt(6+2*P.k-2+P.nL+2*P.n-2+P.nV) = dTwvdt;
dxdt(7+2*P.k-2+P.nL+2*P.n-2+P.nV) = dTvGdt;

if sum(~isreal(dxdt))>0
	error('imaginary results at %g',t);
end

% update P (for vent state)
evalin('base',['LH2Model.VentState = ' num2str(P.VentState) ';']);

% must return a column vector
dxdt = dxdt';


function state = getVentState(P,p)

if P.Venting
	if p < P.pVentLow
		% turn off valve
		state = 0;
	elseif p > P.pVentHigh
		% turn on valve
		state = 1;
	else
		% stays at same value
		state = P.VentState;
	end
else
	state= 0;
end

function mdot = gasFlow(CA,gamma,rho,P1,P2)
% choked/nonchoked flow. 
if P1<P2
	mdot = -gasFlow(CA,gamma,rho,P2,P1);
else
	%assumes P1 always >= P2
	threshold = ((gamma+1)/2)^(gamma/(gamma-1));
	if P1/P2 >= threshold
		% choked flow
		mdot = CA*sqrt(gamma*rho*P1*(2/(gamma+1))^((gamma+1)/(gamma-1)));
	else
		% nonchoked
		mdot = CA*sqrt(2*rho*P1*(gamma/(gamma-1))*((P2/P1)^(2/gamma)-(P2/P1)^((gamma+1)/gamma)));
	end
end