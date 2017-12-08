function drawTank(P,data,i,options)
% plots the tank given simulation data for a certain time point
% P is parameters
% data is data
% i is index (for time) in data
% options is a struct with options
%
% Copyright (c) 2017 United States Government as represented by
%     the Administrator of the National Aeronautics and Space Administration.
%     All Rights Reserved.
%
if nargin<4
	options.colors = 'rgb';
	options.t = 1;
	options.colorbar = 1;
end

% obtain state variables
t = data.t(i);
mL = data.mL(i,:);
TLB = data.TLB(i,:);
TLL = data.TLL(i,:)';
TLI = data.TLI(i,:)';
mv = data.mv(i,:);
mg = data.mg(i,:);
TvB = data.TvB(i,:)';
TvL = data.TvL(i,:)';
TvI = data.TvI(i,:)';
Ts = data.Ts(i,:);
TwL = data.TwL(i,:);
Twv = data.Twv(i,:);
TvG = data.TvG(i,:);

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
nuv = P.mu/rho;					% note, this is not Nusselt number (Nu)
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
VvB(1) = VUllage/P.k;
for i=2:P.n-1
	% bulk layer
	VvB(i) = pi*(P.R-deltav(i))^2*hvi;
	% boundary layer
	VvL(i) = (pi*P.R^2-pi*(P.R-deltav(i))^2)*hvi;
end
VvB(P.k) = VUllage/P.k;
VvL(P.k) = 0; % doing this only to make vectors same length for vector operations below, will ignore this value)

% do the plotting
hold on;

tw = 0.5;

colorMapSize = 256; 10000; % if use more than 256 then colorbar doesn't show colors correctly
if strcmp('rgb',options.colors)
	spectrum = jet(colorMapSize);
else
	spectrum = gray(colorMapSize);
end
minTemp = min(-0.02+[min(data.TLB) min(data.TLL) min(data.TvB) min(data.TvL) min(data.Ts) min(data.TwL) min(data.Twv)]);
maxTemp = max(0.02+[max(data.TLB) max(data.TLL) max(data.TvB) max(data.TvL) max(data.Ts) max(data.TwL) max(data.Twv)]);%max(21,max(data.Tw));
colorGain = colorMapSize/(maxTemp-minTemp);
colormap(spectrum); % set gcf's colormap to mine

% color in tank walls
% in contact with liquid
C = spectrum(max(1,round((TwL-minTemp)*colorGain)),:);
patch([-P.R-tw P.R+tw P.R+tw -P.R-tw],[0-tw 0-tw hL hL],C,'EdgeColor',C);
% in contact with vapor
C = spectrum(max(1,round((Twv-minTemp)*colorGain)),:);
patch([-P.R-tw P.R+tw P.R+tw -P.R-tw],[hL hL P.H+tw P.H+tw],C,'EdgeColor',C);

% set bottom and top deltas to 0 as they should be
deltaL(1) = 0;
deltaL(end) = 0;
deltav(1) = 0;
deltav(end) = 0;

% color in the bulk layers
for i=1:length(TLB)
	C = spectrum(max(1,round((TLB(i)-minTemp)*colorGain)),:);
	patch([-P.R+deltaL(i) P.R-deltaL(i) P.R-deltaL(i) -P.R+deltaL(i)],[hLi*(i-1) hLi*(i-1) hLi*(i) hLi*(i)],C,'EdgeColor',C);
end
for i=1:length(TvB)
	C = spectrum(max(1,round((TvB(i)-minTemp)*colorGain)),:);
	patch([-P.R+deltav(i) P.R-deltav(i) P.R-deltav(i) -P.R+deltav(i)],[hL+hvi*(i-1) hL+hvi*(i-1) hL+hvi*(i) hL+hvi*(i)],C,'EdgeColor',C);
end

% color in the boundary layers
for i=1:length(TLB)
	C = spectrum(max(1,round((TLL(i)-minTemp)*colorGain)),:);
	patch([-P.R -P.R+deltaL(i) -P.R+deltaL(i) -P.R],[hLi*(i-1) hLi*(i-1) hLi*(i) hLi*(i)],C,'EdgeColor',C);
	patch([P.R P.R-deltaL(i) P.R-deltaL(i) P.R],[hLi*(i-1) hLi*(i-1) hLi*(i) hLi*(i)],C,'EdgeColor',C);
end
for i=1:length(TvB)
	C = spectrum(max(1,round((TvL(i)-minTemp)*colorGain)),:);
	patch([-P.R -P.R+deltav(i) -P.R+deltav(i) -P.R],[hL+hvi*(i-1) hL+hvi*(i-1) hL+hvi*(i) hL+hvi*(i)],C,'EdgeColor',C);
	patch([P.R P.R-deltav(i) P.R-deltav(i) P.R],[hL+hvi*(i-1) hL+hvi*(i-1) hL+hvi*(i) hL+hvi*(i)],C,'EdgeColor',C);
end

% draw the boundaries for the bulk layers
for i=1:length(TLB)-1
	plot([-P.R P.R],[hLi*i hLi*i],'k:');
end
for i=1:length(TvB)-1
	plot([-P.R P.R],[hL+hvi*i hL+hvi*i],'k:');
end

% draw the boundaries for the boundary layers
for i=1:length(TLL)-2
	plot([-P.R+deltaL(i) -P.R+deltaL(i)],[hLi*i hLi*(i+1)],'k:');
	plot([P.R-deltaL(i) P.R-deltaL(i)],[hLi*i hLi*(i+1)],'k:');
end
for i=1:length(TvL)-2
	plot([-P.R+deltav(i) -P.R+deltav(i)],[hL+hvi*i hL+hvi*(i+1)],'k:');
	plot([P.R-deltav(i) P.R-deltav(i)],[hL+hvi*i hL+hvi*(i+1)],'k:');
end

% draw interface
plot([-P.R P.R],[hL hL],'k:');
C = spectrum(max(1,round((Ts-minTemp)*colorGain)),:);
plot([-P.R P.R],[hL hL],'Color',C,'LineWidth',3,'LineStyle','--');

% draw tank (inner wall)
line([-P.R P.R],[0 0],'Color','k');
line([-P.R P.R],[P.H P.H],'Color','k');
line([-P.R -P.R],[0 P.H],'Color','k');
line([P.R P.R],[0 P.H],'Color','k');
% draw tank (outer wall)
line([-P.R-tw P.R+tw],[0-tw 0-tw],'Color','k');
line([-P.R-tw P.R+tw],[P.H+tw P.H+tw],'Color','k');
line([-P.R-tw -P.R-tw],[0-tw P.H+tw],'Color','k');
line([P.R+tw P.R+tw],[0-tw P.H+tw],'Color','k');


xlim([-10 10]);
ylim([-tw-0.01 30+tw]);

if options.colorbar
	h = colorbar;%('West');
	labels = minTemp:((maxTemp-minTemp)/5):maxTemp;
	labels = round(labels*100)/100;
	set(h,'YTickLabel',labels);
end
%set(h,'YAxisLocation','left');
%set(h,'Position',[0.27 0.1377 0.0476 0.7556]);

set(gca,'XTick',[],'XColor','w');
set(gca,'YTick',[],'YColor','w','ZColor','w');
set(gcf,'Color','w');

if options.t
	h = text(-2.2*P.R,0,sprintf('%6.2f s',t));
	set(h,'FontName','Courier');
end

hold off;
