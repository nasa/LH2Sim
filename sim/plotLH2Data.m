function plotLH2Data(data)
% plotLH2Data(data)
%	Plots results from data.
%	'data' is a data structure returned from LH2Simulate.
%
% Copyright (c) 2017 United States Government as represented by
%     the Administrator of the National Aeronautics and Space Administration.
%     All Rights Reserved.
%

close all;

specsB = {'k-' 'k--' 'k-.' 'k:' 'k-' 'k--' 'k-.' 'k:' 'k-' 'k--' 'k-.' 'k:'};
specsL = {'c-' 'c--' 'c-.' 'c:' 'c-' 'c--' 'c-.' 'c:' 'c-' 'c--' 'c-.' 'c:'};
specsS = {'m-' 'm-' 'm-' 'm-'};

plotInterface = false;

% plot masses
figure;
set(gcf,'position',get(gcf,'position')+[0 -300 0 300]);
subplot(3,1,1);
plot(data.t,data.mL,specsB{1});
title('Liquid Mass');
xlabel('Time (s)');
ylabel('kg');
axis tight;
subplot(3,1,2);
for i=1:size(data.mv,2)
	plot(data.t,data.mv(:,i),specsB{i});
end
title('Vapor Mass');
xlabel('Time (s)');
ylabel('kg');
axis tight;
subplot(3,1,3);
for i=1:size(data.mg,2)
	plot(data.t,data.mg(:,i),specsB{i});
end
title('Gas Mass');
xlabel('Time (s)');
ylabel('kg');
axis tight;

% plot temperatures
figure;
set(gcf,'position',get(gcf,'position')+[0 -300 0 300]);
%subplot(3,1,1);
subplot(2,1,1);
hold on;
handles = [];
for i=1:size(data.TLB,2)
	h = plot(data.t,data.TLB(:,i),specsB{i});
	if i==1, handles(1) = h; end
end
for i=1:size(data.TLL,2)
	h = plot(data.t,data.TLL(:,i),specsL{i});
	if i==1, handles(2) = h; end
end
handles(3) = plot(data.t,data.TwL,'r--','LineWidth',2);
handles(4) = plot(data.t,data.Twv,'r-.','LineWidth',2);
handles(5) = plot(data.t,data.Ts,'m--','LineWidth',2);
if plotInterface
	for i=1:size(data.TLI,2)
		plot(data.t,data.TLI(:,i),specsS{i});
	end
end
hold off;
title('Liquid Temperature');
xlabel('Time (s)');
ylabel('K');
legend(handles,'T_{i,B}','T_{i,L}','T_{wL}','T_{wv}','T_f');
axis tight;

% subplot(3,1,2);
subplot(2,1,2);
hold on;
handles = [];
for i=1:size(data.TvB,2)
	h = plot(data.t,data.TvB(:,i),specsB{i});
	if i==1, handles(1) = h; end
	%text(data.t(end),data.TvB(end,i),['\leftarrow T_{' num2str(i) ',B}']);
end
for i=1:size(data.TvL,2)
	h = plot(data.t,data.TvL(:,i),specsL{i});
	if i==1, handles(2) = h; end
	%text(data.t(end),data.TvL(end,i),['\leftarrow T_{' num2str(i+1) ',L}']);
end
handles(3) = plot(data.t,data.TwL,'r--','LineWidth',2);
handles(4) = plot(data.t,data.Twv,'r-.','LineWidth',2);
handles(5) = plot(data.t,data.Ts,'m--','LineWidth',2);
if plotInterface
	for i=1:size(data.TvI,2)
		plot(data.t,data.TvI(:,i),specsS{i});
	end
end
handles(6) = plot(data.t,data.TvG,'g--','LineWidth',2);
hold off;
title('Vapor/gas Temperature');
xlabel('Time (s)');
ylabel('K');
legend(handles,'T_{i,B}','T_{i,L}','T_{wL}','T_{wv}','T_f','Average T_v');
axis tight;

% closer plot of liquid temps
figure;
hold on;
handles = [];
for i=1:size(data.TLB,2)
	h = plot(data.t,data.TLB(:,i),specsB{i});
	if i==1, handles(1) = h; end
end
for i=1:size(data.TLL,2)
	h = plot(data.t,data.TLL(:,i),specsL{i});
	if i==1, handles(2) = h; end
end
title('Liquid Temperature');
xlabel('Time (s)');
ylabel('K');
legend(handles,'T_{i,B}','T_{i,L}');
axis tight;
hold off;

% subplot(3,1,3);
% plot(data.t,data.Tw,specs{1});
% title('wall temperature');
% xlabel('Time (s)');
% ylabel('K');
% axis tight;
