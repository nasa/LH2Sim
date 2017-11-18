function plotLH2(name,dataSets,gain)
% plotLH2(name,data,gain), or
% plotLH2(name,{data1,data2,...},gain)
%	Plots results from multiple data sets on the same axes.
%	name is the name of the signal to plot.
%	Each 'data' is a data structure returned from LH2Simulate that contains
%		a field given by name.
%	gain is an optional multiplication factor for units conversion

% set default gain
if nargin<3
	gain = 1;
end

if ~iscell(dataSets)
	dataSets = {dataSets};
end

args = {};
for i=1:length(dataSets)
	data = dataSets{i};
	args{end+1} = data.t;
	args{end+1} = data.(name)*gain;
end
plot(args{:})
xlim([0 data.t(end)]);
xlabel('Time (s)');
grid on;