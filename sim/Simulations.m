% LH2 simulations

% 1g simulation
if 1
	LH2ModelParams;
	data = LH2Simulate(0:1:4000);

	plotLH2Data(data);
	%testPlot;
	makeMovie;
end

% 7g simulation
if 0
	LH2ModelParams;
	LH2Model.g = 7*LH2Model.g;
	data = LH2Simulate(0:1:4000);

	plotLH2Data(data);
	makeMovie;
end

% 1e-5g simulation
if 0
	LH2ModelParams;
	LH2Model.g = 1e-5*LH2Model.g;
	data = LH2Simulate;

	plotLH2Data(data);
	makeMovie;
end