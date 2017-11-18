% Run nominal loading scenario

% 1. initialize parameters
LH2ModelParams;
LH2Model.n = 3;
LH2Model.k = 3;

% 2. run simulation with default timespan
nominal = LH2Simulate(0:1:5000);

% 3. plot results
plotLH2Data(nominal);
