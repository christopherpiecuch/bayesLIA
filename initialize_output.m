% initialize_output
U = nan(NN_thin,N);         % Vector of regional VLM rate
C=nan(NN_thin,K);       % dice rolls
PI=nan(NN_thin,K);      % dice probs
V = nan(NN_thin,N);         % Total VLM vector
ALPHA = nan(NN_thin,1);     % Mean non-GIA VLM rate
OMEGA_2 = nan(NN_thin,1);   % Partial sill of regional VLM rate
RHO = nan(NN_thin,1);       % non-GIA VLM rate inverse range
EPSILON_2 = nan(NN_thin,1); % Variance in tide gauge data biases
