function [bounds,bounds_t,env] = oscillatorSylSeg(envelopes,c,T,lag_correction)
%function [bounds,bounds_t,env] = oscillatorSylSeg(envelopes,c,T)
%
% Performs segmentation on a cell array of amplitude envelopes 
%
% Inputs:
%   envelopes:  cell array of envelopes to segment (at fs = 1000 Hz)
%   c:          damping coefficient of the oscillator
%   T:          period of the oscillator
%
% Outputs:
%   bounds:     segment boundaries in 1000 Hz frames
%   bounds_t:   segment boundaries in seconds
%   env:        oscillator amplitudes corresponding to input envelopes
%
% (c) Okko Rasanen, 2015
% okko.rasanen@aalto.fi

bounds = cell(length(envelopes),1);
bounds_t = cell(length(envelopes),1);

k = 1;   % Fix k = 1, define only mass

% Get mass
b = 2*pi/T;   
m = k/b^2;

% Calculate Q

Q = sqrt(m*k)/c;
fprintf('Oscillator Q-value: %0.4f, center frequency: %0.1f Hz, bandwidth: %0.1f Hz.\n',Q,1/T,1/T/Q);

env = cell(length(envelopes),1);

for signal = 1:length(envelopes)
        
    %e = envelopes{signal};
    e = envelopes;
        
    x = zeros(length(e),1);
    a = zeros(length(e),1);
    v = zeros(length(e),1);
    x(1) = 0;
    a(1) = 0;
        
    for t = 2:length(e)
        
        f_up = e(t);        
        f_down = -k*x(t-1)-c*v(t-1);      
        f_tot = f_up+f_down;
        
        % Get acceleration from force
        a(t) = f_tot./m;
        
        % Get velocity from acceleration
        v(t) = v(t-1)+a(t).*0.001;
        
        % Get position from velocity
        x(t) = x(t-1)+v(t).*0.001;
        
    end
    
    % Do phase shift correction by minimizing L2 norm between oscillator
    % waveform and signal envelope
    
    if nargin <4
        lagvals = 0:1:200;
        
        cor = zeros(length(lagvals),1);
        for j = 1:length(lagvals)
            cor(j) = norm(x-circshift(e,lagvals(j)));
        end
        
        [val,i] = min(cor);
        
        lag_correction = lagvals(i)+1;
    end
    
    
    if(lag_correction ~= 0)
    x = [x(lag_correction:end);zeros(lag_correction-1,1)];
    end    
    
    % Find minima
    [maxtab,mintab] = peakdet(x,0.00001.*(max(x)-min(x)));
        
    if(~isempty(mintab))        
        bounds{signal} = [find(x > 0,1);mintab(:,1)]; % Get boundaries
    else
       bounds{signal} = [1;length(x)]; 
    end
        
    bounds_t{signal} = bounds{signal}./1000;
    env{signal} = x;
end

