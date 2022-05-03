function y = smooth(x , N )

% y = smooth(x , N )
%
% Smooth x with a window of N points... 
% Returned vector is same length as x

%N should be odd
N = round(N/2)*2 + 1;
W = ones(N,1)/N;
y = conv(x,W);

y = y((1:length(x))+(N-1)/2);