function [SpectrumPwr, SpectrumMax, f_max, f_mean, f_median, f_std] = HRV_measures(Fp,FFTpwr,nx,Fmin,Fmax,COR,dim)
% function Output = HRV_measures(Fp,FFTpwr,nx,Fmin,Fmax,COR)

% This program calculates the Spectrum Power (SpectrumPwr), Maximum Spectrum
% (SpectrumMax) Value and the frequency of the Maximum (f_max), mean (f_mean)
% median (f_median) and standard deviation (f_std) % of the power Spectrum FFTpwr 
% distribution between frequency Fmin and Fmax for Samples defined by nx.
% 
% 
% Fp is a vector containing the Frequency
% FFTpwr is a matrices whose column contains Spectrum of different samples
% indiced by column index
% nx is a vector defining samples number for which spectrum measures are computed
% COR defines whether we want to adjust for spectrum power at Fmin and Fmax
% due to discrete nature of frequency values (Fmin/Fmax may falls between
% discrete frequency samples). Adjustment performed is a simple linear interpolation.
% 
% Output results format are row vectors of length (length(nx))
if nargin < 7
    dim = 1;
end
if nargin < 6
    COR = 1;
end

Ndim = ndims(FFTpwr);
order  = [1 2];

if dim == 2 & Ndim > 1
    order = [2 1];
    FFTpwr = permute(FFTpwr,order);
    Fp = permute(Fp,order);
end

if isempty(nx) & Ndim>1
    Nx = size(FFTpwr,Ndim);
    nx = 1:Nx;
elseif isempty(nx)
    nx = 1;
end
[Ni Nj] = size(Fp);
if Ni<Nj;   Fp = Fp'; end

NFFT = length(Fp); 
DFp =  Fp(2) - Fp(1);
% Calculate VLF values
% calculate the values for each sides
rg = find(Fp>=Fmin & Fp<=Fmax); rg1 = rg(1); rg2 = max(rg);
%left
dy1 = COR*(Fp(rg1)-Fmin)/DFp;
y1 = FFTpwr(rg1,nx) - dy1*( FFTpwr(rg1,nx) - FFTpwr(max(1,rg1-1),nx) );
x1 = Fp(rg1) - dy1*(Fp(rg1) - Fp(max(1,rg1-1) ));
%right
dy2 = COR*(Fmax - Fp(rg2) )/DFp;
y2 = FFTpwr(rg2,nx) - dy1*(FFTpwr(rg2,nx) - FFTpwr(min(NFFT,rg2+1),nx) );
x2 = Fp(rg2) - dy1*(Fp(rg2) - Fp(min(NFFT,rg2+1) ));

SpectrumPwr = zeros(1,length(nx));
SpectrumMax = zeros(1,length(nx));
f_max = zeros(1,length(nx));


SpectrumPwr(:) = trapz([x1;Fp(rg);x2],[y1;FFTpwr(rg,nx);y2],1); %in ms^2 
SpectrumMax = max(FFTpwr(rg,nx));
 

for i = 1:length(nx)
    [nmax nx2 dfg] = find( FFTpwr(rg,nx(i)) ==  max( FFTpwr(rg,nx(i)) ) );
    f_max(i) =  Fp( rg(min(nmax)) );
end

%----------- Didier rev 4.1 4/14/2010 ---------------
% Compute frequency related geometric measures (mean, median, std)
% of the FFTpwr distribution
if nargout > 1
    [f_mean, f_median, f_std] = geometric_measures(Fp(rg),FFTpwr(rg,nx),1);
    SpectrumMax = permute(SpectrumMax,order);
    f_mean = permute(f_mean,order);
    f_median = permute(f_median,order);
    f_std = permute(f_std,order);
    
    f_max = permute(f_max,order);
else
    SpectrumMax = [];
    f_mean = [];
    f_median = [];
    f_std = [];
    f_max = [];
end
%----------- Didier  ---------------

    
SpectrumPwr = permute(SpectrumPwr,order);

