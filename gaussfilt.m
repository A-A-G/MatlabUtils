function [ yfilt ] = gaussfilt( xdata, ydata, sigma, xfilt, M )
%GAUSSFILT Non-Uniform 1D-Gauss Filter
%Window size: [t-M,t+M], t in xfilt 
    if nargin < 3
       error('gaussfilt: wrong parameter count') 
    elseif nargin == 3
        xfilt = xdata;
        M = 3 * sigma;
    elseif nargin == 4
        M = 3 * sigma;
    end
    yfilt = zeros(size(xfilt));
    for i = 1:length(xfilt)
        indices = (xdata > (xfilt(i)-M)) & (xdata < (xfilt(i)+M));
        gaussFactors = gaussKernel(xdata(indices)-xfilt(i), sigma);
        yfilt(i) = sum(gaussFactors .* ydata(indices)) ./ sum(gaussFactors);
    end
%     figure
%     plot(xdata,ydata);
%     hold on;
%     plot(xfilt,yfilt);
%     hold off;
end

function [ k ] = gaussKernel( x, sigma)
    k =  1 / (sigma * sqrt(2 * pi)) * exp(-0.5 * (x ./ sigma).^2);
end