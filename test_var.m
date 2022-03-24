function [outputArg1,outputArg2] = test_var(varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if rem(nargin,2) ~= 0
    myTitle = varargin{nargin};
    numPlotInputs = nargin - 1;
else
    myTitle = 'Default Title';
    numPlotInputs = nargin;
end

plot(varargin{1:numPlotInputs})
title(myTitle)
end