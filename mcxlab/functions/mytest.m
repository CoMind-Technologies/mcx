function [outputArg1,outputArg2] = mytest(inputArg1,varargin)
%MYTEST Summary of this function goes here
%   Detailed explanation goes here
x = [1 2 3 4
    4 3 2 1
    6 7 8 9];
outputArg1 = zeros(3,length(varargin));
outputArg1(:,:)=x(1,cell2mat(varargin));
outputArg2=length(varargin);
end
