function FilterBankCSP(data,label,fs,varargin)
% FBCSP 
% Input
%   data:EEG data Chan*Sample*Trial/Chan*Sample(if so, more input required) 
%   label:labels 1*Trial
%   fs:sampling frequency
%   timelen:(optional)time range of each trial, default is [0.5 2]
%   feaNum:(optional)number of features to train,default is 4
%   freqband:(optional)band of frequency,defualt is 8:4:32
%   posmrk:mark of position of each trial

% Usage:
%   FilterBankCSP(data,label,varargin)/ FilterBankCSP(data,label,fs,posmrk,varagin)
%% Check Input
nchan = size(data,1);
nsample = size(data,2);
ntrial = length(label);

if size(data,3)==1
    if nargin<4
        error('position mark needed')
    end
    posmrk = varargin{1};
    if length(posmrk)~=ntrial
        error('different number of trials and position mark')
    end
    if nargin<5
        timelen = [0.5 2];
    else
        timelen = varargin{2};
    end
    if nargin<6
        feaNum = 4;
    else
        feaNum = varargin{3};
    end
    if nargin<7
        flist = transpose(8:4:32);
        freqband = [flist(1:end-1) flist(2:end)];
    else
        freqband = varargin{4};
    end
    
else
    if size(data,3)~=ntrial
        error('number of trials does not match number of labels')
    end
    if nargin<4
        timelen = [0.5 2];
    else
        timelen = varargin{1};
    end
    if nargin<5
        feaNum = 4;
    else
        feaNum = varargin{2};
    end
    if nargin<6
        flist = transpose(8:4:32);
        freqband = [flist(1:end-1) flist(2:end)];
    else
        freqband = varargin{3};
    end
end

%% Filter Bands



%% CSP Features


%% Feature Selection