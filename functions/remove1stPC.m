function [denoisedData] = remove1stPC(data)
%REMOVE1STPC reconstructs data without the first PC.
%   Denoising with principal component analysis 
% 
% EHS::20240410

% author: Elliot H Smith - https://github.com/elliothsmith/seizureCodes

% assumes data matrix is input as channels by samples. 
sz = size(data);
if sz(1)>sz(2); data = data'; end

% PCA
[w,pc,ev] = pca(data);

% reconstruction without 1st component
denoisedData = (pc(:,2:end)*w(:,2:end)');

end

