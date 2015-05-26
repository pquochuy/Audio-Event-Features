function [wfeature] = extractTemkoFeatures(x,fs,frameSize,hop)
   % extract Temko features for 100ms window
   len = length(x);
	if ~exist('frameSize','var') || isempty(frameSize)
		wSize = 0.03; 	% 30ms window size
	end
	
	if ~exist('hop','var') || isempty(hop)
		hop = 0.015;	% 15ms overlap
	end
   
   frameSize = frameSize * fs;    % 20ms
   hop = hop * fs;             % 10ms
   
   %%%%%%%%%% Feature Extraction
   %consider preallocating for speed
   nFrame = size(1:hop:(len - frameSize + 1), 2);
   ff = zeros(nFrame,60);
   currentFrameIndex = 1;   % frame index
   previousFrame = [];
   for i = 1 : hop : (len - frameSize + 1) % length(x) - frameSize
       % Note that we add EPS to prevent divide by 0 errors 
       currentFrame = x(i : i + frameSize - 1) + eps;
       % Hamming windowing
       currentFrame = currentFrame .* hamming(length(currentFrame));
       
       % feature extraction begin here
       ffbe = fFFBE(currentFrame,fs,16);   % 16 frequency-filter log filter-bank energies
       ffbe_delta = deltas(ffbe')'; %note the transpose operator!
       ffbe_delta_delta = deltas(ffbe_delta')';
       ff(currentFrameIndex,1:16) = ffbe;
       ff(currentFrameIndex,17:32) = ffbe_delta;
       ff(currentFrameIndex,33:48) = ffbe_delta_delta;
       ff(currentFrameIndex,49) = fzcr(currentFrame);
       ff(currentFrameIndex,50) = fShortTimeEnergy(currentFrame);
       % log filter-bank energies
       [z] = melfilterbankenergies(currentFrame,fs,16);
       % 4 subband energies
       [sbe,~] = fSubbandEnergies2(z,4);
       ff(currentFrameIndex,51:54) = sbe;
       % 4 subband spectral flux
       nSubband = 4;
       if(isempty(previousFrame))
           ff(currentFrameIndex,55:58) = zeros(1,4)*0.0;
       else
           [prez] = melfilterbankenergies(previousFrame,fs,16);
           subLen = 16/nSubband;
           for k = 1:nSubband
               ff(currentFrameIndex,54+k) = fSpectralFlux2(z(((k-1)*subLen + 1) : k*subLen),...
                   prez(((k-1)*subLen + 1) : k*subLen));
           end
       end
       ff(currentFrameIndex,59) = fSpectralCentroid(currentFrame,fs);
       ff(currentFrameIndex,60) = fSpectralBandwidth(currentFrame,fs);
       
       previousFrame = currentFrame;
       currentFrameIndex = currentFrameIndex + 1;
   end     
   wfeature = zeros(1,120);
   wfeature(1:60) = mean(ff);
   wfeature(61:120) = std(ff);
end

