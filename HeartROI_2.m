% Good Prototype this time, we have a good root
function [HS1, DC, H1] = HeartROI_2(I, A, B, C, D)

 HS1 = zeros(size(I));
 DC = zeros(A, B, C);
 H1 = zeros(A, B, C);
 for slice  = 1:C
     HS0 = zeros(A, B, 1, D);
     for time = 1:D
         % does Fourier transform frame by frame
        Magnitude = (fft2(I(:,:, slice, time)));
        % size(Magnitude) is 216 by 256
        HS0(:,:,1,time) = Magnitude;
        % size(HS0)is 216 by 256 by 1 by 30
     end
    HS1(:,:,slice,:) = abs(ifftn(HS0));
    % size (HS1)is 216 by 256 by 10 by 30
%     subplot(2,5,slice);
%     imshow (HS1(:,:,slice+10), [])
 end

 for temporalfreq = 1:C
     if temporalfreq == 1
         DC = HS1(:,:,:, temporalfreq);
     end
     if temporalfreq == 2
         H1 = HS1(:,:,:, temporalfreq);
     end
 end
 
% figure(2);
% for slice  = 1:C
%     subplot(2,6,slice);
%     imshow (H1(:,:,slice), [])
% end
