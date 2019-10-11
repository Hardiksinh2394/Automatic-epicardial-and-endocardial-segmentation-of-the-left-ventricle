function [polarized] = Centroid_pol(H1_Extracted, DC_AVG)

polarized = zeros(size(DC_AVG));
[A, B, C] = size(DC_AVG);

Centroid_vector = zeros(2,C);
for i = 1:C
% some source code at https://fr.mathworks.com/matlabcentral/answers/28996-centroid-of-an-image
%     Ibw = Storage1(:,:,i);
%     stat = regionprops(Ibw,'centroid');
%  code actually used from https://fr.mathworks.com/matlabcentral/answers/350566-how-to-find-weighted-centroid-of-an-entire-image-in-matlab
% Tried with  centroid and not weighted center, unuseful in our case
binaryImage = true(size(H1_Extracted(:,:,i)));
labeledImage = logical(binaryImage);
measurements = regionprops(labeledImage, H1_Extracted(:,:,i), 'WeightedCentroid');
centerOfMass = measurements.WeightedCentroid;
Centroid_vector(:,i) = centerOfMass;
%     https://fr.mathworks.com/matlabcentral/answers/368901-calculate-center-of-image
%     C = round([A B]/2) ;
%     plot(C(1),C(2),'*r');
%     plot(stat(i).Centroid(1),stat(i).Centroid(2),'ro');
end
Centroid_vector = (Centroid_vector)';

for i = 1:C
    polarized(:,:,i) = ImToPolar(DC_AVG(:,:,i), 0, 1, A, B, Centroid_vector(i,:));
end
