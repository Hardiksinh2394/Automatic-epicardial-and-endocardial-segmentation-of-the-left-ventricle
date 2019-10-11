function [endocardium] = Endocardium(BW, DC_AVG)

endocardium = zeros(size(BW));
epicard = zeros(size(BW));
mask = zeros(size(BW));
[A, B, C] = size(BW);

A2 = BW;
A3 = BW;

Centroid_vector = zeros(2,C);
% for iter = 1:5
for i = 1:C
binaryImage = true(size(BW(:,:,i)));
labeledImage = logical(binaryImage);
measurements = regionprops(labeledImage, BW(:,:,i), 'WeightedCentroid');
centerOfMass = measurements.WeightedCentroid;
Centroid_vector(:,i) = centerOfMass;
%     https://fr.mathworks.com/matlabcentral/answers/368901-calculate-center-of-image
%     C = round([A B]/2) ;
%     plot(C(1),C(2),'*r');
%     plot(stat(i).Centroid(1),stat(i).Centroid(2),'ro');
end
Centroid_vector = (Centroid_vector)';

for k = 1:C    
    m = 1:size(BW, 1);
    l = 1:size(BW, 2);   
    % maybe an if statement
    Dist_hist = sqrt((m.' - Centroid_vector(k,2)) .^ 2 + (l - Centroid_vector(k,1)) .^ 2);%./Storage1(y,x,k);
    mask(:,:,k) = Dist_hist;
    
end

for i = 1:A
    for j = 1:B
        for k = 1:C
            if (mask(i,j,k) > 10)
               mask(i,j,k) = 0;
            else 
               mask(i,j,k) = 1;
            end
            A2(i,j,k) = BW(i,j,k) * mask(i,j,k);
            %Storage_num5(i,j,k) = mask(i,j,k) * Storage_num4(i,j,k);
        end
    end
end
Centroid_vector = (Centroid_vector)';
% end
for i = 1:C
binaryImage = true(size(A2(:,:,i)));
labeledImage = logical(binaryImage);
measurements = regionprops(labeledImage, A2(:,:,i), 'WeightedCentroid');
centerOfMass = measurements.WeightedCentroid;
Centroid_vector(:,i) = centerOfMass;
%     https://fr.mathworks.com/matlabcentral/answers/368901-calculate-center-of-image
%     C = round([A B]/2) ;
%     plot(C(1),C(2),'*r');
%     plot(stat(i).Centroid(1),stat(i).Centroid(2),'ro');
end
Centroid_vector = (Centroid_vector)';

for k = 1:C    
    m = 1:size(A2, 1);
    l = 1:size(A2, 2);   
    % maybe an if statement
    Dist_hist = sqrt((m.' - Centroid_vector(k,2)) .^ 2 + (l - Centroid_vector(k,1)) .^ 2);%./Storage1(y,x,k);
    mask(:,:,k) = Dist_hist;
    
end

for i = 1:A
    for j = 1:B
        for k = 1:C
            if (mask(i,j,k) > 10)
               mask(i,j,k) = 0;
            else 
               mask(i,j,k) = 1;
            end
            A3(i,j,k) = A2(i,j,k) * mask(i,j,k);
            %Storage_num5(i,j,k) = mask(i,j,k) * Storage_num4(i,j,k);
        end
    end
end
Centroid_vector = (Centroid_vector)';
for i = 1:C
binaryImage = true(size(A3(:,:,i)));
labeledImage = logical(binaryImage);
measurements = regionprops(labeledImage, A3(:,:,i), 'WeightedCentroid');
centerOfMass = measurements.WeightedCentroid;
Centroid_vector(:,i) = centerOfMass;
%     https://fr.mathworks.com/matlabcentral/answers/368901-calculate-center-of-image
%     C = round([A B]/2) ;
%     plot(C(1),C(2),'*r');
%     plot(stat(i).Centroid(1),stat(i).Centroid(2),'ro');
end
Centroid_vector = (Centroid_vector)';
for k = 1:C    
    m = 1:size(A3, 1);
    l = 1:size(A3, 2);
    % maybe an if statement
    Dist_hist = sqrt((m.' - Centroid_vector(k,2)) .^ 2 + (l - Centroid_vector(k,1)) .^ 2);%./Storage1(y,x,k);
    mask(:,:,k) = Dist_hist;
    
end

for i = 1:A
    for j = 1:B
        for k = 1:C
            if (mask(i,j,k) > 45)
                mask(i,j,k) = 0;
            else 
                mask(i,j,k) = 1;
            end
            BW(i,j,k) = BW(i,j,k) * mask(i,j,k);
            %Storage_num5(i,j,k) = mask(i,j,k) * Storage_num4(i,j,k);
        end
    end
end 
se = strel('disk',5);
I = imopen(logical(BW), se);
for i = 1:C
    endocardium(:,:,i) = bwpropfilt(I(:,:,i),'FilledArea',1); % 'ConvexArea' 'FilledArea'    
%     if i  <=8
%         endocardium(:,:,i) = bwpropfilt(I(:,:,i),'EquivDiameter',[25 35]); % 'ConvexArea' 'FilledArea'
%     else
%         endocardium(:,:,i) = bwpropfilt(I(:,:,i),'EquivDiameter',[20 30]);
%     %endocardium(:,:,i) = imclose(BW(:,:,i), se);
%     end   
end




