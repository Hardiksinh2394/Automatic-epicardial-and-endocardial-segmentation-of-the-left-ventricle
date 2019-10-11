function [ROI_slices, ROI_H1, Storage_num6, mask] = Heart_detector_3(Img_set, Average)

ROI_slices = zeros(size(Img_set));
ROI_H1 = zeros(size(Img_set));
Storage1 = zeros(size(Img_set));


[A1, B1, C] = size(Img_set);

for i = 1:C
    smoothed = medfilt2(Img_set(:,:,i));
    Storage1(:,:,i) = smoothed;
end


maximum_of3D = max(Storage1(:));
for i = 1:A1
    for j = 1:B1
        for k  = 1:C
        if Storage1(i,j,k) < 0.1 * maximum_of3D
                Storage1(i,j,k) = 0;
        end
        end
    end
end
    

% for slice  = 1:C
%     subplot(2,5,slice);
%     imshow (Storage1(:,:,slice), [])
% end

Centroid_vector = zeros(2,C);

for i = 1:C
    % some source code at https://fr.mathworks.com/matlabcentral/answers/28996-centroid-of-an-image
    %     Ibw = Storage1(:,:,i);
    %     stat = regionprops(Ibw,'centroid');
    %  code actually used from https://fr.mathworks.com/matlabcentral/answers/350566-how-to-find-weighted-centroid-of-an-entire-image-in-matlab
    % Tried with  centroid and not weighted center, unuseful in our case
    binaryImage = true(size(Storage1(:,:,i)));
    labeledImage = logical(binaryImage);
    measurements = regionprops(labeledImage, Storage1(:,:,i), 'WeightedCentroid');
    centerOfMass = measurements.WeightedCentroid;
    Centroid_vector(:,i) = centerOfMass;
    %     https://fr.mathworks.com/matlabcentral/answers/368901-calculate-center-of-image
    %     C = round([A B]/2) ;
    %     plot(C(1),C(2),'*r');
    %     plot(stat(i).Centroid(1),stat(i).Centroid(2),'ro');
end
Centroid_vector = (Centroid_vector)';

% Fitting part
% https://fr.mathworks.com/matlabcentral/answers/405461-how-can-i-create-a-linear-line-of-best-fit-for-a-set-of-3-dimensional-coordinates
% code found in the link uses svd and has been adapted to our case
x = Centroid_vector(:,1);
y = Centroid_vector(:,2);
z = (1:C)';
% May comment this implementation to put process in a while loop
xyz=[x(:),y(:),z(:)];
r = mean(xyz);
xyz=bsxfun(@minus,xyz,r);
[~,~,V]=svd(xyz,0);
x_fit=@(z_fit) r(1)+(z_fit-r(3))/V(3,1)*V(1,1);
y_fit=@(z_fit) r(2)+(z_fit-r(3))/V(3,1)*V(2,1);
A = x_fit(z);
B = y_fit(z);
MA = mean(A);
MB = mean(B);
LLS_centroids = [A(:),B(:)];
%Centro_2D = [x, y];

% figure(5);
% for i = 1:C
%     subplot(2,6,i);
%     imshow(Storage1(:,:,i),[]);
%     hold on;
%     plot(Centroid_vector(i,1), Centroid_vector(i,2), 'r+', 'LineWidth', 1, 'MarkerSize', 16);
%     hold on;
%     plot(LLS_centroids(i,1), LLS_centroids(i,2), 'g+', 'LineWidth', 1, 'MarkerSize', 16);
% %  mean the centroids but what the goal to find the central 3D one
% %  now we have to compute the pixels distances distribution weighted by
% %  their values
% end

Dist_hist = 0;
mask = zeros(A1, B1, C);
% https://fr.mathworks.com/matlabcentral/answers/380796-how-to-find-distance-of-each-pixel-of-image-from-center-of-image

for k = 1:C    
    m = 1:size(Storage1, 1);
    l = 1:size(Storage1, 2);
    c(k,:) = LLS_centroids(k,:);
    % maybe an if statement
    Dist_hist = sqrt((m.' - Centroid_vector(k,2)) .^ 2 + (l - Centroid_vector(k,1)) .^ 2);%./Storage1(y,x,k);
    mask(:,:,k) = Dist_hist;
    
end

Storage_num2 = Storage1;
Storage_num3 = Storage1;
Storage_num4 = Storage1;
Storage_num5 = Storage1;
Storage_num6 = Storage1;

%ROI_image_preview = zeros(size(Img_set));
for i = 1:A1
    for j = 1:B1
        for k = 1:C
            if (mask(i,j,k) > 70)
                mask(i,j,k) = 0;
            else 
                mask(i,j,k) = 1;
            end
            %ROI_slices(i,j,k) = Average(i,j,k) * mask(i,j,k);
            Storage1(i,j,k) = mask(i,j,k) * Storage1(i,j,k);
        end
    end
end 
Centroid_vector = (Centroid_vector)';

for i = 1:C
% some source code at https://fr.mathworks.com/matlabcentral/answers/28996-centroid-of-an-image
%     Ibw = Storage1(:,:,i);
%     stat = regionprops(Ibw,'centroid');
%  code actually used from https://fr.mathworks.com/matlabcentral/answers/350566-how-to-find-weighted-centroid-of-an-entire-image-in-matlab
% Tried with  centroid and not weighted center, unuseful in our case
binaryImage = true(size(Storage1(:,:,i)));
labeledImage = logical(binaryImage);
measurements = regionprops(labeledImage, Storage1(:,:,i), 'WeightedCentroid');
centerOfMass = measurements.WeightedCentroid;
Centroid_vector(:,i) = centerOfMass;
%     https://fr.mathworks.com/matlabcentral/answers/368901-calculate-center-of-image
%     C = round([A B]/2) ;
%     plot(C(1),C(2),'*r');
%     plot(stat(i).Centroid(1),stat(i).Centroid(2),'ro');
end
Centroid_vector = (Centroid_vector)';

for k = 1:C    
    m = 1:size(Storage1, 1);
    l = 1:size(Storage1, 2);
    c(k,:) = LLS_centroids(k,:);
    % maybe an if statement
    Dist_hist = sqrt((m.' - Centroid_vector(k,2)) .^ 2 + (l - Centroid_vector(k,1)) .^ 2);%./Storage1(y,x,k);
    mask(:,:,k) = Dist_hist;
    
end

for i = 1:A1
    for j = 1:B1
        for k = 1:C
            if (mask(i,j,k) > 70)
                mask(i,j,k) = 0;
            else 
                mask(i,j,k) = 1;
            end
            ROI_slices(i,j,k) = Average(i,j,k) * mask(i,j,k);
            Storage_num2(i,j,k) = mask(i,j,k) * Storage1(i,j,k);
        end
    end
end 

Centroid_vector = (Centroid_vector)';
for i = 1:C
binaryImage = true(size(Storage_num2(:,:,i)));
labeledImage = logical(binaryImage);
measurements = regionprops(labeledImage, Storage_num2(:,:,i), 'WeightedCentroid');
centerOfMass = measurements.WeightedCentroid;
Centroid_vector(:,i) = centerOfMass;
%     https://fr.mathworks.com/matlabcentral/answers/368901-calculate-center-of-image
%     C = round([A B]/2) ;
%     plot(C(1),C(2),'*r');
%     plot(stat(i).Centroid(1),stat(i).Centroid(2),'ro');
end
Centroid_vector = (Centroid_vector)';

for k = 1:C    
    m = 1:size(Storage_num2, 1);
    l = 1:size(Storage_num2, 2);
    c(k,:) = LLS_centroids(k,:);
    % maybe an if statement
    Dist_hist = sqrt((m.' - Centroid_vector(k,2)) .^ 2 + (l - Centroid_vector(k,1)) .^ 2);%./Storage1(y,x,k);
    mask(:,:,k) = Dist_hist;
    
end

for i = 1:A1
    for j = 1:B1
        for k = 1:C
            if (mask(i,j,k) > 70)
                mask(i,j,k) = 0;
            else 
                mask(i,j,k) = 1;
            end
            ROI_slices(i,j,k) = Average(i,j,k) * mask(i,j,k);
            %ROI_H1(i,j,k) = mask(i,j,k) * Storage_num3(i,j,k);
            Storage_num3(i,j,k) = mask(i,j,k) * Storage1(i,j,k);
        end
    end
end 

Centroid_vector = (Centroid_vector)';
for i = 1:C
binaryImage = true(size(Storage_num3(:,:,i)));
labeledImage = logical(binaryImage);
measurements = regionprops(labeledImage, Storage_num3(:,:,i), 'WeightedCentroid');
centerOfMass = measurements.WeightedCentroid;
Centroid_vector(:,i) = centerOfMass;
%     https://fr.mathworks.com/matlabcentral/answers/368901-calculate-center-of-image
%     C = round([A B]/2) ;
%     plot(C(1),C(2),'*r');
%     plot(stat(i).Centroid(1),stat(i).Centroid(2),'ro');
end
Centroid_vector = (Centroid_vector)';

for k = 1:C    
    m = 1:size(Storage_num3, 1);
    l = 1:size(Storage_num3, 2);
    c(k,:) = LLS_centroids(k,:);
    % maybe an if statement
    Dist_hist = sqrt((m.' - Centroid_vector(k,2)) .^ 2 + (l - Centroid_vector(k,1)) .^ 2);%./Storage1(y,x,k);
    mask(:,:,k) = Dist_hist;
    
end

for i = 1:A1
    for j = 1:B1
        for k = 1:C
            if (mask(i,j,k) > 70)
                mask(i,j,k) = 0;
            else 
                mask(i,j,k) = 1;
            end
            ROI_slices(i,j,k) = Average(i,j,k) * mask(i,j,k);
            %ROI_H1(i,j,k) = mask(i,j,k) * Storage_num3(i,j,k);
            Storage_num4(i,j,k) = mask(i,j,k) * Storage1(i,j,k);
        end
    end
end 

Centroid_vector = (Centroid_vector)';
for i = 1:C
binaryImage = true(size(Storage_num4(:,:,i)));
labeledImage = logical(binaryImage);
measurements = regionprops(labeledImage, Storage_num4(:,:,i), 'WeightedCentroid');
centerOfMass = measurements.WeightedCentroid;
Centroid_vector(:,i) = centerOfMass;
%     https://fr.mathworks.com/matlabcentral/answers/368901-calculate-center-of-image
%     C = round([A B]/2) ;
%     plot(C(1),C(2),'*r');
%     plot(stat(i).Centroid(1),stat(i).Centroid(2),'ro');
end
Centroid_vector = (Centroid_vector)';

for k = 1:C    
    m = 1:size(Storage_num4, 1);
    l = 1:size(Storage_num4, 2);
    c(k,:) = LLS_centroids(k,:);
    % maybe an if statement
    Dist_hist = sqrt((m.' - Centroid_vector(k,2)) .^ 2 + (l - Centroid_vector(k,1)) .^ 2);%./Storage1(y,x,k);
    mask(:,:,k) = Dist_hist;
    
end

for i = 1:A1
    for j = 1:B1
        for k = 1:C
            if (mask(i,j,k) > 70)
                mask(i,j,k) = 0;
            else 
                mask(i,j,k) = 1;
            end
            ROI_slices(i,j,k) = Average(i,j,k) * mask(i,j,k);
            Storage_num5(i,j,k) = mask(i,j,k) * Storage1(i,j,k);
            %Storage_num5(i,j,k) = mask(i,j,k) * Storage_num4(i,j,k);
        end
    end
end 

Centroid_vector = (Centroid_vector)';
for i = 1:C
binaryImage = true(size(Storage_num5(:,:,i)));
labeledImage = logical(binaryImage);
measurements = regionprops(labeledImage, Storage_num5(:,:,i), 'WeightedCentroid');
centerOfMass = measurements.WeightedCentroid;
Centroid_vector(:,i) = centerOfMass;

%     https://fr.mathworks.com/matlabcentral/answers/368901-calculate-center-of-image
%     C = round([A B]/2) ;
%     plot(C(1),C(2),'*r');
%     plot(stat(i).Centroid(1),stat(i).Centroid(2),'ro');
end
Centroid_vector = (Centroid_vector)';

for k = 1:C    
    m = 1:size(Storage_num5, 1);
    l = 1:size(Storage_num5, 2);
    c(k,:) = LLS_centroids(k,:);
    % maybe an if statement
    Dist_hist = sqrt((m.' - Centroid_vector(k,2)) .^ 2 + (l - Centroid_vector(k,1)) .^ 2);%./Storage1(y,x,k);
    mask(:,:,k) = Dist_hist;
    
end

for i = 1:A1
    for j = 1:B1
        for k = 1:C
            if k <= 8
                if (mask(i,j,k) > 50)
                    mask(i,j,k) = 0;
                else 
                    mask(i,j,k) = 1;
                end
            else
                if (mask(i,j,k) > 35)
                    mask(i,j,k) = 0;
                else 
                    mask(i,j,k) = 1;
                end
            end
            ROI_slices(i,j,k) = Average(i,j,k) * mask(i,j,k);
            ROI_H1(i,j,k) = mask(i,j,k) * Storage_num6(i,j,k);
            %Storage_num5(i,j,k) = mask(i,j,k) * Storage_num4(i,j,k);
        end
    end
end 

end