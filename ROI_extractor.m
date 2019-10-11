function[FINAL_MASK,ROI_slices_ref,ROI_H1_ref] = ROI_extractor(ROI_slices, H1_ROI, Mask_ROI)

FINAL_MASK = zeros(size(ROI_slices));
ROI_slices_ref = zeros(size(ROI_slices));
ROI_H1_ref = zeros(size(ROI_slices));

[A, B, C] = size(ROI_slices);

for i = 1:C
    FINAL_MASK(:,:,i) = activecontour(H1_ROI(:,:,i), Mask_ROI(:,:,i), 100, 'edge',  'ContractionBias', 0.05);
end    

for i = 1:A
    for j = 1:B
        for k = 1:C
            ROI_slices_ref(i,j,k) = ROI_slices(i,j,k) * FINAL_MASK(i,j,k);
            ROI_H1_ref(i,j,k) = FINAL_MASK(i,j,k) * H1_ROI(i,j,k);
            if ROI_slices_ref(i,j,k) ~= 0
                ROI_slices_ref(i,j,k) = ROI_H1_ref(i,j,k) + ROI_slices_ref(i,j,k) ;
            end
        end
    end
end 

end