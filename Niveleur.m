function [nivelee] = Niveleur(Extracted_ROI, DC_AVG)

[A, B, C] = size(DC_AVG)

maxi = max(DC_AVG,[],C);

for i = 1:C
    
    for j = 1:A
        for k = 1:B
            if Extracted_ROI(j,k,i) == 0 
                Extracted_ROI(j,k,i) = maxi(i) * 0.3;
            end
        end
    end
end
                
nivelee = Extracted_ROI;