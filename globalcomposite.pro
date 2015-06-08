function globalcomposite, Img1, Img2

  ; from two input images (2D arrays), this function computes a global gain and bias 
  GainBiasPixels = where(Img1 ne -9999 and Img2 ne -9999, npixels)
  
  if (npixels gt 0) then begin
    
     ; calculate gain from pixels that are good in both images
     G = stddev(Img1[GainBiasPixels])/stddev(Img2[GainBiasPixels])
     
     ; calculuate the bias from pixels that are good in both images
     B = mean(Img1[GainBiasPixels]) - G*mean(Img2[GainBiasPixels])
     
     output = list()
     output -> add, G
     output -> add, B
     return, output
  
  endif else begin 
    
     output = list()
     output -> add, -1
     output -> add, -1
     return, output
     
  endelse 

end 
