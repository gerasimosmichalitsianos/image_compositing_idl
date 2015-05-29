pro pansharpenimage

  ;
  ; This IDL script is an example of how to use 
  ; ENVI functions to read Geotiff files in, and
  ; perform a Gram-Schmidt pan-sharpening. In this
  ; process, the multispectral Geotiff file 
  ; is pan-sharpened to the same dimensions as the 
  ; input panchromatic Geotiff image. To use this 
  ; script, both IDL and ENVI must be installed. 
  ; 
  ; usage for pan sharpening: 
  ; $ idl # start IDL
  ; IDL> .compile pansharpenimage
  ; IDL> pansharpenimage
  ; 
  ; @author: 
  ; Gerasimos Michalitsianos 
  ; NASA/GSFC, Science Systems and Applications, Inc. 
  ; May 2015 
  ;

  ; go to appropriate directory location for files, start ENVI
  cd,'J:\pansharpeningTesting'
  envi 

  ; get our Geotiff filenames 
  ; (1) panchromatic geotiff file (1 band)
  ; (2) multispectral geotiff file (4 bands, RGB, NIR)
  
  msfname = '29JAN07QB020800007JAN29102107-M1BS-005597176010_01_P001_________GA_E0AAAAAAGAAH0_TOA-Multispec.tif'
  panfname = '29JAN07QB020800007JAN29102107-P1BS-005597176010_01_P001_________AAE_0AAAAABAABD0_TOA-Pan.tif'

  ; open up the coarse multispectral file, get its image and map parameters
  envi_open_data_file, msfname, r_fid = fidlow, /tiff
  if (fidlow eq -1) then return 
  envi_file_query, fidlow, nb = nblow, dims = dimslow

  ; open up the high-resolution panchromatic image file, get its image and map parameters
  envi_open_data_file, panfname, r_fid = fidhigh, /tiff
  if (fidhigh eq -1) then return 
  envi_file_query, fidhigh, nb = nbhigh, dims = dimshigh

  ; create an outfile name
  outname = msfname.replace('.tif','_PanSharpenedGramSchmidt.tif')
  outname = outname.trim()
  pos=[0,1,2,3]

  ; now we are going to perform some pan-sharpening
  envi_doit, 'envi_gs_sharpen_doit', $

    fid = fidlow, $
    dims = dimslow, $
    pos = pos, $
    method=0, $
    out_name = outname, $
    interp=2, $
    hires_fid = fidhigh, $
    hires_pos = [0], $ 
    hires_dims = dimshigh

end 
