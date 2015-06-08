pro compositescenes, Img1, Img2, dim

  ;
  ; This IDL script is an example of how to use IDL and ENVI 
  ; functions to composite two 1-band Geotiff images. For example, 
  ; Landsat 7 data had a scanline malfunction in 2003. This 
  ; resulted in linear null data gap features in all Landsat 7 
  ; imagery after May 31st, 2003. The purpose of this code is to 
  ; use one image to fill in the null-data gaps in the other 
  ; image. So while pixels in one image may be null, pixel values
  ; in the other image may contain good pixel values. 
  ; 
  ; This code tiles up both images into 50x50 pixel tiles, and 
  ; computes local gain and bias values so as to perform a 
  ; histogram normalization to fill in the null data gaps, taking 
  ; normalized pixel values from one image and filling in the other 
  ; image. If no null pixels are found in a tile (in the image to 
  ; be filled), a global gain and bias are computed and used to 
  ; compute the values of the filler pixels. Both Geotiff images 
  ; should be from roughly the same time of year to avoid anomalies 
  ; in the outputted composited image. To run this code, for example, 
  ; replace the filename strings that you see below with 
  ; whatever 1-band Geotiff files you wish to composite. 
  ; 
  ; Both IDL and ENVI must be installed to run this script. 
  ; Be sure compositescenes.pro, createtile.pro, and 
  ; globalcomposite.pro are all present in the same directory. 
  ;
  ; usage for compositing two images: 
  ; $ idl 
  ; IDL> envi ; start ENVI 
  ; IDL> .compile compositescenes 
  ; IDL> .compile createtile 
  ; IDL> .compile globalcomposite 
  ; 
  ; @author: 
  ; Gerasimos Michalitsianos
  ; NASA/GSFC, Science Systems and Applications, Inc. 
  ; June 2015 


  ; navigate to appropriate directory where Geotiff files are located so IDL/ENVI may open those files
  cd,'C:\Users\gmichali\Desktop\compositing'

  ; open up first 1-band Geotiff image, extract its map and pixel information 
  envi_open_data_file, 'LE71850332012172ASN00_sr_band1.tif', r_fid = fid1, /tiff
  if (fid1 eq -1) then return 
  envi_file_query, fid1, ns=xsize1, nl=ysize1, nb=nb1
  dims1 = [-1, 0, xsize1-1, 0, ysize1-1]
  proj = envi_get_projection(fid=fid1, pixel_size=ps)
  
  ; open up second 1-band Geotiff image, extract its map and pixel information
  envi_open_data_file, 'LE71850332012188ASN00_sr_band1.tif', r_fid = fid2, /tiff
  if (fid2 eq -1) then return 
  envi_file_query,fid2, ns=xsize2, nl=ysize2, nb=nb2
  dims2=[-1, 0, xsize2-1, 0, ysize2-1]
  
  ; these two images should overlap approximately same area (same Landsat tile)
  ; but they will likewise be slightly different dimensions, so we need to perform a layer-stacking
  
  nb      = nb1+nb2
  fids    = lonarr(nb)
  pos     = lonarr(nb)
  dims    = lonarr(5,nb)
  fids[0] = fid1
  fids[1] = fid2
  pos[*]  = 0
  
  dims[0,0] = dims1
  dims[0,1] = dims2
  
  outproj = proj
  outname = 'composited_band.tif'
  out_ps  = ps
  
  envi_doit, 'envi_layer_stacking_doit', $ 
    fid=fids, pos = pos, dims = dims, $ 
    out_dt = 2, interp=2, out_ps = out_ps, $ 
    out_proj = outproj, r_fid = fidstacked, /in_memory 
    
  ; now get our arrays from the layer stacked virtual image in memory 
  envi_file_query, fidstacked, dims = stack_dims, $ 
    ns = stack_xsize, nl = stack_ysize, nb = stack_nb
  stackMapInfo = envi_get_map_info(fid = fidstacked)
  
  Img1 = envi_get_data(fid = fidstacked, dims = stack_dims, pos=0)
  Img2 = envi_get_data(fid = fidstacked, dims = stack_dims, pos=1)
  
  s = size(Img1,/dimensions)
  ncols = s[0] & nrows = s[1]
  
  ; define a tile size for computing local gain and bias, for compositing images 
  dim=50
  
  ; this is a list of pairs, where each pair is a [row,col] array
  nodes = list()
  for j = 0,ncols-1 do begin  ; cols
    for i = 0,nrows-1 do begin  ; rows
      
      ; add (col, row) to a list of nodes defining each tile
      if (j mod dim eq 0 and i mod dim eq 0) then begin 
        nodes -> add, [j,i] ; [col,row]
      endif    
      
    endfor 
  endfor 
  
  ; allocate memory to create arrays containing cols and rows for each node 
  cols = lonarr(n_elements(nodes))
  rows = lonarr(n_elements(nodes))
  
  ; extract the cols. and rows of each node from node list
  for i = 0,n_elements(nodes)-1 do cols[i] = (nodes[i])[0]
  for i = 0,n_elements(nodes)-1 do rows[i] = (nodes[i])[1]
  
  ; calculate Gain and Bias if no common non-null/non-cld. pixels found in a sub-tile
  GainBias = globalcomposite(Img1, Img2)
  
  if (GainBias[0] ne -1 and GainBias[1] ne -1) then begin 
    Gwhole = GainBias[0] 
    Bwhole = GainBias[1]
  endif else begin
    print, "Warning: gain and bias not successfully computed for whole images"
  endelse
    
  ; using the created list of nodes, print out each subarray
  ntiles = long(n_elements(nodes))
  count = long(0)
  for i = 0,long(ntiles)-1 do begin
  
    ; define sub-tile by node pixel in upper left (by its row & col) (t means tile)
    tcol = cols[i] & trow = rows[i]
    
    ; create the sub-tile from the node defined by tcol and trow
    createtile, $
      Img1, Img2, tcol, trow, dim, s[0], s[1], tileImg1, tileImg2
    
    ; for sub-tiles, find non-null and non-cloudy pixels cloudy to both imgs.
    GainBiasPixels = where(tileImg1 ne -9999 and tileImg2 ne -9999, npixels)
  
    ; for two sub-tiles: if # Gain/Bias pixels > 0, calculate G and B ...
    if (npixels gt 0) then begin
      
      ; calculate the gain and bias from the sub-tiles
      tImg1pixels = tileImg1[GainBiasPixels]
      tImg2pixels = tileImg2[GainBiasPixels]
      
      m1 = moment(tileImg1[GainBiasPixels], sdev = s1)
      m2 = moment(tileImg2[GainBiasPixels], sdev = s2)
      Gtile = s1/s2
      Btile = m1[0] - Gtile*m2[0]
  
      ; create filler scene with same dimensions as small sub-tiles
      FillerScene = Gtile*tileImg2 + Btile
           
      ; for sub-tile, find locations of pixels that need to be filled in
      Fill_Pixels = where(tileImg1 eq -9999 and tileImg2 ne -9999, n2pixels)
  
      ; if # of pixels to be filled in (in the subtile, that is) > 0 ... 
      if (n2pixels gt 0) then begin 
              
        ; fill in pixels in sub-tile
        tileImg1[Fill_Pixels] = FillerScene[Fill_Pixels]
              
        ; replace pixel area in original image with filled-in tile
        Img1[tcol,trow] = tileImg1
  
      endif
  
    ; if no common non-null & non-cldy. pixels found, 
    ; use gain and bias from whole images
    endif else begin
      
      ; use gain and bias calculated from whole image
      FillerScene = Gwhole*tileImg2 + Bwhole
            
      ; for sub-tile, find locations of pixels that need to be filled in
      Fill_Pixels = where(tileImg1 eq -9999 and tileImg2 ne -9999, n3pixels)

      if (n3pixels gt 0) then begin 
              
        ; fill in pixels in sub-tile
        tileImg1[Fill_Pixels] = FillerScene[Fill_Pixels]
                
        ; replace pixel area in original image with filled-in tile
        Img1[tcol,trow] = tileImg1
        
      endif 
  
     endelse
  endfor
  
  ; now write out composited image to a geotiff image file, get
  ; its map information from stacked image (in memory)
  get_lun, u
  openw, u, outname
  writeu, u, Img1
  
  envi_setup_head, fname = outname+'.tif', ns = stack_xsize, $ 
    nl = stack_ysize, nb = 1, data_type=2, interleave=0, $ 
    map_info = stackMapInfo, /write, $ 
    descrip = 'composited image' 
  
  close,u
  free_lun,u

end
