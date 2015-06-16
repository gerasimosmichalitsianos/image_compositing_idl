pro createtile, Img1, Img2, tcol, trow, dim, maxcol, maxrow, tileImg1, tileImg2

  ; 
  ; This IDL script returns two small tiles (2D IDL arrays)
  ; at the same locations and dimensions in both the larger
  ; scenes (arrays) Img1 and Img2. This script makes it so that
  ; Img1 and Img2 are not indexed out of bounds. 
  ;
  
  tcol2 = (tcol + dim - 1) < (maxcol - 1)
  trow2 = (trow + dim - 1) < (maxrow - 1)

  tileImg1 = Img1[tcol:tcol2,trow:trow2]
  tileImg2 = Img2[tcol:tcol2,trow:trow2]

end
