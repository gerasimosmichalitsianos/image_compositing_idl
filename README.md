INSTRUCTIONS: IMAGE COMPOSITING 

The purpose of this code is to perform image compositing. That is, one single band Geotiff image file is used to fill in 
null data gaps in another 1-band Geotiff image file. Null data gaps in Landsat 7 imagery for example, have resulted from a 
spring 2003 scanline malfunction. So for all Landsat 7 imagery after May 31st 2003, have linear features present, going 
across the imagery, whose pixel values are null (no-data, usually a pixel value of -9999). Null data gaps also may have 
been marked manuallywith a cloud mask or cloud-shadow mask. These IDL files just show a quick example of how to composite 
a pair of two 1-band Geotiff image files. Each of these images may for example, be two Band 3 image files (red band for 
Landsat 7) from two separate scenes close to one another in date. 

In this repository are three IDL files. These files have ".pro" extensions. To use this code overall, both IDL (8.1+) and ENVI (5.0+) should both be installed together on your system. To use his code, start IDL, compile the source files, and at the Linux
command line, add the names of your two Geotiff filename input arguments. Please see below: 

usage: 
$ idl 
IDL> envi ; start envi 
IDL> .compile compositescenes
IDL> .compile createtile
IDL> .compile globalcomposite
IDL> exit
$ idl -e "compositescenes" -args fname1.tif fname2.tif 

@author: 
Gerasimos Michalitsianos
Science Systems and Applications, Inc. 
June 2015 
