###### INSTRUCTIONS: IMAGE COMPOSITING 

       The purpose of this code is to perform image compositing. The IDL code files in 
       this repository represent just a simple example using two 1-band Geotiff files. 
       That is, one single band Geotiff image file is used to fill in null data gaps in 
       another 1-band Geotiff image file. Null data gaps in Landsat 7 imagery for example, 
       have resulted from a spring 2003 scanline malfunction. So for all Landsat 7 imagery 
       after May 31st 2003, have linear features present, going across the imagery, whose 
       pixel values are null (no-data, usually a pixel value of -9999). Null data gaps also 
       may have been marked manually with a cloud mask or cloud-shadow mask. These IDL files 
       just show a quick example of how to composite a pair of two 1-band Geotiff image files. 
       Each of these images may for example, be two Band 3 image files (red band for Landsat 7) 
       from two separate scenes close to one another in date. Below is one example of using 
       this code. 

![alt tag](https://i.imgur.com/3we6rtp.png)

       The left panel shows the base image (with null scanline data gaps to be filled in), 
       and the middle panel shows the image that was used to fill in the null data gaps in 
       the image in the left panel. Both the left and middle images came from two different
       Landsat 8 scenes. The right panel shows the filled-in, composited scene, resulting 
       from the histogram normalization and filling of the scene in the left panel. Both of
       these Landsat 8 scenes were over the Landsat tile identified by path 185, row 33
       (western Greece, islands of Kefalonia and Ithaki). In plain view, in the panels above, 
       one can see the cities of Sami and Karavomylos, both on Kefalonia Island, western Greece. 

###### usage: 

       $ idl
       IDL> envi
       IDL> .compile compositescenes
       IDL> .compile createtile
       IDL> .compile globalcomposite
       IDL> exit
       $ idl -e "compositescenes" -args fname1.tif fname2.tif 

###### @author:
       Gerasimos Michalitsianos
       SGT, Inc. 
       Febuary 2018 
