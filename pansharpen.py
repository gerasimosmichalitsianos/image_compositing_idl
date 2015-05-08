#!/usr/bin/env python
import os
import sys
import gc 
import numpy as np
import warnings as warn 
from osgeo import osr, gdal
from scipy.misc import imresize 

class SatelliteImage(object):

    """ 
    This general class contains various methods that act on typical satellite 
    imagery. This class makes use coarser red, green, blue, NIR, and 
    panchromatic bands. Methods include Brovey and PCA pan-sharpening 
    methods, NDVI, and other attributes. 

    @author: 
    Gerasimos Andreas Michalitsianos 
    gerasimosmichalitsianos@gmail.com 
    May 2015 
    """ 
    
    def __init__(self, dspan, dsmulti, red, green, blue, nir, pan):

        # ---- instantiate regular bands arrays, GDAL objects for map info. 
        self.red = red
        self.green = green
        self.blue = blue
        self.nir = nir
        self.pan = pan
        self.dspan = dspan
        self.dsmulti = dsmulti 

    def ndvi(self):
      '''
      Compute NDVI = (NIR-red)/(NIR+red)
      Good for mapping vegetation. 
      '''
        with warn.catch_warnings():
            warn.filterwarnings('ignore',category=RuntimeWarning)
            fltred = np.array(self.red,dtype=np.float32)
            fltnir = np.array(self.nir,dtype=np.float32)
            return (fltnir - fltred)/(fltred + fltnir)

    def saveimg(self, img, outname):
      
        '''
        This method uses GDAL tools to save a 1-band Geotiff. The input argument
        'img' should be a 2D numpy array. Outname should be a string ending in ".tif." 
        '''
      
        driv = gdal.GetDriverByName('GTiff')
        if img.shape == self.pan.shape:
            dst = driv.Create(outname, self.dspan.RasterXSize, self.dspan.RasterYSize, 1, gdal.GDT_Float32)
            dst.SetGeoTransform(self.dspan.GetGeoTransform())
            dst.SetProjection(self.dspan.GetProjection())
            dst.GetRasterBand(1).WriteArray(img)
            dst=None
            del dst
        elif img.shape == self.red.shape:
            dst = driv.Create(outname, self.dsmulti.RasterXSize, self.dsmulti.RasterYSize, 1, gdal.GDT_Float32)
            dst.SetGeoTransform(self.dsmulti.GetGeoTransform())
            dst.SetProjection(self.dsmulti.GetProjection())
            dst.GetRasterBand(1).WriteArray(img)
            dst=None
            del dst

    def pcasharpen(self,band):
      
        """ 
        This method uses the method of principal component analysis (PCA) 
        to perform pan-sharpening of multispectral imagery. The argument 
        'band' should be 'red' if one wants a pan-sharpened red image 
        returned, 'green' if one wants a pan-sharpened green image returned,
        'nir' if one wants a pan-sharpened NIR image returned, and 'blue' if
        one wants a pan-sharpened blue image returned. 
        
        """

        nrows,ncols = self.pan.shape

        redresized   = np.array(imresize(self.red, self.pan.shape, 'nearest'), dtype = np.float32)
        nirresized   = np.array(imresize(self.nir, self.pan.shape, 'nearest'), dtype = np.float32)
        blueresized  = np.array(imresize(self.blue, self.pan.shape, 'nearest'), dtype = np.float32)
        greenresized = np.array(imresize(self.green, self.pan.shape, 'nearest'), dtype = np.float32)
        
        R1d = np.reshape(redresized,   (1,nrows*ncols))
        G1d = np.reshape(greenresized, (1,nrows*ncols))
        B1d = np.reshape(blueresized,  (1,nrows*ncols))
        N1d = np.reshape(nirresized,   (1,nrows*ncols))

        # ---- store 1D vectors in (nrows by 4 column array)
        allbands1d = np.zeros((4,N1d.shape[1]), dtype = np.float32)
        allbands1d[0,:] = R1d
        allbands1d[1,:] = G1d
        allbands1d[2,:] = B1d
        allbands1d[3,:] = N1d

        # ---- compute covariance and eigenvectors/eigenvalues 
        covariance = np.cov(allbands1d)
        eigenvals, eigenvecs = np.linalg.eig(covariance)
        eigDiagonal = np.array([[eigenvals[0],0.0,0.0,0.0],\
            [0.0,eigenvals[1],0.0,0.0],\
            [0.0,0.0,eigenvals[2],0.0],\
            [0.0,0.0,0.0,eigenvals[3]]])

        dat = np.dot(np.dot(eigDiagonal,np.transpose(eigenvecs)), allbands1d)
        dat = np.transpose(dat)

        princ_component1 = np.reshape(dat[:,0], (nrows,ncols))
        princ_component2 = np.reshape(dat[:,1], (nrows,ncols))
        princ_component3 = np.reshape(dat[:,2], (nrows,ncols))
        princ_component4 = np.reshape(dat[:,3], (nrows,ncols))

        # ---- altered panchromatic band, so its mean/std
        # ---- align with to mean/std of first principal component 
        PC1avg = np.mean(princ_component1)
        PC1std = np.std(princ_component1)
        PANavg = np.mean(self.pan)
        PANstd = np.std(self.pan)

        # ---- substitute modified panchromatic band as 1st principal component 
        pan_altered = self.pan*(PC1std/PANstd) + PC1avg - (PC1std/PANstd)*PANavg
        panaltered_1d = np.reshape(pan_altered,(1,nrows*ncols))
        dat[:,0] = panaltered_1d

        w = np.dot(np.dot(eigenvecs,np.linalg.inv(eigDiagonal)), np.transpose(dat))
        w = np.array(np.transpose(w))

        Rfin = np.reshape(w[:,0], (nrows,ncols))
        Gfin = np.reshape(w[:,1], (nrows,ncols))
        Bfin = np.reshape(w[:,2], (nrows,ncols))
        NIRfin = np.reshape(w[:,3], (nrows,ncols))

        # The values obtained from the anti-transformation are not in the interval [0,1]:
        # so we have to interpolate them to bring them back in this range:
        Rfin   = (Rfin-np.min(Rfin))/(np.max(Rfin)-np.min(Rfin))
        Gfin   = (Gfin-np.min(Gfin))/(np.max(Gfin)-np.min(Gfin))
        Bfin   = (Bfin-np.min(Bfin))/(np.max(Bfin)-np.min(Bfin))
        NIRfin = (NIRfin-np.min(NIRfin))/(np.max(NIRfin)-np.min(NIRfin))

        if band == 'red':
            return Rfin
        if band == 'green':
            return Gfin
        if band == 'blue':
            return Bfin
        if band == 'nir':
            return NIRfin
        else:
            return np.zeros(self.pan.shape) 
        
def main():

    # ---- load RGB,NIR,pan images as arrays from Geotiff files 
    dsblue  = gdal.Open('LC81850332013230LGN00_B2_clipped.TIF')
    dsgreen = gdal.Open('LC81850332013230LGN00_B3_clipped.TIF')
    dsred   = gdal.Open('LC81850332013230LGN00_B4_clipped.TIF')
    dsnir   = gdal.Open('LC81850332013230LGN00_B5_clipped.TIF')
    dspan   = gdal.Open('LC81850332013230LGN00_B8_clipped.TIF')

    blue  = dsblue.GetRasterBand(1).ReadAsArray()
    green = dsgreen.GetRasterBand(1).ReadAsArray()
    red   = dsred.GetRasterBand(1).ReadAsArray()
    nir   = dsnir.GetRasterBand(1).ReadAsArray()
    pan   = dspan.GetRasterBand(1).ReadAsArray()

    # ---- initiate a SatelliteImage object
    mysatimg = SatelliteImage(dspan, dsred, red,green,blue,nir,pan)

    # ---- start using SatelliteImage object to do stuff
    
    ndvi = mysatimg.ndvi()  # compute NDVI and save it as a Geotiff 
    mysatimg.saveimg(ndvi, 'ndvi.tif')

    
    redsharp = mysatimg.pcasharpen('red')
    greensharp = mysatimg.pcasharpen('green')
    bluesharp = mysatimg.pcasharpen('blue')
    nirsharp = mysatimg.pcasharpen('nir')

    mysatimg.saveimg(redsharp,'redsharp.tif')
    mysatimg.saveimg(bluesharp,'bluesharp.tif')
    mysatimg.saveimg(greensharp,'greensharp.tif')

if __name__ == '__main__':
    main() 
            
