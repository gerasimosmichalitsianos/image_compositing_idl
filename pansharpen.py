
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

        '''
        This method allows for the initialization of a 
        SatelliteImage object. The dspan and dsmulti arguments 
        should be GDAL objects (as in dspan = gdal.Open('pan.tif') for example). 
        The red, green, blue, nir, and pan arguments should be 2D numpy arrays. 
        '''

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
        This instance method computes NDVI = (NIR - red)/(red + NIR)
        The NDVI image is returned as a 2D numpy array. 
        '''

        with warn.catch_warnings():
            warn.filterwarnings('ignore',category=RuntimeWarning)
            fltred = np.array(self.red,dtype=np.float32)
            fltnir = np.array(self.nir,dtype=np.float32)
            return (fltnir - fltred)/(fltred + fltnir)

    def saveimg(self, img, outname):

        ''' 
        This instance method writes out a Geotiff file. The inputs 
        are an img (a Numpy array) and an outname ('.tif'). 
        ''' 

        # ---- check inputs, crate GDAL driver object for writing out Geotiff 
        driv = gdal.GetDriverByName('GTiff')
        if not outname.endswith('.tif'): 
            outname=outname+'.tif'
        elif not ('str' not in str(type(outname))): 
            raise TypeError('saveimg(img,outname): second argument should be a string (.tif)')
        elif not ('numpy.ndarray' not in str(type(img))): 
            raise TypeError('saveimg(img,outname): first argument should be a Numpy array')
        
        # ---- compute the numpy of data layers 
        if len(img.shape) == 3: 
            nbands = 1 
        else: 
            nbands = img.shape[0]

        if img.shape == self.pan.shape:

            dst = driv.Create(outname, self.dspan.RasterXSize, self.dspan.RasterYSize, nbands, gdal.GDT_Float32)
            dst.SetGeoTransform(self.dspan.GetGeoTransform())
            dst.SetProjection(self.dspan.GetProjection())
            if nbands == 1: 
                dst.GetRasterBand(1).WriteArray(img)
            else: 
                for band in range(0,img.shape[0]+1): 
                    dst.GetRasterBand(band+1).WriteArray(img[band,:,:])

        elif img.shape == self.red.shape:

            dst = driv.Create(outname, self.dsmulti.RasterXSize, self.dsmulti.RasterYSize, nbands, gdal.GDT_Float32)
            dst.SetGeoTransform(self.dsmulti.GetGeoTransform())
            dst.SetProjection(self.dsmulti.GetProjection())

        if nbands == 1: 
            dst.GetRasterBand(1).WriteArray(img)
        else: 
            for band in range(0,img.shape[0]): dst.GetRasterBand(band+1).WriteArray(img[band,:,:])

        dst=None
        del dst

    def pcasharpen(self,band):
        
        ''' 
        This method performs principal component analysis (PCA) 
        pan sharpening. The only argument should specify which 
        band the user wants returned as a Numpy array, as in 
        'red' for the pan-sharpened red Numpy array, 'blue' 
        for the pan-sharpened blue array, 'nir' for the
        pan-sharpened NIR array. 
        ''' 
        
        if 'str' not in str(type(band)): 
            raise TypeError("pcasharpen('red/green/blue/nir') should be a string")
        else: 
            band = band.lower() 

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
            
