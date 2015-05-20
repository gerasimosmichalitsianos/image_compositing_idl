
#!/usr/bin/env python
import os
import sys
import gc
import math 
import numpy as np
import warnings as warn 
from osgeo import osr, gdal
from scipy.misc import imresize 

class SatelliteImage(object):

    """ 
    This general class contains various methods that act on typical satellite 
    imagery. This class makes use coarser red, green, blue, and NIR, bands, as well as 
    the finer panchromatic band. Methods include Brovey and PCA pan-sharpening 
    methods, NDVI, and other attributes. More to come. 
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
        are an image (a Numpy array, 2D or 3D) and an outname (should end in '.tif'). 
        ''' 

        # ---- check inputs, crate GDAL driver object for writing out Geotiff 
        driv = gdal.GetDriverByName('GTiff')
        if not outname.endswith('.tif'): 
            outname=outname+'.tif'
        elif type(outname) is not str: 
            raise TypeError('saveimg(img,outname): second argument should be a string (.tif)')
        elif type(img) is not np.ndarray: 
            raise TypeError('saveimg(img,outname): first argument should be a Numpy array')
        
        # ---- compute the numpy of data layers
        if len(img.shape)==3: 
            nbands = img.shape[0]
        else: 
            nbands = 1

        if img.shape[1:] == self.pan.shape: 
            dst = driv.Create(outname, self.dspan.RasterXSize, self.dspan.RasterYSize, nbands, gdal.GDT_Float32)
            dst.SetGeoTransform(self.dspan.GetGeoTransform())
            dst.SetProjection(self.dspan.GetProjection())

        elif img.shape[1:] == self.red.shape:
            dst = driv.Create(outname, self.dsmulti.RasterXSize, self.dsmulti.RasterYSize, nbands, gdal.GDT_Float32)
            dst.SetGeoTransform(self.dsmulti.GetGeoTransform())
            dst.SetProjection(self.dsmulti.GetProjection())

        if nbands == 1: 
            dst.GetRasterBand(1).WriteArray(img)
        else: 
            for band in range(img.shape[0]): dst.GetRasterBand(band+1).WriteArray(img[band,:,:])
        dst=None
        del dst

    def ihsSharpen(self):

        '''
        This method returns pan-sharpened NDVI resulting from 
        the IHS pan-sharpening of red and NIR bands. 
        
        Source paper: 
        Application of IHS-Sharpening Techniques to Ikonos Images 
        http://ieee.uniparthenope.it/chapter/_private/proc10/24.pdf
        
        '''

        # ---- first we have to resize RGB,NIR bands to same dimensions as panchromatic image using nearest neighbor interpolation
        nrows,ncols = self.pan.shape
        panfloat     = np.array(self.pan, dtype = np.float32)
        redresized   = np.array(imresize(self.red, self.pan.shape, 'nearest'), dtype = np.float32)
        blueresized  = np.array(imresize(self.blue, self.pan.shape, 'nearest'), dtype = np.float32)
        greenresized = np.array(imresize(self.green, self.pan.shape, 'nearest'), dtype = np.float32)

        # ---- [I,v1,v2] = [ 1/3, 1/3, 1/3 ; -sqrt(2)/6, -sqrt(2)/6, 2*sqrt(2)/6; 1/sqrt(2), 1/sqrt(2), 0] * [R; G; B]
        
        a = float(1)/3
        b = -1*math.sqrt(2.0)/float(6)
        c = 2.0*math.sqrt(2.0)/float(6)
        d = float(1)/math.sqrt(2)
        e = math.sqrt(2.0)

        I  = (a*redresized) + (a*greenresized)  + (c*blueresized)
        v1 = (b*redresized) + (b*greenresized)  + (c*blueresized)
        v2 = (d*redresized) + (-d*greenresized) + (0.0*blueresized)
        
        # ---- compute fused R,G,B
        redfused   = (1.0*panfloat) - (d*v1) + (d*v2)
        greenfused = (1.0*panfloat) - (d*v1) - (d*v2)
        bluefused  = (1.0*panfloat) + (e*v1) + (0.0*v2)

        # ---- return pan-sharpened RGB in 3-layer numpy array
        outsharpened = np.zeros((3,self.pan.shape[0],self.pan.shape[1]),dtype=np.float32)
        outsharpened[0,:,:] = bluefused
        outsharpened[1,:,:] = greenfused
        outsharpened[2,:,:] = redfused
        return outsharpened

    def fihsSharpen(self):

        '''
        This method performs a generalized fast IHS transformation to 
        pansharpen the RGB, NIR bands using the panchromatic band. 
        
        Source paper: 
        An IHS-Based Fusion for Color Distortion Reduction and Vegetation Enhancement in IKONOS Imagery 
        http://pgembeddedsystems.com/download/matlab/An%20IHS-Based%20Fusion%20for%20Color%20Distortion%20Reduction.pdf
        
        '''

        panfloat     = np.array(self.pan, dtype = np.float32)
        redresized   = np.array(imresize(self.red, self.pan.shape, 'nearest'), dtype = np.float32)
        nirresized   = np.array(imresize(self.nir, self.pan.shape, 'nearest'), dtype = np.float32)
        blueresized  = np.array(imresize(self.blue, self.pan.shape, 'nearest'), dtype = np.float32)
        greenresized = np.array(imresize(self.green, self.pan.shape, 'nearest'), dtype = np.float32)

        alpha = float(3)/4
        I = (redresized+greenresized+blueresized+nirresized)/float(4)
        redfused   = redresized   + (alpha*(panfloat - I))
        greenfused = greenresized + (alpha*(panfloat - I))
        bluefused  = blueresized  + (alpha*(panfloat - I))
        nirfused   = nirresized   + (alpha*(panfloat - I))

        # ---- return pan-sharpened RGB,NIR bands in 4-layer Numpy array
        outsharpened = np.zeros((4,self.pan.shape[0],self.pan.shape[1]),dtype=np.float32)
        
        outsharpened[0,:,:] = bluefused
        outsharpened[1,:,:] = greenfused
        outsharpened[2,:,:] = redfused
        outsharpened[3,:,:] = nirfused
        return outsharpened
        
    def pcaSharpen(self):
        
        ''' 
        This method returns a pan-sharpened NDVI using the pan-sharpening 
        technique of principal component analysis. 
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

        # ---- return 3D numpy array holding pan-sharpened RGB,NIR bands (4 layers)
        outsharpened = np.zeros((4,Rfin.shape[0],Rfin.shape[1]), dtype = np.float32)
        outsharpened[0,:,:] = Bfin
        outsharpened[1,:,:] = Gfin
        outsharpened[2,:,:] = Rfin
        outsharpened[3,:,:] = NIRfin
        return outsharpened

def main():

    panfname = '24DEC05QB020800005DEC24100239-P1BS-005682967010_01_P002_________AAE_0AAAAABAAAI0_TOA-Pan.tif'
    msfname = '24DEC05QB020800005DEC24100239-M1BS-005682967010_01_P002_________GA_E0AAAAAAGAAC0_TOA-Multispec.tif'
    
    dsms = gdal.Open(msfname)
    dspan = gdal.Open(panfname)

    nir = dsms.GetRasterBand(4).ReadAsArray() 
    red = dsms.GetRasterBand(3).ReadAsArray()
    green = dsms.GetRasterBand(2).ReadAsArray()
    blue = dsms.GetRasterBand(1).ReadAsArray()
    pan = dspan.GetRasterBand(1).ReadAsArray()

    satobj = SatelliteImage(dspan, dsms, red, green, blue, nir, pan)
    
    del nir
    del red
    del green
    del blue
    del pan
    gc.collect() 

    fihsSharpened = satobj.fihsSharpen()
    satobj.saveimg(fihsSharpened,'fihsSharpened.tif')
    
if __name__ == '__main__':
    main() 
