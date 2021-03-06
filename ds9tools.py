'''
This program contains various routines to by used in conjunction with ds9.
'''
from __future__ import division
import numpy
import pylab
import pyfits
import sys
import matplotlib.cm as cm
import matplotlib.colors

def readregions(regfile):
    '''
    regfile = (string) the ds9 region file, assumes that it was written using
              'ds9' Format and 'image' Coordinate System

    Currently this function only works on circles, ellipse, and box regions
    '''
    # find all the circle regions
    circ = numpy.fromregex(regfile,r"circle\(([0-9]*\.?[0-9]+),([0-9]*\.?[0-9]+),([0-9]*\.?[0-9]+)",[('xc',numpy.float),('yc',numpy.float),('rc',numpy.float)])

    # find all the elliptical regions
    ellip = numpy.fromregex(regfile,r"ellipse\(([0-9]*\.?[0-9]+),([0-9]*\.?[0-9]+),([0-9]*\.?[0-9]+),([0-9]*\.?[0-9]+),([0-9]*\.?[0-9]+)",[('xc',numpy.float),('yc',numpy.float),('a',numpy.float),('b',numpy.float),('angle',numpy.float)])

    # find all the box regions
    box = numpy.fromregex(regfile,r"box\(([0-9]*\.?[0-9]+),([0-9]*\.?[0-9]+),([0-9]*\.?[0-9]+),([0-9]*\.?[0-9]+),([0-9]*\.?[0-9]+)",[('xc',numpy.float),('yc',numpy.float),('width',numpy.float),('height',numpy.float),('angle',numpy.float)])

    return circ, ellip, box


def circmask(reg,X,Y,x_ext,y_ext):
    '''
    returns a mask array where all pixels of X,Y coordinates are within the reg
    region.
    '''
    #only examine the area in the vacanity of the region
    mask = numpy.zeros(numpy.shape(X))
    xc = reg[0]
    yc = reg[1]
    rc = reg[2]
    #// by 1 to make valid index the subtract 1 to shift region to zero
    #index
    x0 = (xc-rc)//1-1
    x1 = (xc+rc)//1-1
    y0 = (yc-rc)//1-1
    y1 = (yc+rc)//1-1
    #correct for the possability that a region might extend outside
    #image
    if x0 < 0: x0=0
    if x1 > x_ext-1: x1=x_ext-1
    if y0 < 0: y0=0
    if y1 > y_ext-1: y1=y_ext-1
    #calculate radial separation
    r = numpy.sqrt((X[y0:y1,x0:x1]+1-xc)**2+(Y[y0:y1,x0:x1]+1-yc)**2)
    mask[y0:y1,x0:x1] = r <= rc
    return mask == 1

def ellipmask(reg,X,Y):
    if reg[4] != 0 and reg[4] != 360:
        print 'Error: program does not currently handle ellipses with angles'
    r = numpy.sqrt((X+1-reg[0])**2/(reg[2]**2*1.0)+(Y+1-reg[1])**2/(reg[3]**2*1.0))
    mask = r <= 1
    return mask == 1

def boxmask(box,X,Y):
    if box[4] != 0 and box[4] != 360:
        print 'Error: program does not currently handle boxes with angles'
    mask_xmax = (X+1)-box[0] <= box[2]/2.0
    mask_xmin = (X+1)-box[0] >= -box[2]/2.0
    mask_ymax = (Y+1)-box[1] <= box[3]/2.0
    mask_ymin = (Y+1)-box[1] >= -box[3]/2.0
    mask = mask_xmin*mask_xmax*mask_ymin*mask_ymax
    return mask == 1


def regcount(fitsfile,regfile,cubeindex=None):
    '''
    Calcluates the pixel counts within the ds9 regions.
    Input:
    fitsfile = (string) the fits file name
    regfile = (string) the ds9 region file, assumes that it was written using
              'ds9' Format and 'image' Coordinate System
    cubeindex = [int] If the fits file just contains one fits file then this
                input should be None.  If the fist file contains a Data Cube
                (i.e. various bundeled fits files then this input specifies
                which fits file of the Data Cube to operate on. Note that the
                first fits file of the Data Cube is cubeindex = 1.
    '''
    # import the fits file data
    hdulist = pyfits.open(fitsfile)
    img = hdulist[0].data
    N_axes = numpy.size(numpy.shape(img)) #number of axes in fits file
    # define the extents of the image
    if N_axes > 2 and cubeindex == None:
        print 'regcount: Error the fits file contains a Data Cube of multiple fits images. Please specify which fits image to operate on with cubeindex. Exiting.'
        sys.exit()
    elif N_axes > 2:
        x_ext = numpy.shape(img)[2]
        y_ext = numpy.shape(img)[1]
        cubeindex = cubeindex-1 #convert cubeindex to zero index notation
    elif N_axes == 2 and cubeindex != None:
        print 'regcount: Warning cubeindex != None yet fits file does not contain a Data Cube, continuing.'
        x_ext = numpy.shape(img)[1]
        y_ext = numpy.shape(img)[0]
    elif N_axes == 2:
        x_ext = numpy.shape(img)[1]
        y_ext = numpy.shape(img)[0]
    x_range = numpy.arange(x_ext)
    y_range = numpy.arange(y_ext)
    # Create the mesh grid that will define the coordinates of the pixels
    X,Y = pylab.meshgrid(x_range,y_range)

    # Read in all the region properties
    circ, ellip, box = readregions(regfile)

    # Sum pixel counts in all circle regions
    if numpy.shape(circ)[0] != 0:
        print '{0} == {1} in {2}'.format('Circle Region', 'Counts', 'Pixels')
        for i in numpy.arange(numpy.shape(circ)[0]):
            mask = circmask(circ[i],X,Y,x_ext,y_ext)
            if N_axes == 2:
                count = numpy.sum(img[mask])
            else:
                count = numpy.sum(img[cubeindex,mask])
            pixels = numpy.sum(mask)
            print '{0:3} == {1:10} in {2:10}'.format(i,count,pixels)
    # Sum pixel counts in all Elliptical regions
    if numpy.shape(ellip)[0] != 0:
        print '{0} == {1} in {2}'.format('Ellipse Region', 'Counts', 'Pixels')
        for i in numpy.arange(numpy.shape(ellip)[0]):
            mask = ellipmask(ellip[i],X,Y)
            if N_axes == 2:
                count = numpy.sum(img[mask])
            else:
                count = numpy.sum(img[cubeindex,mask])
            pixels = numpy.sum(mask)
            print '{0:3} == {1:10} in {2:10}'.format(i,count,pixels)
    # Sum pixel counts in all Box regions
    if numpy.shape(box)[0] != 0:
        print '{0} == {1} in {2}'.format('Box Region', 'Counts', 'Pixels')
        for i in numpy.arange(numpy.shape(box)[0]):
            mask = boxmask(box[i],X,Y)
            if N_axes == 2:
                count = numpy.sum(img[mask])
            else:
                count = numpy.sum(img[cubeindex,mask])
            pixels = numpy.sum(mask)
            print '{0:3} == {1:10} in {2:10}'.format(i,count,pixels)

def regcenter(fitsfile,regfile,cubeindex=None):
    '''
    Calcluates the centroid and error on the centroid within the ds9 regions.
    Input:
    fitsfile = (string) the fits file name
    regfile = (string) the ds9 region file, assumes that it was written using
              'ds9' Format and 'image' Coordinate System
    cubeindex = [int] If the fits file just contains one fits file then this
                input should be None.  If the fist file contains a Data Cube
                (i.e. various bundeled fits files then this input specifies
                which fits file of the Data Cube to operate on. Note that the
                first fits file of the Data Cube is cubeindex = 1.
    '''
    def centcalc(img,X,Y):
        '''
        img is the fits image data array
        X, Y are the coordinate arrays from numpy.meshgrid
        '''
        #total number of counts in the region
        N_counts = numpy.sum(img)
        #total number of pixels in the region
        N_pix = numpy.size(img)
        #calculate the x and y centroids
        c_x = numpy.sum(X*img)/N_counts
        c_y = numpy.sum(Y*img)/N_counts
        #calculate the mean pixel count
        mu = N_counts/N_pix
        #calculate the error on the x and y centroids
        sigma_c_x = numpy.sqrt(numpy.sum(numpy.abs(X)*(img-mu)**2)/N_pix)
        sigma_c_y = numpy.sqrt(numpy.sum(numpy.abs(Y)*(img-mu)**2)/N_pix)
        return c_x, sigma_c_x, c_y, sigma_c_y

    # import the fits file data
    hdulist = pyfits.open(fitsfile)
    img = hdulist[0].data
    N_axes = numpy.size(numpy.shape(img)) #number of axes in fits file
    # define the extents of the image
    if N_axes > 2 and cubeindex == None:
        print 'regcenter: Error the fits file contains a Data Cube of multiple fits images. Please specify which fits image to operate on with cubeindex. Exiting.'
        sys.exit()
    elif N_axes > 2:
        x_ext = numpy.shape(img)[2]
        y_ext = numpy.shape(img)[1]
        cubeindex = cubeindex-1 #convert cubeindex to zero index notation
    elif N_axes == 2 and cubeindex != None:
        print 'regcount: Warning cubeindex != None yet fits file does not contain a Data Cube, continuing.'
        x_ext = numpy.shape(img)[1]
        y_ext = numpy.shape(img)[0]
    elif N_axes == 2:
        x_ext = numpy.shape(img)[1]
        y_ext = numpy.shape(img)[0]
    x_range = numpy.arange(x_ext)
    y_range = numpy.arange(y_ext)
    # Create the mesh grid that will define the coordinates of the pixels
    X,Y = pylab.meshgrid(x_range,y_range)

    # Read in all the region properties
    circ, ellip, box = readregions(regfile)

    # Sum pixel counts in all circle regions
    if numpy.shape(circ)[0] != 0:
        print '{0} == {1} ; {2}'.format('Circle Region', 'x centroid', 'y centroid')
        for i in numpy.arange(numpy.shape(circ)[0]):
            mask = circmask(circ[i],X,Y,x_ext,y_ext)
            if N_axes == 2:
                img_tmp = img[mask]
            else:
                img_tmp = img[cubeindex,mask]
            c_x, sigma_c_x, c_y, sigma_c_y = centcalc(img_tmp,X[mask],Y[mask])
            print '{0:3} == {1:0.2f}+/-{2:0.2f} ; {3:0.2f}+/-{4:0.2f}'.format(i,c_x,sigma_c_x,c_y,sigma_c_y)
    # Sum pixel counts in all Elliptical regions
    if numpy.shape(ellip)[0] != 0:
        print '{0} == {1} ; {2}'.format('Ellipse Region', 'x centroid', 'y centroid')
        for i in numpy.arange(numpy.shape(ellip)[0]):
            mask = ellipmask(ellip[i],X,Y)
            if N_axes == 2:
                img_tmp = img[mask]
            else:
                img_tmp = img[cubeindex,mask]
            c_x, sigma_c_x, c_y, sigma_c_y = centcalc(img_tmp,X[mask],Y[mask])
            print '{0:3} == {1:0.2f}+/-{2:0.2f} ; {3:0.2f}+/-{4:0.2f}'.format(i,c_x,sigma_c_x,c_y,sigma_c_y)
    # Sum pixel counts in all Box regions
    if numpy.shape(box)[0] != 0:
        print '{0} == {1} ; {2}'.format('Box Region', 'x centroid', 'y centroid')
        for i in numpy.arange(numpy.shape(box)[0]):
            mask = boxmask(box[i],X,Y)
            if N_axes == 2:
                img_tmp = img[mask]
            else:
                img_tmp = img[cubeindex,mask]
            c_x, sigma_c_x, c_y, sigma_c_y = centcalc(img_tmp,X[mask],Y[mask])
            print '{0:3} == {1:0.2f}+/-{2:0.2f} ; {3:0.2f}+/-{4:0.2f}'.format(i,c_x,sigma_c_x,c_y,sigma_c_y)

def regminmax(fitsfile,regfile,cubeindex=None):
    '''
    Calcluates the pixel counts within the ds9 regions.
    Input:
    fitsfile = (string) the fits file name
    regfile = (string) the ds9 region file, assumes that it was written using
              'ds9' Format and 'image' Coordinate System
    cubeindex = [int] If the fits file just contains one fits file then this
                input should be None.  If the fist file contains a Data Cube
                (i.e. various bundeled fits files then this input specifies
                which fits file of the Data Cube to operate on. Note that the
                first fits file of the Data Cube is cubeindex = 1.
    '''
    def minmax(img,X,Y):
        '''
        img is the fits image data array
        X, Y are the coordinate arrays from numpy.meshgrid
        finds the pixel coordinates of the local min and max
        '''
        max_value = numpy.max(img)
        min_value = numpy.min(img)
        max_id = numpy.argmax(img)
        min_id = numpy.argmin(img)
        x_max = X[max_id]+1
        y_max = Y[max_id]+1
        x_min = X[min_id]+1
        y_min = Y[min_id]+1
        return min_value,x_min,y_min, max_value,x_max,y_max

    # import the fits file data
    hdulist = pyfits.open(fitsfile)
    img = hdulist[0].data
    N_axes = numpy.size(numpy.shape(img)) #number of axes in fits file
    # define the extents of the image
    if N_axes > 2 and cubeindex == None:
        print 'regcenter: Error the fits file contains a Data Cube of multiple fits images. Please specify which fits image to operate on with cubeindex. Exiting.'
        sys.exit()
    elif N_axes > 2:
        x_ext = numpy.shape(img)[2]
        y_ext = numpy.shape(img)[1]
        cubeindex = cubeindex-1 #convert cubeindex to zero index notation
    elif N_axes == 2 and cubeindex != None:
        print 'regcount: Warning cubeindex != None yet fits file does not contain a Data Cube, continuing.'
        x_ext = numpy.shape(img)[1]
        y_ext = numpy.shape(img)[0]
    elif N_axes == 2:
        x_ext = numpy.shape(img)[1]
        y_ext = numpy.shape(img)[0]
    x_range = numpy.arange(x_ext)
    y_range = numpy.arange(y_ext)
    # Create the mesh grid that will define the coordinates of the pixels
    X,Y = pylab.meshgrid(x_range,y_range)

    # Read in all the region properties
    circ, ellip, box = readregions(regfile)

    # Sum pixel counts in all circle regions
    if numpy.shape(circ)[0] != 0:
        print 'Circle Regions'
        print 'min/max = value @ (x, y)'
        for i in numpy.arange(numpy.shape(circ)[0]):
            mask = circmask(circ[i],X,Y,x_ext,y_ext)
            if N_axes == 2:
                img_tmp = img[mask]
            else:
                img_tmp = img[cubeindex,mask]
            min_value,x_min,y_min, max_value,x_max,y_max = minmax(img_tmp,X[mask],Y[mask])
            print 'Region {0}'.format(i)
            print 'min = {0} @ ({1}, {2})'.format(min_value,x_min,y_min)
            print 'max = {0} @ ({1}, {2})'.format(max_value,x_max,y_max)
    # Sum pixel counts in all Elliptical regions
    if numpy.shape(ellip)[0] != 0:
        print 'Ellipse Regions'
        print 'min/max = value @ (x, y)'
        for i in numpy.arange(numpy.shape(ellip)[0]):
            mask = ellipmask(ellip[i],X,Y)
            if N_axes == 2:
                img_tmp = img[mask]
            else:
                img_tmp = img[cubeindex,mask]
            min_value,x_min,y_min, max_value,x_max,y_max = minmax(img_tmp,X[mask],Y[mask])
            print 'Region {0}'.format(i)
            print 'min = {0} @ ({1}, {2})'.format(min_value,x_min,y_min)
            print 'max = {0} @ ({1}, {2})'.format(max_value,x_max,y_max)
    # Sum pixel counts in all Box regions
    if numpy.shape(box)[0] != 0:
        print 'Box Regions'
        print 'min/max = value @ (x, y)'
        for i in numpy.arange(numpy.shape(box)[0]):
            mask = boxmask(box[i],X,Y)
            if N_axes == 2:
                img_tmp = img[mask]
            else:
                img_tmp = img[cubeindex,mask]
            min_value,x_min,y_min, max_value,x_max,y_max = minmax(img_tmp,X[mask],Y[mask])
            print 'Region {0}'.format(i)
            print 'min = {0} @ ({1}, {2})'.format(min_value,x_min,y_min)
            print 'max = {0} @ ({1}, {2})'.format(max_value,x_max,y_max)


def makemaskfits(regfile,fitsout,naxis1,naxis2,binfactor=1,crval1=None,
                 crval2=None,crpix1=None,crpix2=None,cd1_1=None,cd1_2=None,cd2_1=None,cd2_2=None,comment=None):
    '''
    This function creates a fits file that represents masked and unmasked
    regions with pixel values of 0 and 1 respectively.

    Input:
    regfile = [string] file name of the region file that contains regions for
              masked areas. Regions should be defined in image coordinates.
    fitsout = [string] name of the output fits file made of zeros and ones
    naxis1 = number of x pixels for the image corresponding to the regfile
    naxis2 = number of y pixels for the image corresponding to the regfile
    binfactor = [integer] the factor by which to bin the pixels for the fitsout
                file. e.g.: if naxis1 = 100 and binfactor = 2 then the fitsout
                file will be an image with 50 pixels along the x-axis
    c_____ = elements of the WCS center and CD matrix for the image
             corresponding to the regfile. The binfactor parameter will be used
             to adjust the WCS for the fitsout file.
    comment = [string] Optional user specified comment to be added to the fits
              header
    '''
    # Read in all the region properties
    circ, ellip, box = readregions(regfile)

    # Create the coordinate array
    xbins = naxis1//binfactor
    ybins = naxis2//binfactor
    xarray = numpy.arange(xbins)
    yarray = numpy.arange(ybins)
    X, Y = numpy.meshgrid(xarray,yarray)

    ###
    ### Create the mask array
    ###
    # circle regions
    mask_circ = numpy.zeros(numpy.shape(X))
    N_circ = numpy.shape(circ)[0]
    if N_circ != 0:
        for i in numpy.arange(N_circ):
            print 'makemaskfits: processing circle region {0} of {1}'.format(i,N_circ)
            #only examine the area in the vacanity of the region. Subtract by 1
            #to make coordinates zero indexed
            mask_temp = numpy.zeros(numpy.shape(X))
            xc = (circ[i][0]-1)/binfactor
            yc = (circ[i][1]-1)/binfactor
            rc = circ[i][2]/binfactor
            #// by 1 to make valid index. Add a 100% buffer to the radius to
            #make sure that we don't round down to inside the mask radius
            x0 = (xc-2*rc)//1
            x1 = (xc+2*rc)//1
            y0 = (yc-2*rc)//1
            y1 = (yc+2*rc)//1
            #correct for the possability that a region might extend outside
            #image
            if x0 < 0: x0=0
            if x1 > xbins-1: x1=xbins-1
            if y0 < 0: y0=0
            if y1 > ybins-1: y1=ybins-1
            #calculate radial separation
            r = numpy.sqrt((X[y0:y1,x0:x1]-xc)**2+(Y[y0:y1,x0:x1]-yc)**2)
            mask_temp[y0:y1,x0:x1] = r <= rc
            mask_circ += mask_temp
    # Elliptical regions
    mask_ellip = numpy.zeros(numpy.shape(X))
    N_ellip = numpy.shape(ellip)[0]
    if N_ellip != 0:
        for i in numpy.arange(N_ellip):
            print 'makemaskfits: processing ellipse region {0} of {1}'.format(i,N_eillip)
            if ellip[i][4] != 0 and ellip[i][4] != 360:
                print 'Error: program does not currently handle ellipses with angles'
            r = numpy.sqrt((X+1/binfactor-ellip[i][0]/binfactor)**2/((ellip[i][2]/binfactor)**2*1.0)+(Y+1/binfactor-ellip[i][1]/binfactor)**2/((ellip[i][3]/binfactor)**2*1.0))
            mask_ellip += r <= 1
    # Box regions
    mask_box = numpy.zeros(numpy.shape(X))
    N_box = numpy.shape(box)[0]
    if N_box != 0:
        for i in numpy.arange(N_box):
            print 'makemaskfits: processing box region {0} of {1}'.format(i,N_box)
            if box[i][4] != 0 and box[i][4] != 360:
                print 'Error: program does not currently handle boxes with angles'

            #only examine the area in the vacanity of the region. Subtract by 1
            #to make coordinates zero indexed
            mask_temp = numpy.zeros(numpy.shape(X))
            xc = (box[i][0]-1)/binfactor
            yc = (box[i][1]-1)/binfactor
            xw = box[i][2]/binfactor
            yw = box[i][3]/binfactor

            #// by 1 to make valid index. Add a 100% buffer to the radius to
            #make sure that we don't round down to inside the mask radius
            x0 = (xc-xw)//1
            x1 = (xc+xw)//1
            y0 = (yc-yw)//1
            y1 = (yc+yw)//1
            #correct for the possability that a region might extend outside
            #image
            if x0 < 0: x0=0
            if x1 > xbins-1: x1=xbins-1
            if y0 < 0: y0=0
            if y1 > ybins-1: y1=ybins-1

            mask_xmax = X[y0:y1,x0:x1]-xc <= xw/2
            mask_xmin = X[y0:y1,x0:x1]-xc >= -xw/2
            mask_ymax = Y[y0:y1,x0:x1]-yc <= yw/2
            mask_ymin = Y[y0:y1,x0:x1]-yc >= -yw/2

            #mask_xmax = (X+1)-box[i][0]/binfactor <= box[i][2]/2.0/binfactor
            #mask_xmin = (X+1)-box[i][0]/binfactor >= -box[i][2]/2.0/binfactor
            #mask_ymax = (Y+1)-box[i][1]/binfactor <= box[i][3]/2.0/binfactor
            #mask_ymin = (Y+1)-box[i][1]/binfactor >= -box[i][3]/2.0/binfactor

            mask_temp[y0:y1,x0:x1] = mask_xmin*mask_xmax*mask_ymin*mask_ymax
            mask_box += mask_temp

    mask = mask_circ+mask_ellip+mask_box

    #bin the mask array by the binfator
    X_flat = numpy.reshape(X,(xbins*ybins,))
    Y_flat = numpy.reshape(Y,(xbins*ybins,))
    mask_flat = numpy.reshape(mask,(xbins*ybins,))

    mask_binned, xedges, yedges = numpy.histogram2d(Y_flat,X_flat,bins = (ybins,xbins),weights=mask_flat)
    #note that we reversed the order of the x,y to match fits indexing

    # Create the image array
    map_array = mask_binned == 0
    map_array = map_array*1.0 #convert bol to float for fits file

    # Create the fits file
    hdu = pyfits.PrimaryHDU(map_array)
    history = 'This fits file was generate by makemaskfits in ds9tools.py using the following inputs: makemaskfits(regfile={0},fitsout={1},naxis1={2},naxis2={3},binfactor={4},crval1={5},crval2={6},crpix1={7},crpix2={8},cd1_1={9},cd1_2={10},cd2_1={11},cd2_2={12},comment={13}'.format(regfile,fitsout,naxis1,naxis2,binfactor,crval1,crval2,crpix1,crpix2,cd1_1,cd1_2,cd2_1,cd2_2,comment)
    hdu.header.add_history(history)
    # Write optional comment line to the header
    if comment != None:
        hdu.header.add_comment(comment)
    # Write optional WCS to header
    if crval1!=None and crval2!=None and crpix1!=None and crpix2!=None and cd1_1!=None and cd1_2!=None and cd2_1!=None and cd2_2!=None:
        hdu.header.update('ctype1', 'RA---TAN')
        hdu.header.update('ctype2', 'DEC--TAN')
        hdu.header.update('crval1', crval1)
        hdu.header.update('crval2', crval2)
        hdu.header.update('crpix1', crpix1/binfactor)
        hdu.header.update('crpix2', crpix2/binfactor)
        hdu.header.update('cd1_1', cd1_1*binfactor)
        hdu.header.update('cd1_2', cd1_2*binfactor)
        hdu.header.update('cd2_1', cd2_1*binfactor)
        hdu.header.update('cd2_2', cd2_2*binfactor)
        #if verbose:
            #print 'The header of the output fits file is as follows:'
            #print hdu.header.ascardlist()
    hdu.writeto(fitsout,clobber=True)



def pointregions(prefix,ra,dec,style='diamond',color='green',size=11,objid=None,
                 wcs=True):
    '''
    Creates a ds9.reg file where each object input is represented by a point.
    prefix = [string] the prefix associated with the output file
    ra = [1D array of floats; units = degrees] RA of the objects
    dec = [1D array of floats; units=degrees] Dec of the objects
    style = ['arrow', 'box', 'boxcircle', 'circle', 'cross', 'diamond', 'x']
       the shape of the points
    color = ['black', 'white', 'red' , 'green', 'blue', 'cyan', 'magenta',
       'yellow'] color of the point
    size = [integer; units=pixels] the size of the point
    objid = [array of integers] the object id of each object, will be added to
       the text portion of each point
    wcs = [True or False] if True than ra and dec should be input as degrees,
       else they should be input with image pixel coordinates
    '''
    outputname = prefix+'_points.reg'
    F = open(outputname,'w')
    F.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1'+'\n')
    if wcs:
        F.write('fk5'+'\n')
    else:
        F.write('image'+'\n')
    for i in numpy.arange(numpy.size(ra)):
        if objid!=None:
            F.write('point({0:1.5f},{1:1.5f}) # point={2} {3} color={4} text='.format(ra[i],dec[i],style,size,color)+'{'+'{0:0.2f}'.format(objid[i])+'}\n')
        else:
            F.write('point({0:1.5f},{1:1.5f}) # point={2} {3} color={4}\n'.format(ra[i],dec[i],style,size,color))
    F.close()


def pointregions_scale(prefix,ra,dec,z,colormap=cm.jet,style='diamond',size=11,
                       label=True):
    '''
    Creates a ds9.reg file where each object input is represented by a point.
    The points are assigned a color based on z and the specified colormap.
    prefix = [string] the prefix associated with the output file
    ra = [1D array of floats; units = degrees] RA of the objects
    dec = [1D array of floats; units=degrees] Dec of the objects
    z = [1D array of floats] Some scaler to determine the color of the point
       based on the specified colormap.
    style = ['arrow', 'box', 'boxcircle', 'circle', 'cross', 'diamond', 'x']
       the shape of the points
    colormap = a matplotlib.cm
    size = [integer; units=pixels] the size of the point
    label = [True or False] will label the points according to their z value
    '''
    outputname = prefix+'_points.reg'
    #normalize the scalar
    z_norm = (z-numpy.min(z))/numpy.max(z-numpy.min(z))
    F = open(outputname,'w')
    F.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1'+'\n')
    F.write('fk5'+'\n')
    for i in numpy.arange(numpy.size(ra)):
        color = matplotlib.colors.rgb2hex(colormap(z_norm[i]))
        if label:
            F.write('point({0:1.5f},{1:1.5f}) # point={2} {3} color={4} text='.format(ra[i],dec[i],style,size,color)+'{'+'{0:0.2f}'.format(z[i])+'}\n')
        else:
            F.write('point({0:1.5f},{1:1.5f}) # point={2} {3} color={4}\n'.format(ra[i],dec[i],style,size,color))
    F.close()

def circleregions(prefix,ra,dec,radius,color='green',objid=None):
    '''
    Creates a ds9.reg file where each object input is represented by a circle.
    prefix = [string] the prefix associated with the output file
    ra = [1D array of floats; units = degrees] RA of the objects
    dec = [1D array of floats; units=degrees] Dec of the objects
    radius = [1D array of floats; units=arcsec] Radius of the circle
    color = ['black', 'white', 'red' , 'green', 'blue', 'cyan', 'magenta',
       'yellow'] color of the point
    size = [integer; units=pixels] the size of the point
    objid = [array of integers] the object id of each object, will be added to
       the text portion of each point
    '''
    outputname = prefix+'_circles.reg'
    F = open(outputname,'w')
    F.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1'+'\n')
    F.write('fk5'+'\n')
    for i in numpy.arange(numpy.size(ra)):
        if objid!=None:
            F.write('circle({0:1.5f},{1:1.5f},{2:1.2f}") # color={3} text='.format(ra[i],dec[i],radius[i],color)+'{'+'{0:0.2f}'.format(objid[i])+'}\n')
        else:
            F.write('circle({0:1.5f},{1:1.5f},{2:1.2f}") # color={3}\n'.format(ra[i],dec[i],radius[i],color))
    F.close()

def ellipseregions(prefix,ra,dec,e1,e2,radius,color='green',objid=None,cd=((-1,0),(0,1)),wcs=True):
    '''
    Creates a ds9.reg file where each object input is represented by a circle.
    prefix = [string] the prefix associated with the output file
    ra = [1D array of floats; units = degrees or pixels] RA of the objects, if
       wcs=True then it is assumed the units are degrees, else pixels.
    dec = [1D array of floats; units = degrees or pixels] Dec of the objects, if
       wcs=True then it is assumed the units are degrees, else pixels.
    e1 = [1D array of floats] e1 ellipticity component. Ellipticities should
       be defined with respect to the standard x,y coordinate system and it is
       assumed that ra and dec axes are orientated such that +x=-RA and +y=+Dec
       unless a cd matrix is input, otherwise the cd option should be specified
    e2 = [1D array of floats] e2 ellipticity component.
    radius = [1D array of floats; units=arcsec or pixels] Radius of the object,
       e.g. FWHM in units of arcsec if wcs=True, else pixels.
    color = ['black', 'white', 'red' , 'green', 'blue', 'cyan', 'magenta',
       'yellow'] color of the point
    size = [integer; units=pixels] the size of the point
    objid = [array of integers] the object id of each object, will be added to
       the text portion of each point
    cd = [((float,float),(float,float))] This is the CD matrix typically found
         in WCS header definitions which provides the scale and orientation of
         the RA,Dec coordinate system with the x,y pixel coordinate system such
         that [[x];[y]] = [[CD_11, CD_12];[CD_21, CD_22]] [[RA]; [Dec]]. For the
         purposes of this program the pixel scale is not important as we are
         only after the respective orientation of the two coordinates.
    wcs = [True or False] if True than ra and dec should be input as degrees,
         else they should be input with image pixel coordinates
    '''
    # Convert the ellipticity components to the reference RA,Dec frame
    cd = numpy.array(cd)
    # angle that will rotate N to +y and E to -x
    theta = numpy.arctan(cd[1,0]/cd[0,0])
    e1_wcs = e1*numpy.cos(-2*theta) + e2*numpy.sin(-2*theta)
    e2_wcs = -e1*numpy.sin(-2*theta) + e2*numpy.cos(-2*theta)

    # Calculate the position angle of the ellipse
    phi = numpy.arctan2(e2_wcs,e1_wcs)*180/numpy.pi
    phi /= 2.

    # Calculate the size of the ellipse components
    # calculate the ellipticity magnitude
    e_mag = numpy.sqrt(e1_wcs**2+e2_wcs**2)
    # calculate a and b for the ellipse
    a = radius*(1+e_mag)
    b = radius*(1-e_mag)

    # Save the ellipses to file
    outputname = prefix+'_ellipses.reg'
    F = open(outputname,'w')
    F.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1'+'\n')
    if wcs:
        F.write('fk5'+'\n')
        for i in numpy.arange(numpy.size(ra)):
            if objid!=None:
                F.write('ellipse({0:1.5f},{1:1.5f},{2:1.3f}",{3:1.3f}",{4:1.1f}) # color={5} text='.format(ra[i],dec[i],a[i],b[i],phi[i],color)+'{'+'{0:0.2f}'.format(objid[i])+'}\n')
            else:
                F.write('ellipse({0:1.5f},{1:1.5f},{2:1.3f}",{3:1.3f}",{4:1.1f}) # color={5}\n'.format(ra[i],dec[i],a[i],b[i],phi[i],color))
    else:
        F.write('image'+'\n')
        for i in numpy.arange(numpy.size(ra)):
            if objid!=None:
                F.write('ellipse({0:1.5f},{1:1.5f},{2:1.3f},{3:1.3f},{4:1.1f}) # color={5} text='.format(ra[i],dec[i],a[i],b[i],phi[i],color)+'{'+'{0:0.2f}'.format(objid[i])+'}\n')
            else:
                F.write('ellipse({0:1.5f},{1:1.5f},{2:1.3f},{3:1.3f},{4:1.1f}) # color={5}\n'.format(ra[i],dec[i],a[i],b[i],phi[i],color))
    F.close()
    print 'Finished making ellipse regions.'

def cropcontours(file_in, file_out, range_ra, range_dec):
    '''
    Will crop contours in a ds9 .con file to specified the ra and dec range.
    Input:
    file_in = [string], the path/name of the input ds9 contour file, should be
        saved using WCS degree options
    file_out = [string], the path/name of the cropped ds9 contour file
    range_ra = [(float, float); units=degrees] (min, max) ra range
    range_dec = [(float, float); units=degrees] (min, max) dec range
    '''
    import pandas
    # read the file into a pandas dataframe
    con = pandas.read_table(file_in,delimiter=' ', header=None,
                            skip_blank_lines=False)

    # Find all values that are outside the range_ra and range_dec
    mask_con_ra = numpy.logical_or(con[1]<range_ra[0],con[1]>range_ra[1])
    mask_con_dec = numpy.logical_or(con[2]<range_dec[0],con[2]>range_dec[1])
    mask_con_range = numpy.logical_or(mask_con_ra,mask_con_dec)
    # Convert the values outside the crop bounds into nulls.
    # This is a necessary process, rather than masking, so that all the
    # contours do not get connected together, and so that the the ends of
    # cropped contours do not get connected
    con[1][mask_con_range] = numpy.nan
    con[2][mask_con_range] = numpy.nan

    # due to formatting issues this pandas dataframe can't simply be written
    # using .to_csv
    array_ra = numpy.array(con[1])
    array_dec = numpy.array(con[2])
    F = open(file_out,'w')
    nanrow = True # a flag to see if
    for i in numpy.arange(numpy.shape(con)[0]):
        if numpy.isnan(array_ra[i]) and nanrow:
            # then this nan row was preceeded by a nan row or is the first row
            # that is also a nan row and shouldn't be written
            continue
        elif numpy.isnan(array_ra[i]):
            # then this is the first nan row in a series, write a blank line
            F.write('\n')
            nanrow = True
        else:
            # reset the nanrow switch
            nanrow = False
            # write the coordinates to the file
            F.write(' {0:1.8e} {1:1.8e} \n'.format(array_ra[i],array_dec[i]))
    F.close()







##debug
#fitsname = '/Users/dawson/Documents/Research/Filaments/StackWL/5to10a1_wNden_Nden.fits'
#regname = '/Users/dawson/Documents/Research/Filaments/StackWL/temp.reg'
#regcount(fitsname,regname,1)

#regname = '/Users/dawson/Documents/Research/Filaments/ExclusionRegions/full_field/150.1.0/F2R.150.1.0_reformat.reg'

#output = '/Users/dawson/Documents/Research/Filaments/ExclusionRegions/full_field/150.1.0/temp50_wcs.fits'
#naxis1 = 28132
#naxis2 = 28132
##makemaskfits(regname,output,naxis1,naxis2,binfactor=10)
#crval1=139.885
#crval2=30
#crpix1=14066
#crpix2=14066
#cd1_1=0.
#cd1_2=-7.13889E-05
#cd2_1=-7.13889E-05
#cd2_2=0.
#makemaskfits(regname,output,naxis1,naxis2,binfactor=50,crval1=crval1,crval2=crval2,crpix1=crpix1,crpix2=crpix2,cd1_1=cd1_1,cd1_2=cd1_2,cd2_1=cd2_1,cd2_2=cd2_2,comment=None)

##Debug regcenter
#fit = '/Users/dawson/Documents/Research/DLSCL09162953/WeakLensing/Centriod/Subaru/DLSCL09162953_053_100boot_15bin.fits'
#reg = '/Users/dawson/Documents/Research/DLSCL09162953/WeakLensing/Centriod/Subaru/temp.reg'
#ind = 1
##regcenter(fit,reg,ind)
#regminmax(fit,reg,ind)

## debug ellipsereg()
#import tools
#catalog_sub = '/Users/dawson/OneDrive/Research/ShapeComparison/PureCatalogs/SubaruPureCat_revB.txt'
#prefix = 'test'
#cat_sub = tools.readcatalog(catalog_sub)
#key_sub = tools.readheader(catalog_sub)
#ra_sub_id = 'alpha'
#dec_sub_id = 'delta'
#fwhm_sub_id = 'FWHM_IMAGE'
#e1_sub_id = 'e1'
#e2_sub_id = 'e2'
#ra = cat_sub[:,key_sub[ra_sub_id]]
#dec = cat_sub[:,key_sub[dec_sub_id]]
#e1 = cat_sub[:,key_sub[e1_sub_id]]
#e2 = cat_sub[:,key_sub[e2_sub_id]]
#size = cat_sub[:,key_sub[fwhm_sub_id]]
#pixscale_sub = 0.2
#size = size*pixscale_sub
#ellipseregions(prefix,ra,dec,e1,e2,size,color='magenta')

# # debug cropcontours()
# confile = '/Users/dawson/Git/Toothbrush-Dynamics-Paper/RedshiftAnalysis/TB_radio_GMRT610_debug.con'
# outfile = 'temp.con'
# ra_range = (90.604923243686727, 91.055076756313269)
# dec_range = (42.060555333333333, 42.393888666666662)
# cropcontours(file_in=confile, file_out=outfile, range_ra=ra_range, range_dec=dec_range)
