from __future__ import division
import sys
import numpy
'''
Meant to contain misc. tools
'''

def sort(array,sortcolumn,axis=0):
    '''
    Sorts the input 'array' according to the values in 'sortcolumn' number.
    By default it sorts the rows according to the values in a given column
    '''
    index = numpy.argsort(array,axis=axis)
    return array[index[:,sortcolumn],:]

def nearest(array,value):
    '''
    Returns the member of the 'array' that is nearest the specified 'value'.
    '''
    idx=(numpy.abs(array-value)).argmin()
    return array[idx]

def argnearest(array,value):
    '''
    Returns the index of the member of the 'array' that is nearest the specified
    'value'.
    '''
    idx=(numpy.abs(array-value)).argmin()
    return idx

def filtercat(cat,param_index,min_param,max_param,min_inc=True,max_inc=True,
              verbose=True):
    '''
    Filters the input catalog (cat) based on some minumum (min_param) and
    maximum (max_param) cluster pair parameter with column number (param_index)
    Returns a filterd catalog
    if min_inc == True: mask_min = cat[:,param_index] >= min_param
    elif min_inc == False: mask_min = cat[:,param_index] > min_param
    similar for max_inc
    '''
    Nint = len(cat)
    if min_param !=None and max_param != None:
	if min_param >= max_param:
	    print 'filtercat: Error, the minimum bound is greater than or equal to the maximum bound, exiting.'
	    sys.exit()
    if min_param != None:
        if min_inc == True:
            mask_min = cat[:,param_index] >= min_param
        elif min_inc == False:
            mask_min = cat[:,param_index] > min_param
    else:
        mask_min = True
    if max_param != None:
        if max_inc == True:
            mask_max = cat[:,param_index] <= max_param
        elif max_inc == False:
            mask_max = cat[:,param_index] < max_param
    else:
        mask_max = True
    if min_param != None or max_param != None:
        mask = mask_min*mask_max
        cat = cat[mask,:]
    if verbose:
        Nfin = len(cat)
        Ncut = Nint-Nfin
        print 'filtercat: {0} rows were removed from the catalog with {1} initial rows, leaving {2} rows'.format(Ncut,Nint,Nfin)
    return cat

def matchcat(cat0,index0,cat1,index1,nullval=-99):
    '''
    Use: matchcat(cat0,index0,cat1,index1,nullval=-99)
    
    Concatinates two catalogs by matching the rows of cat1 with those of cat0
    based on the values of cat0 column number index0 cat1 column number index1.
    cat0 is the primary array and all of its rows will be represented in the
    concatinated output array. Only the first row of cat1 that satisfies the 
    matching criteria will be concatinated with that row from cat0.
    
    cat0 = [1D or 2D array], primary base array
    index0 = [integer], 0 indexed column number of cat0 that contains the values
             to be matched with cat1
    cat1 = [1D or 2D array], secondary array
    index0 = [integer], 0 indexed column number of cat1 that contains the values
             to be matched with cat0
    nullval = [integer or float] the numerical value to assign to the columns of
              the output catalog, pertaining to cat1 columns, if there is no
              match for a given row
    Note that if a catalog is vector array then the param_index for that catalog
    should be entered as zero
    
    Output:
    Concatinated array with shape (Nrow0,Ncol0+Ncol1), the first columns will be
    those of cat0 followed by the columns of cat1.
    '''
    import numpy
    #check if either catalog is a vector array
    if numpy.shape(cat0)[0] == numpy.size(cat0):
        cat0 = numpy.reshape(cat0,(numpy.size(cat0),1))
    if numpy.shape(cat1)[0] == numpy.size(cat1):
        cat1 = numpy.reshape(cat1,(numpy.size(cat1),1))
    Nrow0 = numpy.shape(cat0)[0]
    Ncol0 = numpy.shape(cat0)[1]
    Nrow1 = numpy.shape(cat1)[0]
    Ncol1 = numpy.shape(cat1)[1]
    cat2 = numpy.ones((Nrow0,Ncol0+Ncol1))*nullval
    for i in numpy.arange(Nrow0):
        cat2[i,:Ncol0]=cat0[i,:]
        for j in numpy.arange(Nrow1):
            if cat0[i,index0] == cat1[j,index1]:
                cat2[i,Ncol0:]=cat1[j,:]
                break
    return cat2



def catfilter(min_param,max_param,cat,param_index,min_inc=True,max_inc=True,
              verbose=True):
    '''
    Filters the input catalog (cat) based on some minumum (min_param) and
    maximum (max_param) cluster pair parameter with column number (param_index)
    Returns a filterd catalog
    if min_inc == True: mask_min = cat[:,param_index] >= min_param
    elif min_inc == False: mask_min = cat[:,param_index] > min_param
    similar for max_inc
    '''
    Nint = len(cat)
    if min_param != None:
        if min_inc == True:
            mask_min = cat[:,param_index] >= min_param
        elif min_inc == False:
            mask_min = cat[:,param_index] > min_param
    else:
        mask_min = True
    if max_param != None:
        if max_inc == True:
            mask_max = cat[:,param_index] <= max_param
        elif max_inc == False:
            mask_max = cat[:,param_index] < max_param
    else:
        mask_max = True
    if min_param != None or max_param != None:
        mask = mask_min*mask_max
        cat = cat[mask,:]
    if verbose:
        Nfin = len(cat)
        Ncut = Nint-Nfin
        print 'catfilter: {0} rows were removed from the catalog with {1} initial rows, leaving {2} rows'.format(Ncut,Nint,Nfin)
    return cat

def angdist(ra1,dec1,ra2,dec2):
    '''
    Takes two spherical coordinates and determines the angluar distance between
    them.
    Input:
    ra1,dec1 [degrees] angular coordinates of the first point
    ra2,dec2 [degrees] angular coordinates of the second point
    Output:
    distance [degrees] angular distance between the two points
    '''
    from math import pi,sin,cos,acos
    d2r = pi / 180.0
    if ra1 < ra2:
        ra1,ra2 = ra2,ra1
    cosdist = cos((ra1-ra2)*d2r)*cos(dec1*d2r)*cos(dec2*d2r)+sin(dec1*d2r)*sin(dec2*d2r)
    distance = acos(cosdist) / d2r
    return distance


# UNITS: arcsec (output), degrees (input)
def angcomp(ra1,dec1,ra2,dec2,method=3):
    """
    Return the delta_RA and delta_dec (delta = 1-2) components of the angular
    distance between objects.  This is simply an alternate output of the
    angular_distance function above. Distance is returned as degrees, and method chooses a more or less accurate way of determining the distance (with 1 being the fastest/least accurate).
    from astorlib.py
    """
    DEGRAD = pi/180.
    import scipy
    if scipy.isscalar(ra1) and scipy.isscalar(ra2):
	from math import cos,sin,sqrt
	if ra1-ra2>180:
		ra1 -= 360.
	elif ra2-ra1>180:
		ra2 -= 360.
    else:
	from scipy import cos,sin,sqrt
	if scipy.isscalar(ra1):
		t = ra2.copy()
		ra2 = ra2*0. + ra1
		ra1 = t.copy()
		t = dec2.copy()
		dec2 = dec2*0. + dec1
		dec1 = t.copy()
		del t
	
	ra1 = scipy.where(ra1-ra2>180,ra1-360.,ra1)
	ra2 = scipy.where(ra2-ra1>180,ra2-360.,ra2)

    ra1 = ra1*DEGRAD
    dec1 = dec1*DEGRAD
    ra2 = ra2*DEGRAD
    dec2 = dec2*DEGRAD

    if method==1:
	deltadec = dec1-dec2
	deltara = (ra1-ra2)*cos(dec2)
    else:
	div = 1.
	if method==3:
		div = sin(dec2)*sin(dec1)+cos(dec2)*cos(dec1)*cos((ra1-ra2))
	deltara = cos(dec2)*sin(ra1-ra2)/div
	deltadec = -(sin(dec2)*cos(dec1)-cos(dec2)*sin(dec1)*cos(ra1-ra2)/div)
	#if sum(div == 0) != 0: #attempt at making array compatable but doesn't
	#work for single integers. Note that could just remove this section of
	#the code since it only here for QA
	if div == 0:
	    import sys
	    print 'tools: div = 0, exiting'
	    sys.exit()

    return deltara/DEGRAD, deltadec/DEGRAD

def angendpt(ra1,dec1,d_arcmin,pa):
    '''
    Given an endpoint coordinate (degrees), a distance (arcmin), and a position
    angle (degrees) ccw from the +dec axis it returns the coordinates of the 
    second endpoint.
    '''
    from math import pi,sin,cos,asin,acos
    d2r = pi / 180.0
    d = d_arcmin / 60.0
    phi = pa-90
    delta_dec = -asin(sin(phi*d2r)*sin(d*d2r))/d2r
    if d_arcmin >= 0: 
	if sin(pa*d2r) > 0: sign = 1
	else: sign = -1
    else:
	if sin(pa*d2r) > 0: sign = -1
	else: sign = 1
    delta_ra = sign*acos(cos(d*d2r)/cos(delta_dec*d2r))/d2r
    ra2 = ra1+delta_ra
    dec2 = dec1+delta_dec
    return (ra2,dec2)

def readheader(catalog):
    '''
    This function extracts the #ttype indexed header of a catalog.
    Input:
    catalog = [string], Name (perhaps including path) of the catalog
        that contains all of the data (e.g. x,y,e1,e2,...). Must include
        ttype header designations for the columns e.g.:
        #ttype0 = objid
        #ttype1 = x
    Output:
    dic = dictonary that contains the {ttype string,column #}.  
    '''
    import numpy
    import sys
    header = numpy.fromregex(catalog,r"ttype([0-9]+)(?:\s)?=(?:\s)?(\w+)",
                                 [('column',numpy.int64),('name','S20')])
    # Determine if the catalog is 0 or 1 indexed and if 1 indexed then change to 0
    if header['column'][0] == 1:
        header['column']-=1
    elif header['column'][0] != 0:
        print 'readheader: the catalog is not ttype indexed, please index using format ttype(column#)=(column name), exiting'
        sys.exit()
    for i in range(len(header)):
        if i == 0:
            dic = {header[i][1]:header[i][0]}
        else:
            dic[header[i][1]]=header[i][0]
    return dic

def readcatalog(catalog,verbose=True):
    '''
    This function extracts the #ttype indexed header of a catalog.
    Input:
    catalog = [string], Name (perhaps including path) of the catalog
        that contains all of the data (e.g. x,y,e1,e2,...). Must include
        ttype header designations for the columns e.g.:
        #ttype0 = objid
        #ttype1 = x
    Output:
    cat = [numpy array], numpy array formatted catalog.
    '''
    import numpy
    # Read in the catalog
    if verbose:
        print 'readcatalog: reading in '+catalog
    cat = numpy.loadtxt(catalog)
    if verbose:
        print 'readcatalog: read in '+catalog+' containing '+str(numpy.shape(cat)[0])+' rows and '+str(numpy.shape(cat)[1])+' columns of data'
    return cat

def distancematch(objid_1,ra_1,dec_1,objid_2,ra_2,dec_2,outputfile,tolerance = 2,stampsize=10):
    '''
    Given two catalogs, for each object in catalog 2 it will try to find the 
    nearest neighbor in catalog 1.  Then output a catalog of the matches.
    INPUT:
    objid_1 = [1D array] unique identifiers for the objects in catalog 1
    ra_1 = [1D array; units=degrees] Right Ascension of the objects in catalog 1
    dec_1 = [1D array; units=degrees] Declination of the objects in catalog 1
    objid_2 = [1D array] unique identifiers for the objects in catalog 2
    ra_2 = [1D array; units=degrees] Right Ascension of the objects in catalog 2
    dec_2 = [1D array; units=degrees] Declination of the objects in catalog 2
    outputfile = [string] name of the output file
    tolerance = [float; units=arcsec] any catalog 1 objects within this distance
        of the catalog 2 object will be considered a possible match
    stampsize = [float; units=arcsec] the square stamp size for which to
        consider catalog 1 objects as possible matches of the catalog 2 object.
	This is primarily meant to speed up the code.
    OUTPUT:
    text file of the matched objects
    '''
    #Create the ouput file and write header information
    fh = open(outputfile,'w')
    fh.write('#This catalog was created by the distancematch function of tools.py and matches objects in catalog 2 \n')
    fh.write('#with objects in catalog 1 based on their angular separation (matchdelta; units=arcsec) and user input.\n')
    fh.write('#ttype0 = objid_1\n')
    fh.write('#ttype1 = ra_1\n')
    fh.write('#ttype2 = dec_1\n')
    fh.write('#ttype3 = objid_2\n')
    fh.write('#ttype4 = ra_2\n')
    fh.write('#ttype5 = dec_2\n')
    fh.write('#ttype6 = matchdelta\n')
    
    N_2 = numpy.size(objid_2)
    for i in numpy.arange(N_2):
	# Do a quick trim of the catalog to reduce calculation time
	ramin = ra_2[i]-stampsize/(60.**2*2.*numpy.cos(dec_2[i]*numpy.pi/180.))
	ramax = ra_2[i]+stampsize/(60.**2*2.*numpy.cos(dec_2[i]*numpy.pi/180.))
	decmin = dec_2[i]-stampsize/(60.**2*2.)
	decmax = dec_2[i]+stampsize/(60.**2*2.)
	mask_ra = numpy.logical_and(ra_1>=ramin, ra_1<=ramax)
	mask_dec = numpy.logical_and(dec_1>=decmin, dec_1<=decmax)
	mask_stamp = numpy.logical_and(mask_ra,mask_dec)
	objid_flt = objid_1[mask_stamp]
	ra_flt = ra_1[mask_stamp]
	dec_flt = dec_1[mask_stamp]
	# determine the number of catalog 1 objects in the stamp
	N = numpy.size(objid_flt)
    
	# Calculated the angular separation between all catalog 1 objects and
	# the catalog 2 object
	j=0
	delta = numpy.zeros(N)
	for k in range(N):
	    delta[k] = numpy.abs(angdist(ra_flt[k],dec_flt[k],ra_2[i],dec_2[i])*60**2)
	    if delta[k] < tolerance:
		j+=1
	if j==1:
	    #there was a single match satisfying the tolerence
	    if N == 1:
		match_id = objid_flt[0]
		match_ra = ra_flt[0]
		match_dec = dec_flt[0]
		match_delta = delta[0]		
	    else:
		match_id = objid_flt[delta<tolerance][0]
		match_ra = ra_flt[delta<tolerance][0]
		match_dec = dec_flt[delta<tolerance][0]
		match_delta = delta[delta<tolerance][0]
	elif j == 0:
	    # then there were no objects in the postage stamp that were within
	    # the tollerance separation of the catalog 2 object
	    if numpy.size(delta) != 0:
		# then present the user with the nearest objects within the stamp
		# sort match_delta smallest to largest
		index = numpy.argsort(delta)
		delta=delta[index]
		objid_flt = objid_flt[index]
		ra_flt = ra_flt[index]
		dec_flt = dec_flt[index]
		print 'distancematch: No catalog matches were found for this catalog 2 object that satisfy the specified tollerance.'
		print 'Catalog2 objid: {0}'.format(objid_2[i])
		print 'The closest catalog 1 objects to the catalog 2 object are:'
		print 'Object\tobjid\tRA\t\tdec\tSeparation (arcsec)'
		for k in range(numpy.size(delta)):
		    print '{0}\t{1}\t{2:0.5f}\t{3:0.4f}\t{4:0.3f}'.format(k,objid_flt[k],ra_flt[k],dec_flt[k],delta[k])
		print '{0}\tSelect none.'.format(numpy.size(delta))
		selection = raw_input('Enter the number of the correct object match: ')
		if numpy.size(numpy.arange(k+1)==int(selection))==0:
		    selection = rawinput("Input invalid. Please enter a valid number.: ")
		if selection == str(numpy.size(delta)):
		    # Don't associate the trace with an object
		    match_id = match_ra = match_dec = match_delta = -99
		elif numpy.size(numpy.arange(k+1)==int(selection))!=0:
		    selection=int(selection)
		    match_id = objid_flt[selection]
		    match_ra = ra_flt[selection]
		    match_dec = dec_flt[selection]
		    match_delta = delta[selection]
	    else:
		# notify the user that there were no catalog 1 objects within
		# the stamp area
		print 'distancematch: No catalog 1 matches were found for this catalog 2 object in the stamp area, consider increasing the stamp size.'
		print 'Catalog2 objid: {0}'.format(objid_2[i])
		match_id = match_ra = match_dec = match_delta = -99
	elif j > 1:
	    objid_flt = objid_flt[delta<tolerance]
	    ra_flt = ra_flt[delta<tolerance]
	    dec_flt = dec_flt[delta<tolerance]
	    delta = delta[delta<tolerance]
	    print 'distancematch: More than one matches satisfy the separation tolerence.'
	    print 'Catalog2 objid: {0}'.format(objid_2[i])
	    print 'Match\tRA\t\tdec\tSeparation (arcsec)'
	    for k in range(j):
		print '{0}\t{1}\t{2:0.5f}\t{3:0.4f}\t{4:0.3f}'.format(k,objid_flt[k],ra_flt[k],dec_flt[k],delta[k])
	    print '{0}\tSelect none.'.format(j)
	    selection = raw_input('Enter the number of the correct match: ')
	    if numpy.size(numpy.arange(k+1)==int(selection))==0:
		selection = rawinput("Input invalid. Please enter a valid number.: ")
	    if selection == str(j):
		# Don't associate the catalog 2 object with a catalog 1 object
		match_id = match_ra = match_dec = match_delta = -99
	    elif numpy.size(numpy.arange(k+1)==int(selection))!=0:
		selection=int(selection)
		match_id = objid_flt[selection]
		match_ra = ra_flt[selection]
		match_dec = dec_flt[selection]
		match_delta = delta[selection]
	fh.write('{0:0.0f}\t{1:0.6f}\t{2:0.5f}\t{3:0.0f}\t{4:0.6f}\t{5:0.5f}\t{6:0.2f}\n'
	        .format(match_id,match_ra,match_dec,objid_2[i],ra_2[i],dec_2[i],match_delta))
    fh.close()
    print 'distancematch: process complete'

"""
Some celestial units tools.... taken from astrolib.py given by Dave L.
"""
from math import pi
# Check if the coordinate is in degrees...really just makes sure it's a float.
def is_degree(comp):
        if type(comp)==float:
                return True
        return False

# Convert string ra to degrees
def ra2deg(ra):
        comp = ra.split(" ")
        if comp[0]==ra:
                comp = comp[0].split(":")
        deg = float(comp[0])+float(comp[1])/60.+float(comp[2])/3600.
        return deg*15.

# Convert string declination to degrees
def dec2deg(dec):
        comp = dec.split(" ")
        if comp[0]==dec:
                comp = comp[0].split(":")
        if comp[0][0]=="-":
                comp[0] = comp[0][1:]
                sign = -1.
        else:
                sign = 1.
        return sign*(float(comp[0])+float(comp[1])/60.+float(comp[2])/3600.)

# Convert decimal ra to HMS format
def deg2ra(ra,sep=" "):
	from math import floor
        ra /= 15.
        h = floor(ra)
        res = (ra-h)*60.
        m = floor(res)
        s = (res-m)*60.
        if sep=="hms":
                sep1 = "h"
                sep2 = "m"
                sep3 = "s"
        else:
                sep1 = sep
                sep2 = sep
                sep3 = ""
        return "%02d%s%02d%s%06.3f%s" % (h,sep1,m,sep2,s,sep3)

# Convert decimal declination to DaMaS format
def deg2dec(dec,sep=" ",addsign=False):
	from math import floor
        if dec<0:
                sign = -1.
                dec = abs(dec)
        else:
                sign = 1.
        d = floor(dec)
        res = (dec-d)*60.
        m = floor(res)
        s = (res-m)*60.
        if sep=="dms":
                sep1 = "d"
                sep2 = "m"
                sep3 = "s"
        else:
                sep1 = sep
                sep2 = sep
                sep3 = ""
        if sign==-1:
                return "-%02d%s%02d%s%06.3f%s" % (d,sep1,m,sep2,s,sep3)
	if addsign:
		return "+%02d%s%02d%s%06.3f%s" % (d,sep1,m,sep2,s,sep3)
        return "%02d%s%02d%s%06.3f%s" % (d,sep1,m,sep2,s,sep3)

def addwcs(fitsfile,option,parentfits=None,ra_bounds=None,dec_bounds=None,
           unitcd=None,wcscards=None):
    '''
    fitsfile = [string] name of the fits file you want to add a wcs to
    option = [string], ('parent', 'bounds', or 'user'): 
             'parent' requires parentfits to be defined; it assumes that the two 
             fits files have the same extents and x,y orientation with respect to
             ra,dec and will create a matching wcs scale according to their size
             
             'bounds' requires ra_min, ra_max, dec_min, dec_max, and unitcd to be
             defined; creates a wcs based on input physical bounds and 
             orientation; currently it only allows orthogonal orientations of 
             x,y with respect to ra,dec
             
             'user' requires wcs to be defined; creates wcs header cards mathing
             the user defined values in the wcs dictionary
    parentfits = [string], name of the parent fits file from which to base the
                 current wcs on; only required if option='parent'
    ra_bounds = [(float,float)] {units: (degrees,degrees)}, (min ra, max ra); 
                only required if option='bounds'
    dec_bounds = [(float,float)] {units: (degrees,degrees)}, (min dec, max dec); 
                 only required if option='bounds'
    unitcd = [(int,int,int,int)], a unit cd matrix that defines the orientation of 
             the x,y coordinates with respect to the ra,dec coordinates. Currently only
             orthogonal relations are allowed. e.g. if positive x is right and
             positive y is up then (-1,0,0,1) sets the positive ra axis to the
             left and the positive dec axis up, or (0,1,1,0) sets the positive
             ra axis up and the positive dec axis to the right             
    wcscards = {'crval1':crval1,'crval2':crval2,'crpix1':crpix1,'crpix2':crpix2,
                'cd1_1':cd1_1,'cd1_2':cd1_2,'cd2_1':cd2_1,'cd2_2':cd2_2} a 
                dictionary defining the wcs cards of the header
    '''
    import pyfits
    import numpy
    import sys
    # read in the fits file to add the wcs to
    hdulist = pyfits.open(fitsfile, mode='update')
    hdr = hdulist[0].header
    naxis1 = hdr['naxis1']
    naxis2 = hdr['naxis2']
    if option == 'parent':
        if parentfits == None:
            print "addwcs: error, parentfits must be specified if option=='parent', exiting"
            sys.exit()
        phdulist = pyfits.open(parentfits)
        phdr = hdulist[0].header
        pnaxis1 = phdr['naxis1']
        pnaxis2 = phdr['naxis2']
        pcrval1 = crval1 = phdr['crval1']
        pcrval2 = crval2 = phdr['crval2']
        pcrpix1 = phdr['crpix1']
        pcrpix2 = phdr['crpix2']
        pcd1_1 = phdr['cd1_1']
        pcd1_2 = phdr['cd1_2']
        pcd2_1 = phdr['cd2_1']
        pcd2_2 = phdr['cd2_2']
        #scale the crpix value for
        pixscale1 = naxis1/pnaxis1
        pixscale2 = naxis2/pnaxis2
        crpix1 = pcrpix1*pixscale1
        crpix2 = pcrpix2*pixscale2
        cd1_1 = pcd1_1/pixscale1
        cd1_2 = pcd1_2/pixscale2
        cd2_1 = pcd2_1/pixscale1
        cd2_2 = pcd2_2/pixscale2
    elif option == 'bounds':
        if ra_bounds == None or dec_bounds == None or unitcd == None:
            print "addwcs: ra_bounds, dec_bounds, and unitcd must be specified if option=='bounds', exiting"
            sys.exit()
        ra_min = ra_bounds[0]
        ra_max = ra_bounds[1]
        dec_min = dec_bounds[0]
        dec_max = dec_bounds[1]
        avgdec = (dec_max + dec_min) / 2.0 * numpy.pi / 180.0 
        ra_length = (ra_max - ra_min) * numpy.cos(avgdec)
        dec_length = (dec_max - dec_min)
        crval1 = (ra_max - ra_min) / 2.0 + ra_min
        crval2 = (dec_max - dec_min) / 2.0 + dec_min
        crpix1 = (naxis1+1)/2.
        crpix2 = (naxis2+1)/2.
        cd1_1 = unitcd[0]*ra_length/naxis1
        cd1_2 = unitcd[1]*ra_length/naxis2
        cd2_1 = unitcd[2]*dec_length/naxis1
        cd2_2 = unitcd[3]*dec_length/naxis2
    elif option == 'user':
        if wcscards == None:
            print "addwcs: wcscards must be specified if option=='user', exiting"
            sys.exit()
        crval1 = wcscards['crval1']
        crval2 = wcscards['crval2']
        crpix1 = wcscards['crpix1']
        crpix2 = wcscards['crpix2']
        cd1_1 = wcscards['cd1_1']
        cd1_2 = wcscards['cd1_2']
        cd2_1 = wcscards['cd2_1']
        cd2_2 = wcscards['cd2_2']
    else:
	print "addwcs: invalid option entered, exiting"
	sys.exit()
    hdr.add_history("This wcs was created using Will Dawson's python program tools.addwcs")
    hdr.update('ctype1', 'RA---TAN')
    hdr.update('ctype2', 'DEC--TAN')
    hdr.update('crval1', crval1)
    hdr.update('crval2', crval2)
    hdr.update('crpix1', crpix1)
    hdr.update('crpix2', crpix2)
    hdr.update('cd1_1',cd1_1)
    hdr.update('cd1_2',cd1_2)
    hdr.update('cd2_1',cd2_1)
    hdr.update('cd2_2',cd2_2)
    hdulist.close()

def fitline(x,y):
    '''
    This is a simple least square line fit to points without error. 
    x and y are equal length 1D arrays
    Output:
    m = slope
    b = intercept
    '''
    A = numpy.vstack([x, numpy.ones(len(x))]).T
    m, b = np.linalg.lstsq(A, y)[0]
    return m, b

def coterminal(theta):
    '''
    This finds the minimum positive coterminal angle. For example if theta=370
    then the minimum positive coterminal angle = 10.
    Input:
    theta = [degrees] some angle in units of degrees
    Output:
    theta_min = [degrees] the minimum positive coterminal angle
    '''
    from math import pi
    if theta >= 0:
	k=0
	while k < 100:
	    theta_min = theta-k*360
	    k+=1
	    if theta_min < 360:
		break
	    if k==99:
		print 'coterminal: input angle is too large'
    else:
	k=-1
	while k > -100:
	    theta_min = theta-k*360
	    k-=1
	    if theta_min >= 0:
		break
	    if k == -99:
		print 'coterminal: input angle is too small'
    return theta_min


"""
Copyright (c) 2012, William A. Dawson
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the University of California, Davis nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""
