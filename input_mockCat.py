# coding: utf-8
import mockCat
catalog = '/Users/karenyng/Documents/Research/code/ResearchCode'+\
           '/nfwMCMC_TestData/testshapecat.txt'
M_200_pred = (1e0, 2e0) 
coord = ('ra','dec','z')
ellip_theo = ('e1_theo', 'e2_theo', 'de')
ellip_inststd = 0.
halos = (139, 30, 0.53, 138.916679, 29.9167, 0.53)
filename = 'mockCat_Will_test.txt'
coord = ('ra', 'dec', 'z')
mockCat.makeMockCat(catalog, M_200_pred, coord, ellip_theo, ellip_inststd, 
                    halos, filename)
#reload(mockCat)
#mockCat.makeMockCat(catalog, M_200_pred, coord, ellip_theo, ellip_inststd, halos, filename)
#M_200_pred[h]
#M_200_pred
#M_200_pred[0]
#M_200_pred[1]
#reload(mockCat)
#reload(mockCat)
#reload(mockCat)
#mockCat.makeMockCat(catalog, M_200_pred, coord, ellip_theo, ellip_inststd, halos, filename)
#from profiles import nfwparam
#nfwparam(M_200_pred[0],z_halo)
#reload(mockCat)
#mockCat.makeMockCat(catalog, M_200_pred, coord, ellip_theo, ellip_inststd, halos, filename)
#get_ipython().system(u'clear ')
#reload(mockCat)
#mockCat.makeMockCat(catalog, M_200_pred, coord, ellip_theo, ellip_inststd, halos, filename)
#ellip_meas = ('e1_theo', 'e2_theo', 'de')
#reload(mockCat)
#mockCat.makeMockCat(catalog, M_200_pred, coord, ellip_theo, ellip_inststd, halos, filename)
#reload(mockCat)
#mockCat.makeMockCat(catalog, M_200_pred, coord, ellip_theo, ellip_inststd, halos, filename)
#get_ipython().system(u'gvim mockCat_Will_test.txt')
#get_ipython().magic(u'pwd ')
#get_ipython().magic(u'save input_mockCat.py')
#get_ipython().magic(u'pinfo %save')
#get_ipython().magic(u'save input_mockCat.py n1-60')
