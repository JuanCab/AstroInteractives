import numpy as np

# I discovered Ballestros routines to convert B-V to temperature and back
# (Ballestreros et al. 2012) are used by PyAstronomy (props to Andrew
# Louwagie-Gordon for initial pointing).    In addition to that, I grabbed
# a python routine for getting RGB colors from Temp developed at
# https://stackoverflow.com/questions/21977786/star-b-v-color-index-to-apparent-rgb-color#22630970
# and added code defining that conversion that as well.

class BallesterosBV_T:
    """
      Black-body based conversion between effective temperature and B-V color.
      
        Ballesteros 2012 (EPL 97, 34008) present a conversion between
        effective temperature and B-V color index based on a black body
        spectrum and the filter functions.
    """

    def __init__(self):
        self._a = 0.92
        self._b = 1.7
        self._c = 0.62
        self._T0 = 4600.0

    def bv2T(self, bv):
        """
            Convert B-V color into temperature.
            
            Parameters
            ----------
            bv : float
                    B-V color index [mag]

            Returns
            -------
            T : float
                    Temperature [K]
        """
        T = self._T0*(1.0/(self._a*bv + self._b) + 1.0/(self._a*bv + self._c))
        return T

    def t2bv(self, T):
        """
            Convert temperature into B-V color.

            Parameters
            ----------
            T : float
                    Temperature in K.

            Returns
            -------
            bv : float
                    B-V color index [mag].
        """
        z = T / self._T0
        ap = z*self._a**2
        bp = self._a * self._c*z + self._b * self._a*z - 2.0 * self._a
        cp = self._b * self._c*z - self._c - self._b
        sqrtarg = bp**2 - 4.0 * ap * cp

        # By comparison with a BB it can be verified that
        # the physical solution is this one
        bv1 = (-bp + np.sqrt(sqrtarg))/(2.0*ap)
        return bv1


def bv2temp(bv):
    b = BallesterosBV_T()
    return b.bv2T(bv)
bv2temp.__doc__ = BallesterosBV_T.bv2T.__doc__


def temp2bv(T):
    b = BallesterosBV_T()
    return b.t2bv(T)
temp2bv.__doc__ = BallesterosBV_T.t2bv.__doc__



# In addition to the above code grabbed from PyAstronomy 0.12.0, I have added
# the following code for handling Temp to RGB conversion.  The original version 
# was developed at
# https://stackoverflow.com/questions/21977786/star-b-v-color-index-to-apparent-rgb-color#22630970
# although I modified it to handle nupy arrays as input and improved the fit
# (by updating the coefficients in redcoeff, greencoeff, and bluecoeff), since
# it really sucked at temperatures above 30000K. - Juan 05.24.2018

redcoeff = [5.342704572546931e-65, -1.7530489889997716e-59,
            2.611699966829653e-54, -2.336551831940032e-49,
            1.3988694428944896e-44, -5.910744925400012e-40,
            1.8111137280491087e-35, -4.069728884879208e-31,
            6.70106411740406e-27, -7.986900646511553e-23,
            6.718060628006256e-19, -3.818072040555441e-15,
            1.3634292254897947e-11, -2.681411890077735e-08,
            1.992979567493027e-05, 0.009860757283589523, 240.26480175678807]
greencoeff = [3.085231607019656e-65, -1.0026437222390378e-59,
              1.4740671747456953e-54, -1.2948517180175341e-49,
              7.556995796732612e-45, -3.0797149816066274e-40,
              8.951146405144254e-36, -1.855396384296827e-31,
              2.675257288412904e-27, -2.4876605560741676e-23,
              1.1179033888032595e-19, 3.7085233457054114e-16,
              -8.864398822838709e-12, 5.7326280719389936e-08,
              -0.00019235618674474583, 0.3734155034811879, -171.76671884893784]
bluecoeff = [1.188720685431585e-64, -4.003520387616542e-59,
             6.147034859201976e-54, -5.6968616597947265e-49,
             3.5563375780826955e-44, -1.580336787400662e-39,
             5.150962232882478e-35, -1.2505710564037003e-30,
             2.2739687882654258e-26, -3.0894110337858097e-22,
             3.106328242111663e-18, -2.2729529247194542e-14, 
             1.179761492097702e-10, -4.1754168041643567e-07, 
             0.0009361737030249371, -1.0933024353912841, 483.40830263785773]

redco = np.poly1d(redcoeff)
greenco = np.poly1d(greencoeff)
blueco = np.poly1d(bluecoeff)

def temp2rgb(temp):
    red = 0
    green = 0
    blue = 0

    """ Force values outside temperature bounds into proper range """
    temp = np.array(temp)
    temp[temp < 2000] = 2000
    temp[temp > 40000] = 40000
        
    """ since all lists have to have equal length, this is ok."""
    red = np.array(redco(temp))
    green = np.array(greenco(temp))
    blue = np.array(blueco(temp))

    red[np.round(red)>255] = 255
    red[np.round(red)<0] = 0
    green[np.round(green)>255] = 255
    green[np.round(green)<0] = 0
    blue[np.round(blue)>255] = 255
    blue[np.round(blue)<0] = 0
    
    red = red.astype('int')
    green = green.astype('int')
    blue = blue.astype('int')

    # Return a numpy array of R,G,B values    
    color = np.zeros([red.size,3], dtype='int')
    color[:,0] = red
    color[:,1] = green
    color[:,2] = blue
    return color


# Created a routine to take B-V straight to RGB
def bv2rgb(bv):
    return temp2rgb(bv2temp(bv))


# Convert RGB array (rows of R, G, B values) in np.array of hex codes
def rgb2hex(rgb):
    return np.array([ '#%02x%02x%02x' % (r,g,b)  for r, g, b in rgb ])
                     
"""
Code for making spectral plots in color is based on code found at
http://codingmess.blogspot.com/2009/05/conversion-of-wavelength-in-nanometers.html

I modified it to handle array inputs.
"""

def wav2rgb(wavelength):
    """
    Input wavelength is in nanometers (either scalar or integer)

    Return R, G, and B values in three column array color.
    """
    w = np.round(np.array(wavelength))
    color = np.zeros([w.size,3], dtype='float')
    SSS = np.zeros(w.size, dtype='float')
    
    # Assign RGB values over various ranges of wavelength
    range1 = (w >= 380) & (w < 440)
    color[range1, 0] = - (w[range1] - 440) / (440 - 350)
    color[range1, 1] = 0
    color[range1, 2] = 1

    range2 = (w>=440) & (w<490)
    color[range2, 0] = 0
    color[range2, 1] = (w[range2] - 440) / (490 - 440)
    color[range2, 2] = 1
    
    range3 = (w>=490) & (w<510)
    color[range3, 0] = 0
    color[range3, 1] = 1
    color[range3, 2] = - (w[range3] - 510) / (510 - 490)

    range4 = (w>=510) & (w<580)
    color[range4, 0] = (w[range4] - 510) / (580 - 510)
    color[range4, 1] = 1
    color[range4, 2] = 0
    
    range5 = (w>=580) & (w<645)
    color[range5, 0] = 1
    color[range5, 1] = -(w[range5] - 645) / (645 - 580)
    color[range5, 2] = 0

    range6 = (w>=645) & (w<=780)
    color[range6, 0] = 1
    color[range6, 1] = 0
    color[range6, 2] = 0

    # intensity corrections
    range1 = (w >= 380) & (w < 420)
    SSS[range1] = 0.3 + 0.7*(w[range1] - 350) / (420 - 350)    
    range2 = (w >=420) & (w <= 700)
    SSS[range2] = 1
    range3 = (w > 700) & (w <= 780)
    SSS[range3] = 0.3 + 0.7*(780 - w[range3]) / (780 - 700)

    return np.array(color*255*SSS[:,np.newaxis], dtype='int')

def wav2hex(wavelength):
    return rgb2hex(wav2rgb(wavelength))