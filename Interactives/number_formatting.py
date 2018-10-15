#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
number_formatting

This library containsfunctions for format numbers into various useful formats.

Functions included
------------------
exp2HTML(num, sig_fig=2):
    Converts number into scientific notation string (and HTML string).
SigFig(num, sig_fig=2):
    Rounds a number into a set number of significant figures.

Created on Thu May 24 12:00:00 2018

@author: Sam Holen and Juan Cabanela
"""
import numpy as np


def exp2LaTeX(num, sig_fig=2):
    """
    This function takes any integer, float, or string number,
    so long as it is entered in a format that can be converted into
    a float and converts it into two strings, one for non-LaTeX scientific
    notation (e.g. 6.02 * 10^2) and one of the equivalent LaTeX.

    By default, rounds numbers to 2 sig figs. Can enter the desired number of
    sig figs as an optional argument.

    Example: To display the LaTeX version of a mole in a ipywidget Label,
    one can do the following.
        r'\({}\)'.format(exp2LaTeX(6.02e23)[1])

    Parameters
    ----------
    num : float
           a number (must be convertable to float with call to float(num))
    sig_fig : int
               the number of signficant figures to keep.

    Returns
    -------
    [new_num, new_num_LaTeX] : list of two strings
            new_num is the string equivalent of the number in scientific
            notation. new_num_LaTeX is the LaTeX string version of new_num.
    """

    if type(num) == str and not(num.isdigit()):
        TypeError("Invalid input num, must be a float or integer")
    else:
        num = float(num)

    # Return zero if that is the value
    if (num == 0):
        return [0, '0', '0']

    # Now deal with non-zero results
    if ((abs(num) < 1e-3) or (abs(num) >= 100000.0)):
        # Write out in scientific notation
        power = int(np.floor(np.log10(abs(num))))
        base = str(SigFig(num/pow(10, power), sig_fig))

        # Generate output
        new_num = base + ' * 10^' + str(power)
        new_num_LaTeX = base + ' \\times 10^{' + str(power) + '}'
        new_num_HTML = base + ' &times; 10<sup>' + str(power) + '</sup>'
    else:
        # Don't bother with scientific notation
        new_num = str(SigFig(num, sig_fig))
        new_num_LaTeX = new_num
        new_num_HTML = new_num
        
    return [new_num, new_num_LaTeX, new_num_HTML]


def SigFig(num, sig_fig=2):
    """
    This function takes any integer, float, or string number,
    so long as it is entered in a format that can be converted into
    a float and rounds it to 'digits' significant figures.

    Parameters
    ----------
    num : float
           a number (must be convertable to float with call to float(num))
    sig_fig : int
               the number of signficant figures to keep. (default 2)

    Returns
    -------
    new_num : float
               returns the number rounded to the proper number of significant
               figures.
    """

    if type(num) == str and not(num.isdigit()):
        TypeError("Invalid input num, must be a float or integer")
    else:
        num = float(num)
        
    # Confirm sig_fig value acceptable
    if (sig_fig < 1):
        sig_fig = 2

    # Compute log of number and use that to separate out the exponent and
    # coefficients
    num = float(num)  # Convert input to float

    # Return zero if that is the value
    if (num == 0):
        return num

    exponent = np.floor(np.log10(abs(num)))
    order = pow(10, np.floor(exponent))
    coeff = num/order

    # Create the format string for the coeff
    format_str = "{0:."+str(int(sig_fig-1))+"f}"
    new_num = float(format_str.format(coeff))*order

    if ((exponent > 13) or (exponent < -13)):
        # Returns scientific notation
        return float(new_num)
    elif (sig_fig <= exponent):
        return int(new_num)
    elif (exponent < 0):
        return float(round(new_num, sig_fig+int(abs(exponent))))
    else:
        return float(new_num)
