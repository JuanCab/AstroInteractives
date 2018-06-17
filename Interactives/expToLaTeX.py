#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
expToLaTeX

This library contains a function to convert a number into a equivalent
scientific notation string (and LaTeX string) to display.

Functions included
------------------
exp2LaTeX(num, sig_fig=2):
    Converts number into equivalent scientific notation string (and LaTeX string).

Created on Thu May 24 12:00:00 2018

@author: Sam Holen
"""

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
        print("Invalid input")
    else:
        num = float(num)
    if 'e' in str(num):
        num = str(num)
        stop = num.find("e")
        if '+' in num:
            start = num.rfind("+")+1
        elif '-' in num:
            start = num.rfind('-')
        else:
            start = num.rfind('e')+1
        base = num[:stop]
        power = num[start:]
        if len(base) >= sig_fig:
            base = ('{'+':'+'.'+str(sig_fig-1)+'f'+'}').format(float(base))
        new_num = base + ' * 10^' + power
        new_num_LaTeX = base + ' \cdot 10^{' + power + '}'
    elif num >= 100000.0:
        num = str(num)
        stop = len(num) - 1
        base = num[0] + '.' + num[1:num.find('.')] + num[num.find('.')+1:]
        power = str(len(num[1:num.find('.')]))
        if len(base)-1 >= sig_fig:
            base = ('{'+':'+'.'+str(sig_fig-1)+'f'+'}').format(float(base))
        new_num = base + ' * 10^' + power
        new_num_LaTeX = base + ' \cdot 10^{' + power + '}'
    elif num < 1:  
        if num == 0:
            new_num = str(num)
        else:
            num2 = str(num)
            i = 2
            while num2[i] == '0':
                i += 1
            start = i
            if len(num2[start:]) <= sig_fig:
               if len(num2[start:]) == sig_fig:
                   new_num = num2[:start+sig_fig]
               else:
                   new_num = num2
               
            elif num2[start+sig_fig] >= '5':    
                new_num = num2[:start+1] + str(int(num2[start+(sig_fig-1)])+1)
            else:
                new_num = num2[:start+sig_fig]            
        new_num_LaTeX = new_num
    else:
        num2 = str(num)
        if len(num2[:num2.find('.')]) > sig_fig:
            length = len(num2[:num2.find('.')])
            if num2[sig_fig] < '5':     
                new_num = str(round(num))
            else:
                new_num = num2[:sig_fig-1] + str(int(num2[sig_fig-1])+1) + ('0'*(length-sig_fig))
        elif len(num2[:num2.find('.')]) == sig_fig:
            new_num = ('{:.0f}').format(num)

        else:
            if len(num2) - 1 >= sig_fig:    
                new_num = ('{'+':'+'.'+str(sig_fig-len(num2[:num2.find('.')]))+'f'+'}').format(num)
            else:
                new_num = num2
        new_num_LaTeX = new_num
    return [new_num, new_num_LaTeX]

# Test this library
if __name__ == '__main__':
    tests = [526553.5,.002,12.79,2.561e17,0]
    for i in tests:
        print(exp2LaTeX(i,12))




