
# coding: utf-8

# In[32]:


## Function to convert a number such as 6.02e23 to LaTeX display i.e. 6.02 * 10^2 ##

def exp2LaTeX(num, sig_fig=2):
    '''
    This function takes any integer, float, or string number,
    so long as it is entered in decimal notation or like 3.9e8, and 
    converts it to an equivalent string that displays in LaTeX
    scientific notation. By default, rounds numbers to 2 sig figs.
    Can enter the desired number of sig figs as an optional 
    argument.
    
    ## READ BEFORE YOU USE ##
    This function outputs a list of two values. The first is a universal
    text version of the number (e.g. 6.02 * 10^2) while the second is the
    code necessary to display the value in a widget in LaTeX notation.
    ## To use the LaTeX notation, one must do the following.
    r'\({}\)'.format(exp2LaTeX(6.02e23)[1])
    Adding this to a description or a label of a widget will result in the
    number being dislayed in LaTeX notation.
    '''
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

if __name__ == '__main__':
    tests = [526553.5,.002,12.79,2.561e17,0]
    for i in tests:
        print(exp2LaTeX(i,12))




