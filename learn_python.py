# Great! You've already learned something: comments in Python are made with a
# '#' character.
# I'll try to write this to be as general as possible, but there will be geared
# a bit towards people who are used to IDL.

# The first thing we'll do is import some useful modules into our Python 
# environment. A module is just a collection of functions and/or classes (we'll
# talk about what a class is in a bit). You simply import them with an "import"
# statement:
import pyfits
# you can also import them and give them a different name, as is commonly done 
# with Numpy:
import numpy as np # I could just have easily said: "import numpy as poop"
# You can also import only specific portions of a module:
import matplotlib.pyplot as plt # pyplot is a class within the matplotlib module


# Ok now let's define a function to do some cool stuff. We'll start with pure
# Python code and add in the cool science-related stuff in a bit

def greatFunction(num, power, inverse=False): #keywords work just like IDL!
    '''Another way to comment python code is like this. It's useful for
    multi-line comments. 

    It is good pratice to put a docstring like this at the start of every 
    function that explains inputs and outputs or w/e.

    This function takes in a number and raises it to some power. If inverse
    is True then the number will be raised to 1/power
    '''

    # let's see an if/else statement 
    # NOTE: Whitespace in Python is Very Important. It is the only way to
    # indicate the end of an if statement, loop, or whatever.
    if inverse:
        output = num**(1./power)
    
    else:
        output = num**power # this took me a long time to get used to. '^' is
                            # bitwise XOR in python so using won't raise any
                            # errors

    return output


# Ok. That was pretty fun, but what if we have a big groups of numbers we want
# to do that to. Python Lists are one way to do this, so let's look at them.
# The syntax is similar to IDL, but keep in mind that Python lists are VERY
# DIFFERNET than IDL lists in some very important ways.

# define a list

x = [1,2,3,4,5,6,7]

# slice it

x0 = x[0] # the first element of x
xl = x[-1] # the last element of x
x_sub = x[1:4] # another list that is [x[1],x[2],x[3]] NOTE the last slice
               # index is NOT included


