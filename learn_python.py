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
#
# One thing to say right away is that, unlike IDL, capitalization matters in
# Python: variable != Variable

def great_function(num, power, inverse=False): #keywords work just like IDL!
    """Another way to comment python code is like this. It's useful for
    multi-line comments. 

    It is good pratice to put a docstring like this at the start of every 
    function that explains inputs and outputs or w/e.

    This function takes in a number and raises it to some power. If inverse
    is True then the number will be raised to 1/power

    """

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


# Lists can contain any variable type you want, and even a mix of them.
# Here's a fun and important fact: EVERYTHING in Python can be assigned a
# variable. This means that lists can contain numbers, strings, characters,
# functions, classes, modules, lists, etc. This is very exciting.
#
# This is a good time to talk about loops in Python, which can be quite
# differnt from loops in other languages. In the most basic sense, a Python
# loops iterates over the elements in a list. Here's a very basic example that
# also shows off the versitility of lists:

def loop_test():
    """a simple function to show how Python loops work. Returns the contents
    of a list

    """
    # define the list, let's put multiple types in it
    l = [1,'yes',1e4,True,False]

    # now loop over that list
    for element in l:
        print element

    return

# Ok, that's cool, but what if you want to loop in the more "traditional
# sense? Easy! Use the range() function to generate a list and then loop over
# it:

def fib(n):
    """generates the first n fibbboonnaacci (sp?) numbers while teach about
    loops and lists! Any fibboonnnnnnnnaci program worth it's salt should be
    recursive, but w/e.

    Returns: a list containing the first n fibbonaciii (sp?) numbers and list
    with the output of range(n)

    """
    
    outlist = [0,1] # we have to cheat a little bit
    countlist = [] # we'll use this one to illustrate how range() works
    
    for i in range(n-2): # the n-2 is because we start with the first 2
                         # elements already
        
        countlist.append(i)
        outlist.append(outlist[-1] + outlist[-2])

    return outlist, countlist

# what's up with that return statement? It's True! You can return as many
# things as you want. There are few different ways to capture this output:

fiblist, rangelist = fib(5)

listtup = fib(5)

# listtup is what we call a tuple (don't ask how to pronounce it) Tuples are a
# lot like lists, but they defined with paranthesis rather than brackets, e.g.

l = [1,2,3] # a list
t = (1,2,3) # a tuple

l[1] == t[1] # True!

# There are some subtlties between lists and tuples that are not important
# right now. In general you should be using lists beacuse they have more
# features than tuples.
