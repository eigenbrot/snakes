# Great! You've already learned something: comments in Python are made with the
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
import numpy as np # I could just have easily said: "import numpy as poop" You
# can also import only specific portions of a module although this is
# generally Very Bad:
import matplotlib.pyplot as plt # pyplot is a class within the matplotlib module



#################################
#         PYTHON BASICS         #
#################################
# Ok now let's define a function to do some cool stuff. We'll start with pure
# Python code and add in the cool science-related stuff in a bit
#
# One thing to say right away is that, unlike IDL, capitalization matters in
# Python: variable != Variable

def great_function(num, power, inverse=False): #keywords work just like IDL!
    """Another way to comment python code is like this. It's useful for
    multi-line comments. 

    It is good practice to put a docstring like this at the start of every 
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
                            # bit-wise XOR in python so using won't raise any
                            # errors

    return output

print "According to great_function, 4 to the 3rd power is {}\n".\
    format(great_function(4,3))
# Note the use of str.format() and \ to break up a line

# Ok. That was pretty fun, but what if we have a big groups of numbers we want
# to do that to. Python Lists are one way to do this, so let's look at them.
# The syntax is similar to IDL, but keep in mind that Python lists are VERY
# DIFFERENT than IDL lists in some very important ways.

# define a list

x = [1,2,3,4,5,6,7]

# slice it

x0 = x[0] # the first element of x
xl = x[-1] # the last element of x
x_sub = x[1:4] # another list that is [x[1],x[2],x[3]] NOTE the last slice
               # index is NOT included

# Before we talk about what makes Python lists awesome, let's talk about what
# they CAN'T do. The big one for us is that they can't be used like IDL
# lists. For example, here's what happens if I multiply x by a number:

print "x is {}\nand x*2 is {}\n".format(x,x*2)
# Hmm, maybe not what we expected. Don't worry, Numpy will save us from this
# in a moment.
#
# OK. Now on to why lists are awesome. Lists can contain any variable type you
# want, and even a mix of them.  Here's a fun and important fact: EVERYTHING
# in Python can be assigned a variable. This means that lists can contain
# numbers, strings, characters, functions, classes, modules, lists, etc. This
# is very exciting.
#
# This is a good time to talk about loops in Python, which can be quite
# different from loops in other languages. In the most basic sense, a Python
# loops iterates over the elements in a list. Here's a very basic example that
# also shows off the versatility of lists:

def loop_test():
    """a simple function to show how Python loops work. Returns the contents
    of a list

    """
    # define the list, let's put multiple types in it
    lst = [1,'yes',1e4,True,False]

    # now loop over that list
    for element in lst:
        print element

    return

print "loop_test() output:"
loop_test()
print # just a newline

# Ok, that's cool, but what if you want to loop in the more "traditional"
# sense? Easy! Use the range() function to generate a list and then loop over
# it:

def range_demo(n):
    """Returns the sum of the first n integers

    """
    
    output = 0
    countlist = [] # we'll use this one to illustrate how range() works
    
    for i in range(n):

        countlist.append(i)
        output += i

    return output, countlist

# what's up with that return statement? It's True! You can return as many
# things as you want. There are few different ways to capture this output:

result, rangelist = range_demo(5)

result_tup = range_demo(5)
print "the results of range_demo(5) are:\noutput: {}\ncountlist: {}\n".\
    format(*result_tup) # you can pass a list of arguments as a list or tup
                        # with a *

# result_tup is what we call a tuple (don't ask how to pronounce it). Tuples
# are a lot like lists, but they defined with parenthesis rather than
# brackets, e.g.

lst = [1,2,3] # a list
tup = (1,2,3) # a tuple

lst[1] == tup[1] # True!

# There are some subtleties between lists and tuples that are not important
# right now. In general you should be using lists because they have more
# features than tuples. A lot of numpy functions take tuples as arguments.




######################################################
#           A DIGRESSION: NAMESPACES                 #
######################################################
# Let's talk about namespaces and related ilk for a while.
#
# The Zen of Python states:
# "Namespaces are one honking great idea -- let's do more of those!"
# 
# The way I think of it is that every variable in Python (remember, EVERYTHING
# is a variable in Python) is scoped as narrowly as possible unless you
# explicitly decided to change that. For example, you might notice that I've
# defined a list called "lst" twice in this code; one within the function
# loop_test() and one outside of any function. This is totally OK (although
# maybe not the best practice) because they are in different
# namespaces. Variables assignment inside a function will never leave that
# function and variables defined within a module stay associated with that
# module.
#
# For example: To run the functions in this module (called learn_python.py)
# you will need to import this module into your Python session:
#
# >>>import learn_python as lp
#
# now everything in this file is contained in the lp namespace. Check it out:
# 
# --------
# >>>lst = ['not','the same list', 'as defined above']
# >>>lp.lst
# [1,2,3]
# >>>lst
# ['not','the same list', 'as defined above']
# lst == lp.lst
# False
# lp.loop_test()
# 1
# yes
# 10000.0
# True
# False
# >>>lp.loop_test.lst
# AttributeError: 'function' object has no attribute 'lst'
# ---------
#
# Notice how the lst defined in this module stays in this module. This is
# nice.  
# 
# IMPORTANT NOTE: It is possible to import individual items from a module into
# the CURRENT namespace and it is almost always a Bad Idea. For example, you
# could say:
# 
# >>>from learn_python import *
# 
# and now everything (functions, variables, etc.) in learn_python is in your
# "working" namespace. This means you don't need to prepend the functions with
# anything:
# 
# >>>range_demo(5)
#
# would work just fine. This is not good practice, however, because there
# could be anything in that module you're importing and who knows what
# functions and variables it might overload (re-define). In general you should
# only ever import modules.
#
# There are however, some acceptable exceptions that exist mostly for
# convenience and conformation with the community. For example, at the start
# of this module we imported the pyplot Class from the matplotlib
# Module. Importing a Class from a module is generally OK because Classes are
# their own separate namespaces (as we'll see later). Another type of exception
# is when you need only one function from a module. For example, I use the
# glob() function a lot to get lists of files and whatnot. The glob() function
# lives inside the glob module, but I don't want to type glob.glob() every time
# I use the function so I import it with:
from glob import glob

# If you want to be a strict follower of the "only import modules" camp (not a
# bad idea) but still want to have concise code you can just assign whatever
# function/class you want to use to a variable. For example:
#
# import glob
# glob = glob.glob
#
# will result in the exact same usage of glob for the programmer, but is a
# little more correct regarding namespaces and the like.
#
#
# One more (important) note. As we saw above, variables defined in a function
# stay in that function, but a function can still modify existing variables
# that are passed to it:

def inplace_demo(LL):
    """Designed to show how functions can modify variables outside of their
    namespace.

    Inputs: And python list
    Output: None
    
    """
    LL[-1] = 'poop'


# wow, now that is a concise program. It doesn't even need a return statement;
# the whitespace takes care of that. Let's try it out

lst1 = [4,5,6]
print "lst1 is now " + str(lst1) # another way to print things!
inplace_demo(lst1)
print "and after inplace_demo(lst1) lst1 is " + str(lst1) + '\n'

# will give us [4,5,'poop']. This is a potentially very useful but dangerous
# feature. Know about it. (I think IDL might do this too?)



###########################################
#             ENTER NUMPY                 #
###########################################
# OK. That was a nice little digression into some Python nuts and bolts, but we
# still haven't talked about how to apply our great_function() over a whole
# list!. The purely Python way might be fairly obvious at this point, just loop
# over the list:

def python_func(LL, power, inverse=False):
    """Same as great_function, but with the ability to operate on lists of
    numbers.
    
    Input: 
        LL - python list of numbers
        power - number
        inverse - If True then LL will be raised to 1/power
        
    Output:
        LL ** power (unless inverse=True)
        
    """
    
    # This function will not operate on LL inplace, so we need an output list
    
    output = []
    
    for elem in LL:
        if inverse:
            output.append(elem**(1./power))
        else:
            output.append(elem**power)
            
    return output

print "The output of our python_func(range(lst),3) is"
print python_func(lst,3)

# Alright, that works fine, but for long lists, or multi-dimensional lists it
# will take a long time. This is why Numpy was invented, it allows us to
# operate on every element in an array with single commands, and it happens
# much faster than it would if we used Python lists.
#
# The central component of Numpy is the ndarray type. These ndarrays (I'll
# call them arrays for short) act a lot like Python lists, except they let you
# use all of Numpy's power. Let's take a look:
# 
# The most basic way to make an array is to pass a Python list to np.array
array = np.array(lst)
# or make one up on the fly:
array2 = np.array([56,78,13])

# and now we don't even need that python_func()!

funcdarray = array**3 # this will cube every element of the array
print "...compare to this much faster numpy version:"
print funcdarray, '\n'

# THIS IS IMPORTANT! In general, if you find yourself looping over every
# element in an array (like python_func() ) then you are probably doing
# something wrong!

# There many many Numpy functions that use ndarray operations to make your
# lives easier. I'll mention the few that I use the most below:

r = np.arange(12) # just like Python's range(), but produces a ndarray
print "the size of r is {}".format(r.size) # ndarray.size will tell you how
                                           # many elements are in that array
r2d = r.reshape((4,3)) # creates a 3x3 array from r
print "and the shape of r2d is {}\n".format(r2d.shape) # this returns a tuple
                                                     # with the size of each
                                                     # dimension
idx = np.where(r > 3) # much like IDL's where
idx2 = np.where((r > 2) & (r < 7)) # note the slightly strange syntax. All of
                                   # those parenthesis are necessary. Also
                                   # note the use of & rather than the Python
                                   # 'and'
rmax = np.max(r) # np.min() also exists
rmax2 = r.max() # ndarry has min/max built in
rmean = np.mean(r)
rmedian = np.median(r)
rstd = np.std(r)
rsum = np.sum(r)
r2sum1 = np.sum(r2d,axis=0) # most of these commands have the axis keyword,
                            # which will collapse the array down along that
                            # axis

print "r2d is:\n", r2d # another way to print!!
print "the total sum is {}".format(np.sum(r2d))
print "and the sum just along axis 0 is: {}\n".format(r2sum1)

# Note that r.size = r2d.size, but r.shape != r2d.shape
#
# Remember N_ELEMENTS() in IDL? If you've ever written:
#
# FOR i=0, N_ELEMENTS(data) - 1 DO BEGIN
#
# FUCK THAT!
#
# you can use:
print "with range:"
for i in range(r.size):
    print r[i]

# but that also sucks! Remember, ndarrays are a lot like Python lists, and in
# Python for loops just iterate over any list. So try this:
print "with straight iteration:"
for i in r:
    print i

# nice!

# There a few important differences between ndarrays and Python lists.
# Probably the most obvious one from a coding standpoint is that very few
# Numpy functions operate on ndarrys in place. For example:
#
# >>> pythonlist.append(1)
#
# will actually change pythonlist. There are no functions in Numpy that do
# this. Instead, each function returns something, so if you want to change an
# array in place you need to redefine it. For example:
#
# >>> myarray = np.append(myarray, 1)
#
# This is not a concise as pure Python, but it does make it harder for you to
# unknowingly fuck up your data.



###################################
#      READING/WRITING DATA       #
###################################
# OK! We've learned the basics of ndarrays and the basics of how we can use
# them to manipulate our data. But how do we get our data into these arrays?
# Let's look at the two that are (probably?) the most common:
#
# From txt files (csv)
# --------------------
#
############### Reading
#
# np.loadtxt is your friend. Check >>>help(np.loadtxt) for more info but
# here's a basic primer:
#
# Most basic usage:
data = np.loadtxt('my_awesome_data.dat')

# Note that data is a 2D array, so that the first column in your data is
# accessed by data[0,:].
#
# By default, lines that start with '#' will be ignored and columns are
# separated by any whitespace. But maybe your data file's format uses
# something besides '#' as a comment and something besides whitespace as a
# delimiter. No worries:
data = np.loadtxt('my_unusual_data.poop',comments='>',delimiter=',')

# The final few features that I use often are selecting which columns you want
# to import and importing each column as a separate variable. The latter is
# particularly useful.
col1, col4, col5 = np.loadtxt('my_data.txt',usecols=(1,4,5),unpack=True)

# The unpack keyword is what allows us to assign multiple variables. Notice
# that usecols is a tuple! Tuples can be very useful when you want to allow
# the user to pass an arbitrary number of arguments to a function. In this
# case, the author of loadtxt had no idea how many columns a user might want
# to import, so she made that keyword a tuple, and tuples can be any length.
#
########### Writing
#
# To write to a csv file we essentially print each row of our columns to a
# file one at a time. There are few different ways to do this, but both first
# involve opening a file object:

f = open('writeout.txt','w') # The 'w' is for write
# First let's write some comments:

f.write('# This is a comment line\n'
        +'#\n'
        +'# {:6}= '.format('c1')+'column 1\n'
        +'# {:6}= '.format('c4')+'column 4\n'
        +'# {:6}= '.format('c5')+'column 5\n'
        +'#\n'
        +'#{:>10}{:>10}{:>10}\n'.format('c1','c4','c5')
        +'#{:>10}{:>10}{:>10}\n'.format(*(np.arange(3)+1))
        )

# Now we can actually write our data in two different ways. With Numpy's
# ndarray.tofile() method, or with python's str.format() method.
#
# Here's Numpy's way:

for i in range(col1.size):
    np.array([col1[i],
              col4[i],
              col5[i]]).tofile(f,sep='  ',format='%4.3f')

# and python's

for i in range(col1.size):
    f.write('{:4.3f}{:4.3f}{:4.3f}\n'.format(col1[i],
                                             col4[i],
                                             col5[i]))

# both are fine.


#
# From a FITS file
# ----------------
#
########### Reading
#
# The more an more I learn about FITS files and pyfits, the more and more I
# store my data in FITS files. Maybe you will too!
#
# First, let's just open a FITS file:
hdus = pyfits.open('my_great_fits_file.fits')

# hmm, what exactly have we got here?
print "pyfits generates data of type"+type(hdus)

# it's a pyfits.HDUlist!, which is really just a python list of individual
# Header Data Units (HDU). An HDU is the same as a FITS "extension" if you're
# familiar with IRAF. They always have two things, a header and some data. We
# can access either easily enough:
header = hdus[0].header
fitsdata = hdus[0].data # notice how we index the HDUList just like a python
                        # list

print "the data in that hdu is of type"+type(fitsdata)

# great, and now we've some normal Numpy data that we can use
#
# The header can also be useful. It acts a lot like a python dictionary. We
# haven't covered dictionaries yet, but in short they are like python lists
# that can be indexed with anything (not just integers that range from 0 to
# one less than the length of the list). The header is a dictionary of all the
# FITS header keywords that is indexed by those keywords. For example, let's
# say you have a keyword in your FITS header called 'EXPTIME' that tells you
# the integration time for that image. You can access this information via:
exptime = header['EXPTIME']

# Neat! Ok, those are the most basic things you might want to do with a FITS
# file in python. STSCI maintains Pyfits and they have a great document called
# "The_PyFITS_Handbook.pdf" if you want more info.
#
# A final note: Sometimes you just want that FITS data RIGHT NOW! and don't
# want to deal with HDULists and that other shit. Easy, just put all that
# stuff together:
fitsdata2 = pyfits.open('my_great_fits_file.fits')[0].data

# Final note v2: pyfits.open() ALWAYS returns an HDUList. Even when there is
# only one extension. So you will ALWAYS have to index the output of
# pyfits.open(). i.e. pyfits.open('foo.fits').data will NEVER work, use
# pyfits.open('foo.fits')[0].data.

# Final note v3: For reasons that I'm sure are very smart and complicated,
# Pyfits and DS9 disagree about row/column order. For example, if you are
# looking at a FITS file in DS9 and find a cool feature at (x,y) = (100,200)
# you access that pixel in Python via fitsdata[200,100].
#
############ Writing
#
# Just like Pyfits reads in lists of HDUs, it also writes lists of HDUs. If
# all you want to do is just dump some data to a fits file, the easiest way is:

pyfits.PrimaryHDU(fitsdata2).writeto('fitsoutput.fits')

# What you're actually doing here is first creating a HDU with
# pyfits.PrimaryHDU() and then calling its writeto() method. You can break up
# these steps if you want to, for example, add some header values:

out_HDU = pyfits.PrimaryHDU(fitsdata2)
out_HDU.header.update('DANK',True,comment='Is this file dank or not?')
out_HDU.writeto('fitsoutput.fits',clobber=True)

# There are a lot more you can do with fits files, and the Pyfits Handbook
# gives a clear description of most of them. Keep in mind that the newest
# version of pyfits is v.3.x and has some very significant changes over the UW
# astro dept's v.2.x.


########################
#       PLOTTING       #
########################
#
# Ok, last thing we should look at is how to plot all of that awesome
# data you just imported. To do this we use matplotlib.pyplot. We
# already imported this way at the beginning of this file as plt,
# which is the pyplot Class. There are a few different ways to plot
# things, and you will probably find examples that use all of
# them. Trust me when I say that the way outlined below is absolutely
# the best way to do it. It might seem like a lot of extra work now,
# but it will save you many many headaches in the future if you keep
# your plotting paradigm consistent, and this way allows for the most
# flexibility. Don't be fooled by the apparent easiness of other
# methods; they are the path to sin and frustration!

# So what is The Plotting Way? First we make a figure:

fig = plt.figure()
# which is a pyplot Figure class. Now we're going to create a
# pyplot.axis class that will be attached to fig:

ax = fig.add_subplot(111)
# add_subplot() adds, of all things, a subplot to fig. The 111 option
# tells fig that this particular subplot will take up the entire
# figure area, but this doesn't have to be the case. ax =
# add_subplot(221) will make a 2x2 grid of subplots in fig and assign
# ax to the first one (the upper left), ax = add_subplot(223) will
# assign ax to the 3rd plot (lower left), etc.
#
# And now we're ready to plot!

x = np.linspace(0,100,1000) # Kind of like arange. Give me 1000
                            # numbers evenly spaced between 0 and 100
y1 = np.exp((-1*(x-20)**2)/20)
y2 = x**0.5

ax.plot(x,y1)
ax.plot(x,y2,'r:') # plot it in red (r) with a dotted line (:)

fig.show() # This is what actually draws the plot
#
# See, that wasn't so bad. Some people might tell you to plot from a
# figure with fig.plot(), and others might (gasp!) tell you to plot
# straight from pyplot with plt.plot(). Both of these methods will
# work, but I strongly discourage them because 1). they do not allow
# for as many options as plotting from an axis, 2).  everything you
# can do with plt.plot and fig.plot you can do with ax.plot, and
# 3). keeping your plotting syntax the same will make your code easier
# to write, easier to read, and totally sexy.
#
# There are a billion + 1 different options to all the plotting
# methods called above, and a billion + 2 ways to customize your
# plots. If you want to do something with a plot and don't know how to
# do it, the first step is to check out the excellent Matplotlib
# gallery at http://matplotlib.org/gallery.html They have a lot of
# cool plots there with code showing you how they are produced (they
# do use plt.plot and fig.plot a lot, though). After that, check Stack
# Overflow or ask the UW Astro Python user's group.



# A. Eigenbrot, Madison, 2013
