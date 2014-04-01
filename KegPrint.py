#! /usr/bin/python

import subprocess
import os
import time

def check_log():
    
    last_line = subprocess.check_output(['tail','-1','/var/log/cups/page_log'])
    printer = last_line.split(' ')[0]
    printID = int(last_line.split(' ')[1])

    return printer, printID

def check_new(printID):

    with open(os.path.expanduser('~/.bk.log'),'r') as f:
        lines = f.readlines()
        lastID = int(lines[-1].split('\n')[0])
    
    if printID > lastID:
        return True
    else:
        return False

def do_print(printer):

    os.system('lpr -P {} -o sides=two-sided-long-edge {}'.\
                  format(printer,os.path.expanduser('~/.bk.pdf')))
    return

def update_ID(printID):

    with open(os.path.expanduser('~/.bk.log'),'a') as f:
        f.write('{}\n'.format(printID+1))

    return

def main():

    printer, printID = check_log()

    if check_new(printID):
        do_print(printer)
        update_ID(printID)
        print "printed something new at {}".format(time.asctime())
    else:
        print "nothing to print"

    return

if __name__ == '__main__':
    main()
             
