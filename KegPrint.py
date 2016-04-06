#! /usr/bin/python

import subprocess
import os
import time
import socket

logfile = os.path.expanduser('~/.kp_{}.log'.format(socket.gethostname()))

def check_log():
    
    last_line = subprocess.check_output(['tail','-1','/var/log/cups/page_log'])
    printer = last_line.split(' ')[0]
    try:
        printID = int(last_line.split(' ')[1])
    except ValueError:
        #The Cupsd.conf PageLogFormat is a little weird, we'll try one alternative
        printID = int(last_line.split(' ')[2])
        
    return printer, printID

def check_new(printID):

    with open(logfile,'r') as f:
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

    with open(logfile,'a') as f:
        f.write('{}\n'.format(printID+1))

    return

def main():

    printer, printID = check_log()

    if not os.path.exists(logfile):
        print "Initializing logfile"
        with open(logfile,'w') as f:
            f.write('{}\n'.format(printID))
            
    print "Waiting..."
    while True:

        if check_new(printID):
            do_print(printer)
            update_ID(printID)
            print "printed something new at {}".format(time.asctime())
        
        time.sleep(10)
        
    return

if __name__ == '__main__':
    main()
             
