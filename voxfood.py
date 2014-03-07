#! /usr/bin/python
import os
import sys
import pwd
import time
import pickle

enc_db = '3ywv3ywivw3imkirfvsx32zs|sqrsq32hexe32hf2hf'
db_filename = ''.join([chr(ord(i) - 4) for i in enc_db])

class VoxOmnom:

    def __init__(self, username):
        
        self.username = username
        self.mjd = int(time.time()/86400.0 + 40587.0)
        self.db_filename = db_filename
        self.command_dict = {'l': self.list_db,
                             '?': self.show_commands,
                             'q': sys.exit,
                             'w': self.vote,
                             '1': self.choice_1,
                             '2': self.choice_2,
                             '3': self.choice_3}
        self.choice_dict = {}
        f = open(self.db_filename,'rb')
        self.db = pickle.load(f)
        f.close()

        if not self.check_legality():
            print 'You already voted today!'
            self.list_db()
            sys.exit()

        self.ask_input()

    def ask_input(self):

        scratch = raw_input("\nPlease enter a command:\n('?' gives list of available commands)\n")
        if scratch.lower() not in self.command_dict.keys():
            print 'the request was made but it was not good'
            self.ask_input()
        else:
            self.command_dict[scratch.lower()]()
            self.ask_input()

    def list_db(self):
        
        try:
            current_dict = self.db[self.mjd]
        except KeyError:
            print "No votes today. Be the first!"
            return
        
        print 'Current votes:\n\t{:10}{:3}'.format('Name','num')
        print '\t'+'-'*15
        for venue in current_dict['votes'].keys():
            print '\t{:10}{:3}'.format(venue,current_dict['votes'][venue])
            
        return

    def show_commands(self):

        print "l: list today's votes\n"\
            + "?: print list of commands\n"\
            + "q: cancel voting and quit\n"\
            + "w: submit votes and quit\n"\
            + "1: vote for first choice\n"\
            + "2: vote for second choice\n"\
            + "3: vote for third choice\n"\
            + "\nYou only get one chance to vote per week, so don't mess up!"

        return

    def choice_1(self):

        venue = str(raw_input('Enter name for your first choice:\n'))
        self.choice_dict[1] = self.clean_input(venue)
        return

    def choice_2(self):

        venue = str(raw_input('Enter name for your second choice:\n'))
        self.choice_dict[2] = self.clean_input(venue)
        return

    def choice_3(self):

        venue = str(raw_input('Enter name for your third choice:\n'))
        self.choice_dict[3] = self.clean_input(venue)
        return
        
    def clean_input(self, venue):

        cleaned  = venue.\
            replace('!','').\
            replace('_',' ').\
            replace('-',' ').\
            replace('@','').\
            replace('#','').\
            replace('$','').\
            replace('%','').\
            replace('^','').\
            replace('&','').\
            replace('*','').\
            replace('(','').\
            replace(')','').\
            replace('+','').\
            replace('=','').\
            replace('{','').\
            replace('[','').\
            replace('}','').\
            replace(']','').\
            replace(':','').\
            replace(';','').\
            replace("'",'').\
            replace('"','').\
            replace('<','').\
            replace(',','').\
            replace('>','').\
            replace('.','').\
            replace('?','').\
            replace('/','').\
            replace('|','').\
            replace('~','').\
            replace('`','').\
            replace('\\','').lower()
        return cleaned

    def vote(self):
        
        print 'your current choices are:'
        for choice in self.choice_dict.keys():
            print '{}: {}'.format(choice,self.choice_dict[choice])
        
        ask = raw_input('Are you sure you want to cast your vote? (Y/n):\n')
        if ask.lower() != 'n':
            try:
                current_day = self.db[self.mjd]
            except KeyError:
                self.db[self.mjd] = {'users':[], 'votes': {}}
                current_day = self.db[self.mjd]

            for choice in self.choice_dict.keys():
                try:
                    current_day['votes'][self.choice_dict[choice].lower()] \
                        += abs(choice - 4)
                except KeyError:
                    current_day['votes'][self.choice_dict[choice].lower()] \
                        = abs(choice - 4)

            current_day['users'].append(self.username)
            f = open(self.db_filename,'wb')
            pickle.dump(self.db,f)
            f.close()
            sys.exit()
        return

    def check_legality(self):

        try:
            current_day = self.db[self.mjd]
            if self.username in current_day['users']:
                return False
        except KeyError:
            return True
        return True

def main():
    uid = os.getuid()
    user = pwd.getpwuid(uid).pw_name

    if not os.path.exists(db_filename):
        print 'Could not load the database file!'
        print 'This program must be run on Astro dept. machines'
        return

    if len(sys.argv) > 1:
        if 's' in sys.argv[1]:
            if len(sys.argv) > 2:
                mjd = sys.argv[2]
            else:
                mjd = int(time.time()/86400.0 + 40587.0)
            show_winner(mjd)
            return

        if 'h' in sys.argv[1]:
            print 'usage: voxomnom [-h|-s]\n'\
                + '\n-h: This text\n'\
                + '-s: Show current standings\n'\
                + 'If no option is given, vote!'
            return
    else:
        VoxOmnom(user)

    return

def create_db():

    d = {}
    f = open(db_filename,'wb')
    pickle.dump(d,f)
    f.close()
    os.system('chmod 777 {}'.format(db_filename))

def get_db():

    f = open(db_filename,'rb')
    d = pickle.load(f)
    f.close()
    return d

def write_db(d):
    
    f = open(db_filename,'wb')
    pickle.dump(d,f)
    f.close()

def show_winner(mjd):

    d = get_db()
    try:
        votes = d[mjd]['votes']
        standings = sorted(votes.items(), key=lambda x: x[1])[::-1]
        
        print 'Votes for {}:\n\t{:10}{:3}'.format(mjd,'Name','num')
        print '\t'+'-'*15
        for venue in standings:
            print '\t{:10}{:3}'.format(venue[0],venue[1])
    except KeyError:
        print "No votes on that date"
    return

if __name__ == '__main__':
    sys.exit(main())
