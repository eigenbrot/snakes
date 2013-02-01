#! /usr/bin/env python

import random, os, sys
        
class Board:
    '''The Board class holds all the information about a current game. It takes 
    palyer moves and records them on the game field. It also displays the 
    current state of the game field'''
    
    def __init__(self,numRows,numCols):
        """Initializes the data used in connect four"""
        self.rows = numRows
        self.cols = numCols
        
        'set up a blank game board (all blanks)'
        d=[]
        for i in range(numRows):
            t=[]
            for j in range(numCols):
                t+=[' ']
            d+=[t]
        self.data = d
        
    def __repr__(self):
        """Prints the game board"""
        b=''
        for i in range(len(self.data)):
            for j in range(len(self.data[0])):
                b+='|'
                b+=self.data[i][j]
            b+='|\n'
        b+='-'*(self.cols*2 +1) + '\n'
        for i in range(self.cols):
            b+=' '+str(i%10)
        return b+'\n'
    
    def addMove(self,c,ox):
        """Adds a piece (either O or X) to the column c"""
        r=self.rows-1
        while self.data[r][c] != ' ':
            r+=-1
        self.data[r][c]=ox
        
    def allowsMove(self,c):
        """Czechs to see if moving in column c is legal. Illegalities
            arise when the column is full or non-existant.
        """
        if c not in range(self.cols) or self.data[0][c] != ' ':
            return False
        else:
            return True
            
    def isFull(self):
        """Czechs to see if the board is full of czechers. Returns true
            if the board is full.
        """
        for i in range(self.cols):
            if self.allowsMove(i) == True:
                return False
        return True
        
    def delMove(self,c):
        """Deletes a czecher from column c."""
        r=0
        while self.data[r][c] == ' ':
            r+=1
            if r not in range(self.rows): return
        self.data[r][c]=' '
        
    def winsFor(self,ox):
        """Returns true if player ox (either 'O' or 'X') has four pieces
            in a row.
        """
        for i in range(self.rows):
            if self.horz(ox,i):return True
        for j in range(self.cols):
            if self.vert(ox,j):return True
        for i in range(self.rows):
            for j in range(self.cols):
                if self.agywag(ox,i,j): return True
                if self.vagyvag(ox,i,j):return True
        return False
    
    def vert(self,ox,c):
        """Czechs for a vertical win."""
        z=0
        i=0
        while i in range(self.rows):
            if self.data[i][c]==ox:
                z+=1
            if self.data[i][c]!=ox:
                z=0
            if z>=4:return True
            i+=1
        return False
        
    def horz(self,ox,r):
        """Czechs for a horizontal win."""
        z=0
        i=0
        while i in range(self.cols):
            if self.data[r][i]==ox:
                z+=1
            if self.data[r][i]!=ox:
                z=0
            if z>=4:return True
            i+=1
        return False
        
    def agywag(self,ox,i,j):
        """Czechs for a diagonal win."""
        z=0
        r=i
        c=j
        while r in range(self.rows):
            while c in range(self.cols):
                if self.data[r][c]==ox:
                    z+=1
                if self.data[r][c]!=ox:
                    z=0
                if z>=4: return True
                else:break
            c+=1
            r+=1
        return False
        
    def vagyvag(self,ox,i,j):
        """Czechs for a diagonal win, but differently"""
        z=0
        r=i
        c=j
        while r in range(self.rows):
            while c in range(self.cols):
                if self.data[r][c]==ox:
                    z+=1
                if self.data[r][c]!=ox:
                    z=0
                if z>=4: return True
                else:break
            c+=-1
            r+=1
        return False
    
    def gameOver(self):
        """Returns true if the board represents a completed game of 
            connect 4. Completion occurs with victory or a full board.
        """
        if self.winsFor('X')==True or self.winsFor('O')==True:
            return True
        if self.isFull()==True: return True
        else:return False

    def playGame(self,p1,p2):
        """Sets up a game of connect four b/t two players, p1 and p2.
            Each player must be a Player class unless the player's
            ply=='HUMAN', in which case that player will allow for user
            input.
        """
        print self
        c=1
        while True:

            while c==1:
                if p2.ply=='HUMAN' and p1.ply!='HUMAN':
                    print 'Thinking...\n'
                moveCol=p1.nextMove(self)

                while self.allowsMove(moveCol)==False:
                    print "That column is either full or non-existant. "\
                        +"Please try again."
                    moveCol=p1.nextMove(self)

                self.addMove(moveCol,p1.ox)

                print self
                if p2.ply=='HUMAN' and p1.ply!='HUMAN':
                    print 'I choose '+str(moveCol)

                if self.winsFor(p1.ox):
                    print p1.ox+" wins!!"
                    return
                c=2

            while c==2:
                if p1.ply=='HUMAN' and p2.ply!='HUMAN':
                    print 'Thinking...\n'
                moveCol=p2.nextMove(self)

                while self.allowsMove(moveCol)==False:
                    print "That column is either full or non-existant."\
                        +"Please try again."
                    moveCol=p2.nextMove(self)

                self.addMove(moveCol,p2.ox)

                print self
                if p1.ply=='HUMAN' and p2.ply!='HUMAN':
                    print 'I choose '+str(moveCol)

                if self.winsFor(p2.ox):
                    print p2.ox+" wins!!"
                    return
                c=0

            while c==0:
                if self.isFull():
                    print 'Tie'
                    return 
                else:c=1
            
                
                
class Player:
    
    def __init__(self,ox,tbt,ply,totalply):
        """initializes the data for the player class"""
        self.ox=ox
        self.tbt=tbt
        self.ply=ply
        self.totalply=totalply
        
    def __repr__(self):
        """Displays the characteristics of the player class"""
        print "Player for "+self.ox
        print "  with tiebreak: "+self.tbt
        print "  and ply == "+str(self.ply)
        return ''
        
    def oppChar(self):
        """Returns the oponent's connect four character ('O' or 'X')"""
        if self.ox == 'X': return 'O'
        else: return 'X'
        
    def scoreOneBoard(self,b):
        """Takes a board, b, as input and returns 100 if the player won,
            0.0 if the opponent won, and 50.0 if no one has won.
        """
        if b.winsFor(self.ox):return 100.0
        if b.winsFor(self.oppChar()):return 0.0
        else: return 50.0
        
    def tiebreakMove(self,scores):
        """Takes in a list of colum scores and moves in the column with the
            highest score. If there are more than one high score then 
            the player will chose one depending on the tbt type defined in 
            the class.
        """
        t=[]
        m=max(scores)
        while m in scores:
            '''scores is changed every iteration'''
            h=self.findhi(scores)
            t.append(h)
            scores[h]=-0.5
        
        '''if there is only one high score then len(t) == 1'''
        if self.tbt =='LEFT':
            return min(t)
        if self.tbt =='RIGHT':
            return max(t)
        if self.tbt =='RANDOM':
            return random.choice(t)
            
            
    def findhi(self,H):
        """Finds the list index of the element with the highest value.
            Input: ([List])
        """
        if len(H) == 1:
            return 0
        elif max(H) == H[0]:
            return 0
        else:
            return self.findhi(H[1:]) + 1
        
    def scoresFor(self,b):
        """Returns a list of goodness scores with each element in the list
            coresponding to a column on the game board.
        """
        '''base case; just score the board'''
        if self.ply==0:
            Z=[]
            for i in range(b.cols):
                if b.allowsMove(i):
                    Z.append(self.scoreOneBoard(b))
                else:
                    Z.append(-1.0)
            return Z
        
        else:
            '''we need to recurse more until we run out of ply-ness'''
            L=[]
            for i in range(b.cols):
                if b.allowsMove(i)==True:
                    b.addMove(i,self.ox)

                    if b.gameOver()==False:
                        '''the next few lines are the heart of the AI.
                        Understand them and you will understand your opponent.
                        Understand your opponent and you will understand
                        yourself.'''
                        badguy = Player(self.oppChar(),self.tbt,
                                      self.ply-1,self.totalply)
                        
                        '''weight each score by how many moves it took
                        to get there'''
                        if self.ply - 1 > 100:
                            badscore = [s * (float(self.ply - 1)\
                                                 /(self.totalply))
                                        for s in badguy.scoresFor(b)]
                        else:
                            badscore = badguy.scoresFor(b)

#                        print self.ply,i,badscore
                        L.append(100.0-max(badscore))

                    else: L.append(100.0)
                    
                    '''we need to remove any recursive plays'''
                    b.delMove(i)

                else: L.append(-1.0)

#            if self.ply == self.totalply:
#            print "at ply {0:n} L is {1}".format(self.ply,L)
            return L
            
    def nextMove(self,b):
        """Creates the next move for the computer player, or asks the user
            for input.
        """
        if self.ply!='HUMAN':
            return self.tiebreakMove(self.scoresFor(b))
        else:
            x=input(self.ox+"'s move: ")
            return x

def main():

    '''get information about the size of the terminal window'''
    rows, columns = os.popen('stty size', 'r').read().split()
    os.system('clear')
    
    print '*'*int(columns)+str('{:*^'+columns+'}')\
        .format(' WELCOME TO CONNECT FOUR ')+'*'*int(columns)\
        +'\nWhat board size do you want?'
    
    r=input('(Rows): ')
    c=input('(Columns): ')

    print "Will the X player be a...\n(1) Human\n(2) Computer"
    p1=input()
    
    while p1 != 1 and p1 != 2:
        print "Sorry, {} is not yet a supported player type".format(
            random.choice(('Parrot','Amoeba','Gazebo','God','Turtle','Snake'))) 
        print "Will the X player be a...\n(1) Human\n(2) Computer"
        p1 = input()

    if p1 == 1: ply1='HUMAN'    
    elif p1 == 2:
        print "How hard should the computer be?\n"\
            +"(0) Super Easy\n"\
            +"(1) Pretty Easy\n"\
            +"(2) A Tad Hard\n"\
            +"(3) Kinda Hard (can take some time to think)\n"\
            +"(4) Very Thoughtful (Takes a long time to think)"
        ply1=input()


    print "Will the O player be a...\n(1) Human\n(2) Computer"
    p2=input()

    while p2 != 1 and p2 != 2:
        print "Sorry, {} is not yet a supported player type".format(
            random.choice(('Parrot','Amoeba','Gazebo','God','Turtle','Snake'))) 
        print "Will the O player be a...\n(1) Human\n(2) Computer"
        p2 = input()

    if p2 != 1:
        print "How hard should the computer be?\n"\
            +"(0) Super Easy\n"\
            +"(1) Pretty Easy\n"\
            +"(2) A Tad Hard\n"\
            +"(3) Kinda Hard (can take some time to think)\n"\
            +"(4) Very Thoughtful (Takes a long time to think)"
        ply2=input()
    else: ply2='HUMAN'

    B=Board(r,c)
    P1=Player('X','RANDOM',ply1,ply1)
    P2=Player('O','RANDOM',ply2,ply2)

    print

    B.playGame(P1,P2)
    return

if __name__ == '__main__':
    sys.exit(main())
