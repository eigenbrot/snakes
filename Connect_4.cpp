#include <stdio.h>
using namespace std;

/*##################################################

BOARD CLASS

##################################################*/

class Board {
public:
    char *data;
    int row, col;
    Board (int,int);
    void repr();
    void addMove(int, char);
    void delMove(int);
    int allowsMove(int);
    int isFull();
    int vert(char, int);
    int horz(char, int);
    int agywag(char, int, int);
    int vagyvag(char, int, int);
    int winsFor(char);
    int gameOver();
    void playGame();
};

Board::Board (int r, int c) {
    row = r;
    col = c;
    data = new char[r*c];
    for (int rr = 0; rr < row; rr++) {
	for (int cc = 0; cc < col; cc++) {
	    data[rr*col+cc] = ' ';
	};
    };
};

void Board::repr() {
    for (int rr = 0; rr < row; rr++){
	for (int cc = 0; cc < col; cc++){
	    printf("|%c",data[rr*col+cc]);
	};
	printf("|\n");
    };
    for (int cs = 0; cs < col*2+1; cs++){printf("-");};
    printf("\n");
    for (int cl = 0; cl < col; cl++){
	printf(" %d",cl);
    };
    printf("\n\n");
    return;
};

void Board::addMove(int c, char ox) {
    
    int r = row - 1;
    while (data[r*col+c] != ' ') {
	r += -1;
    };
    data[r*col+c] = ox;
};

void Board::delMove(int c) {
    
    int r = 0;
    while (data[r*col+c] == ' ') {
	r += 1;
	if (r >= row) {return;};
    };
    data[r*col+c] = ' ';
};

int Board::allowsMove(int c) {
    
    if (c < 0 || c >= col) 
	return 0;
    else if (data[0*col+c] != ' ') 
	return 0;
    else 
	return 1;
};

int Board::isFull() {
    
    for (int cc = 0; cc < col; cc++) {
	if (allowsMove(cc))
	    return 0;
    };
    return 1;
};

int Board::vert(char ox, int c) {
    
    int z = 0;
    int r = 0;
    
    while (r < row) {
	if (data[r*col+c] == ox)
	    z += 1;
	else if (data[r*col+c] != ox)
	    z = 0;
	if (z >= 4)
	    return 1;
	r += 1;
    };
    return 0;
};

int Board::horz(char ox, int r) {
    
    int z = 0;
    int c = 0;
    
    while (c < col) {
	if (data[r*col+c] == ox)
	    z += 1;
	else if (data[r*col+c] != ox)
	    z = 0;
	if (z >= 4)
	    return 1;
	c += 1;
    };
    return 0;
};

int Board::agywag(char ox, int i, int j) {
    
    int z = 0;
    int r = i;
    int c = j;
    
    while (c < col && r < row) {
	if (data[r*col+c] == ox)
	    z += 1;
	else if (data[r*col+c] != ox)
	    z = 0;
	if (z >= 4)
	    return 1;
	c += 1;
	r += 1;
    };
    return 0;
};

int Board::vagyvag(char ox, int i, int j) {
    
    int z = 0;
    int r = i;
    int c = j;
    
    while (c < col && r < row) {
	if (data[r*col+c] == ox)
	    z += 1;
	else if (data[r*col+c] != ox)
	    z = 0;
	if (z >= 4)
	    return 1;
	c += -1;
	r += 1;
    };
    return 0;
};

int Board::winsFor(char ox) {
    
    for (int rr = 0; rr < row; rr++) {
	if (horz(ox,rr))
	    return 1;
	if (agywag(ox, rr, 0))
	    return 1;
	if (vagyvag(ox, rr, col - 1))
	    return 1;
    };
    for (int cc = 0; cc < col; cc++) {
	if (vert(ox,cc))
	    return 1;
	if (agywag(ox,0,cc))
	    return 1;
	if (vagyvag(ox,0,cc))
	    return 1;
    };
    return 0;
};

/*##################################################

PLAYER CLASS

##################################################*/

class Player {
public:
    char ox, tbt;
    int ply;
    Player(char, char, int);
    void repr();
    char opChar();
    float scoreOneBoard(Board);
    int tiebreakMove(float[]);
    int findHi(float[]);
    float *scoresFor(Board);
    int nextMove(Board);
};

Player::Player(char c, char tb, int pl) {  
    ox = c;
    tbt = tb;
    ply = pl;
};

void Player::repr() {
    
    printf("Player for %c\n\twith tiebreak %c\n\tand ply = %d\n",ox,tbt,ply);

};


int main () {
    
    Board b1 (4,4);
    b1.addMove(2,'O');
    b1.addMove(2,'X');
    b1.addMove(2,'O');
    printf("Allows? %d\n",b1.allowsMove(2));
    b1.repr();
    b1.addMove(2,'X');
    b1.addMove(0,'O');
    b1.addMove(0,'X');
    b1.addMove(0,'O');
    b1.addMove(1,'X');
    b1.addMove(1,'O');
    b1.addMove(1,'X');
    b1.addMove(1,'X');
    b1.addMove(3,'O');
    b1.addMove(0,'X');
    b1.addMove(3,'O');
    b1.addMove(3,'O');
    b1.addMove(3,'O');
    b1.repr();
    printf("Allows? %d\n",b1.allowsMove(2));
    printf("full? %d\n",b1.isFull());
    printf("agy? %d\n",b1.agywag('X',0,0));
    printf("vagy? %d\n",b1.vagyvag('O',0,3));
    printf("win bot %d\n",b1.horz('O',3));
    printf("winsfor X? %d\n",b1.winsFor('X'));
    printf("winsfor O? %d\n",b1.winsFor('O'));
    Player p1 ('X','r',3);
    p1.repr();
    delete[] b1.data;
    return 0;
};
