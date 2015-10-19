#include <stdio.h>
#include <sys/ioctl.h>
#include <cfloat>
#include <climits>
#include <random>
using namespace std;

class Player;
class Board;
float floatMax(float[], int);
float floatMin(float[], int);
int floatWhere(float[], int, float);
float *floatCopy(float[], int);
int floatIn(float[], int);

int intMax(int[], int);
int intMin(int[], int);
int *intCopy(int[], int);

float floatMax(float *list, int len) {
    
    float max = FLT_MIN;
    for (int i = 0; i < len; i++) {
	if (list[i] > max)
	    max = list[i];
    };
    
    return max;
};

float floatMin(float *list, int len) {
    
    float min = FLT_MAX;
    for (int i = 0; i < len; i++) {
	if (list[i] < min)
	    min = list[i];
    };

    return min;
};

int floatWhere(float *list, int len, float elem) {

    if (len == 1)
	return 0;
    else if (list[0] == elem)
	return 0;
    else
	return floatWhere(list + 1, len - 1, elem) + 1;
};

float *floatCopy(float *list, int len) {

    float *copy = new float[len];
    
    for (int i = 0; i < len; i++) {
	copy[i] = list[i];
    };
    return copy;
};

int floatIn(float *list, int len, float elem) {
    
    for (int i = 0; i < len; i++) {
	if (list[i] == elem)
	    return 1;
    };
    
    return 0;
};

int intMax(int *list, int len) {
    
    int max = INT_MIN;
    for (int i = 0; i < len; i++) {
	if (list[i] > max)
	    max = list[i];
    };
    
    return max;
};

int intMin(int *list, int len) {
    
    int min = INT_MAX;
    for (int i = 0; i < len; i++) {
	if (list[i] < min)
	    min = list[i];
    };

    return min;
};
	
int *intCopy(int *list, int len) {

    int *copy = new int[len];
    
    for (int i = 0; i < len; i++) {
	copy[i] = list[i];
    };
    return copy;
};


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
    void playGame(Player, Player);
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

int Board::gameOver() {
    
    if (winsFor('X') || winsFor('O') || isFull())
	return 1;
    else
	return 0;
};

/*##################################################

PLAYER CLASS

##################################################*/

class Player {
public:
    char ox, tbt;
    mt19937 mt;
    int ply;
    Player(char, char, int);
    void repr();
    char oppChar();
    float scoreOneBoard(Board *);
    int tiebreakMove(float[], int);
    int findHi(float[]);
    float *scoresFor(Board *);
    int nextMove(Board *);
};

Player::Player(char c, char tb, int pl) {  
    ox = c;
    tbt = tb;
    ply = pl;
};

void Player::repr() {
    
    printf("Player for %c\n\twith tiebreak %c\n\tand ply = %d\n",ox,tbt,ply);
    
};

char Player::oppChar() {
    
    if (ox == 'X')
	return 'O';
    else
	return 'X';
};

float Player::scoreOneBoard(Board * b) {

    if (b->winsFor(ox)) 
	return 100.0;
    if (b->winsFor(oppChar())) 
	return 0.0;
    else 
	return 50.0;
};

int Player::tiebreakMove(float *scores, int len) {

    float *tmp = floatCopy(scores, len);
    int *t = new int[len];
    float m = floatMax(scores, len);
    int i = 0;
    int h;
    int out;
    
    while (floatIn(tmp,len,m)) {
	h = floatWhere(tmp, len, m);
	t[i] = h;
	tmp[h] = -99.9;
	i += 1;
    }
    
    if (tbt == 'l')
	out = intMin(t,i);
    else if (tbt == 'r')
	out = intMax(t,i);
    else {
	mt();
	uniform_int_distribution<int> R(0,i-1);
	out = t[R(mt)];
    };

    delete[] t;
    delete[] tmp;
    
    return out;
};

/*##############################*/

float *Player::scoresFor(Board * b) {
    
    if (ply <= 0) {
	float *Z = new float[b->col];
	for (int c = 0; c < b->col; c++) {
	    if (b->allowsMove(c))
		Z[c]  = scoreOneBoard(b);
	    else
		Z[c] = -1.0;
	};
	return Z;
    }
    
    else {
	float *L = new float[b->col];
	for (int c = 0; c < b->col; c++) {
	    if (b->allowsMove(c)) {
	        b->addMove(c,ox);
		
		if (b->winsFor(ox)) 
		    L[c] = 100.0;
		
		else if (b->winsFor(oppChar())) 
		    L[c] = 0.0;
		
		else if (b->isFull()) 
		    L[c] = 50.0;
		
		else {
		    Player badguy (oppChar(),tbt,ply-1);
		    float *badscore = badguy.scoresFor(b);
		    L[c] = 100.0 - floatMax(badscore,b->col);
		    delete[] badscore;
		}
		b->delMove(c);
	    }
	    else 
		L[c] = -1.0;
	};
	
	return L;
    }
};

int Player::nextMove(Board * b) {
    
    if (ply < 0) {
	int x;
	printf("%c's move: ",ox);
	scanf("%d", &x);
	return x;
    }
    
    else 
	return tiebreakMove(scoresFor(b), b->col);
    
};

void Board::playGame(Player p1, Player p2) {
    
    int moveCol;
    repr();
    
    while (1) {
	
	//P1
	if (p1.ply >= 0 & p2.ply < 0)
	    printf("Thinking...\n");
	moveCol = p1.nextMove(this);
	while (allowsMove(moveCol) == 0) {
	    printf("That column is either full or non-existant.\nPlease try again\n");
	    moveCol = p1.nextMove(this);
	};
	addMove(moveCol,p1.ox);
	if (p1.ply >= 0 & p2.ply < 0)
	    printf("I choose %d\n",moveCol);		  
	repr();
	
	if (winsFor(p1.ox)) {
	    printf("%c wins!\n",p1.ox);
	    return;
	};
	
	//P2
	if (p2.ply >= 0 & p1.ply < 0)
	    printf("Thinking...\n");
	moveCol = p2.nextMove(this);
	while (allowsMove(moveCol) == 0) {
	    printf("That column is either full or non-existant.\nPlease try again\n");
	    moveCol = p2.nextMove(this);
	};

	addMove(moveCol,p2.ox);	
	if (p2.ply >= 0 & p1.ply < 0)
	    printf("I choose %d\n",moveCol);
	repr();
	
	if (winsFor(p2.ox)) {
	    printf("%c wins!\n",p2.ox);
	    return;
	};	
	
	if (isFull()) {
	    printf("Tie\n");
	    return;
	};
    };
};


int main () {
    
    struct winsize w;
    int r;
    int c;
    int p1 = 99;
    int p2 = 99;
    int ply1;
    int ply2;

    ioctl(0, TIOCGWINSZ, &w);

    for (int cc=0; cc < w.ws_col; cc++){printf("*");};
    printf("\n");
    int cp = (w.ws_col - 25) / 2;
    
    for (int cc=0; cc < cp; cc++){printf("*");};
    printf(" WELCOME TO CONNECT FOUR ");
    for (int cc=0; cc < cp; cc++){printf("*");};
    if ((w.ws_col - 25) % 2 == 1)
	printf("*\n");
    else
	printf("\n");

    printf("What size board do you want?\n(Rows): ");
    scanf("%d", &r);
    printf("(Columns): ");
    scanf("%d", &c);
    
    while (p1 != 1 & p1 != 2) {
	printf("\nWill the X player be a...\n(1) Human\n(2) Computer\n");
	scanf("%d",&p1);
    };

    if (p1 == 1)
	ply1 = -1;
    else {
	printf("How hard should the computer be?\n(0) Super Easy\n(1) Pretty Easy\n(2) A Tad Hard\n(3) Kinda Hard\n(4) Quite Thoughtful\n");
	scanf("%d",&ply1);
    };

    while (p2 != 1 & p2 != 2) {
	printf("\nWill the O player be a...\n(1) Human\n(2) Computer\n");
	scanf("%d",&p2);
    };

    if (p2 == 1)
	ply2 = -1;
    else {
	printf("How hard should the computer be?\n(0) Super Easy\n(1) Pretty Easy\n(2) A Tad Hard\n(3) Kinda Hard\n(4) Quite Thoughtful\n");
	scanf("%d",&ply2);
    };

    Board B(r,c);
    
    Player P1('X','?',ply1);
    Player P2('O','?',ply2);
    
    P1.repr();
    P2.repr();

    B.playGame(P1,P2);

    delete[] B.data;

    return 0;
};
