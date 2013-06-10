#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;

float fast(float x) {
  float xhalf = 0.5f * x;
  int i = *(int*)&x;
  i = 0x5f3759df - (i >> 1);
  x = *(float*)&i;
  x = x*(1.5f - (xhalf*x*x));
  return x;
}

int main(int argc, char * argv[]){
    
    int count = argc - 2;
    char * output = argv[1];
    ofstream fout;

    fout.open(output);
    fout<<setw(10)<<"Input"<<setw(10)<<"Output"<<endl;
    fout<<setiosflags(ios::fixed);

    for (int i = 0; i < count; i++) {
	float input = atof(argv[i+2]);
	fout<<setprecision(2)<<setw(10)<<input
	    <<setprecision(7)<<setw(10)<<fast(input)<<endl;
    }
    fout.close();

    return 0;
}
      

// int main(int argc, char *argv[]) {

//     int count = argc - 2;
//     char * output = argv[1];
//     FILE * ofile;
    
//     ofile = fopen(output,"w");
//     fprintf(ofile, "%10s%10s\n", "Input","Output");

//     for (int i = 0; i < count; i++) {
// 	float input = atof(argv[i+2]);
// 	fprintf(ofile,"%10.2f%10.7f\n",input,fast(input));
//     }

//     fclose(ofile);
//     return 0;
// }

// int main(int argc, char * argv[]) {
    
//     float x;
//     char * input;

//     printf("Enter you want to find the inv. sqrt. of ('q' to quit)\n");
//     scanf("%f", &x);
    
//     while (x) {
// //	x = (float) atof(input);
// //	printf("%10.2f\n",x);
// 	printf("%10.2f\n",fast(x));
// 	scanf("%f", &x);
//     }

//     return 0;
// }
