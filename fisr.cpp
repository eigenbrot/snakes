#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;

float fast(float x) {
    // An implementation of the fast inverse square root algorithm
       
    float xhalf = 0.5f * x;
    int i = *(int*)&x;          //interpret the bits of x as an int
    i = 0x5f3759df - (i >> 1);  //Crazee constant minus half of i
    x = *(float*)&i;            //Re-interpret i as a float
    x = x*(1.5f - (xhalf*x*x)); //One iteration of Newton's method for good
                                // measure
    return x;
}

// int main(int argc, char * argv[]){
//     /* Computes the inverse square root on a bunch of inputs.

//        Usage:
//             ./a.out outputfile float1 float2 float3 ...

//        outputfile will contain all of the input numbers and their inverse
//        square roots.

//        This version of main uses c++ output streams to write the output. I'm
//        not crazy about how the formatting for this type of output works.

//     */

//     int count = argc - 2;
//     char * output = argv[1];
//     ofstream fout;

//     fout.open(output);
//     fout<<setw(10)<<"Input"<<setw(10)<<"Output"<<endl;
//     fout<<setiosflags(ios::fixed);

//     for (int i = 0; i < count; i++) {
// 	float input = atof(argv[i+2]);
// 	fout<<setprecision(2)<<setw(10)<<input
// 	    <<setprecision(7)<<setw(10)<<fast(input)<<endl;
//     }
//     fout.close();

//     return 0;
// }
      

// int main(int argc, char *argv[]) {
//     /* Computes the inverse square root on a bunch of inputs.

//        Usage:
//             ./a.out outputfile float1 float2 float3 ...

//        outputfile will contain all of the input numbers and their inverse
//        square roots.
       
//        This version of main uses normal c format strings and fprintf for the
//        output. I like it better because it's not totally stupid.

//     */

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

int main(int argc, char * argv[]) {
    /* Computes the inverse square root on a bunch of inputs.
       
       Usage:
           ./a.out
              
       This version of main does not output anything. Instead it asks the
       user to keep providing values to compute until the use inputs 0.

    */

    float x;
    char * input;

    printf("Enter you want to find the inv. sqrt. of ('q' to quit)\n");
    scanf("%f", &x);
    
    while (x) {
//	x = (float) atof(input);
//	printf("%10.2f\n",x);
	printf("%10.2f\n",fast(x));
	scanf("%f", &x);
    }

    return 0;
}
