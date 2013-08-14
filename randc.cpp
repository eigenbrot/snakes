#include <stdio.h>
#include <stdlib.h>
#include <math.h>
using namespace std;

float rand_float(void);                //makes a random float b/t 0 and 1
void fill_array(float [], int);        //populates array with floats
float *fast_array_sqrt(float [], int); //computes inverse sqrt quickly
float *slow_array_sqrt(float [], int); //computes inverse sqrt slowely
float fast(float);                     //computes inverse sqrt on single number

int main(int argc, char **argv) {

    int size;
    float *sqrt_arr;
    char *sizearg = argv[2];
    char method = *argv[1];
    size = atoi(sizearg);

    float *test;
    test = new float[size];
    
    fill_array(test,size);
    if (method == 'f') sqrt_arr = fast_array_sqrt(test,size);
    else if (method == 's') sqrt_arr = slow_array_sqrt(test,size);
    
    // for (int i = 0; i < size; i++) {
    // 	printf("%9.7f - %9.7f\n",test[i],sqrt_arr[i]);
    // }

    delete[] test;
    delete[] sqrt_arr;

    return 0;
}

float rand_float(void) {
    //Returns a random floating point between 0 and 1

    float r = (float)rand()/(float)RAND_MAX;
    
    return r;
}

void fill_array(float *vec, int size) {
    //Fills an array with random floating point numbers between 0 and 1

    for (int i = 0; i < size; i++) {
	vec[i] = rand_float();
    }
}

float *fast_array_sqrt(float *vec, int size) {
    //Takes an array of floats and returns an array with the inverse square
    //root of each element
    
    float *result;
    result = new float[size];
    
    for (int i = 0; i < size; i++) {
	result[i] = fast(vec[i]);
    }

    return result;
}   

float *slow_array_sqrt(float *vec, int size) {
    //Takes in an array of floats and computes the inverse square root,
    //slowly

    float *result;
    result = new float[size];
    
    for (int i = 0; i < size; i++) {
	result[i] = (float)pow(vec[i],-0.5);
    }

    return result;
}

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
