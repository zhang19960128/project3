#include<iostream>
#include "parameter.h"
#include "atom.h"
#include <vector>
#include <list>
#include <fstream>
int main(){
    int size=20;
    int N=size*size;
    std::vector<atom> atomall(N);
    int temp;
    double r_verlet=1.2*r_cut;
    for(size_t i=0;i<size;i++)
     for(size_t j=0;j<size;j++){
        temp=i*size+j;
        atomall[temp].setx(i*r_min);
        atomall[temp].sety(j*r_min);
        atomall[temp].setr(0.1*r_min);
       }
    updatelist(atomall,N,r_verlet);
}
