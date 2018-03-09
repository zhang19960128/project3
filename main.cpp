#include<iostream>
#include "parameter.h"
#include "atom.h"
#include <cmath>
#include <vector>
#include <list>
#include <fstream>
int main(){
    int size=20;
    int N=size*size;
    double lamda=2;
    double maxdis=0.0;
    std::vector<atom> atomall(N);
    int temp;
    double r_verlet=1.2*r_cut;
    double r_shell_initial=r_verlet-r_cut;
    double r_shell=r_shell_initial;
    for(size_t i=0;i<size;i++)
     for(size_t j=0;j<size;j++){
        temp=i*size+j;
        atomall[temp].setx(i*r_min);
        atomall[temp].sety(j*r_min);
        atomall[temp].setr(0.1*r_min);
       }
    updatelist(atomall,r_cut);
    int count=0;
    double e_end=0.0;
    double e_start=0.0;
    do{
        e_start=e_end;
        /*force has already been updated on the updateallposition*/
        maxdis=updateallposition(atomall,lamda);
        updatelist(atomall,r_cut);
        if(maxdis*2>r_shell){
            updatelist(atomall,r_verlet);
            r_shell=r_shell_initial;
        }
        else{
            r_shell=r_shell-2*maxdis;
        }
        e_end=allpotential(atomall);
        count++;
        std::cout<<count<<std::endl;
    }while(fabs(e_end-e_start)>1e-10);
    std::cout<<"There are "<<count<<" steps"<<std::endl;
    std::cout<<"The final energy is: "<<e_end<<std::endl;
    std::fstream fs;
    fs.open("atominfo.txt",std::fstream::out);
    for(size_t i=0;i<N;i++){
        fs<<atomall[i]<<std::endl;
    }
}
