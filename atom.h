#ifndef atom_h
#define atom_h
#include "parameter.h"
#include <list>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
class atom{
	public:
		atom()=default;
		atom(double x1,double x2,double r):x(x1),y(x2),radius(r){};
		void setx(double x1);
		void sety(double x2);
		void setr(double ra);
		double getx();
		double gety();
		double getr();
		void setstress_tensor(std::vector<double> a);
		std::vector<double> getstress();
		void printneighbor();
		void printstress();
        void printforce();
		void printinfo();
        void updateforce(std::vector<atom>& atomall);
		double updateposition(double);
        friend double distance(atom&,atom&);
		friend double potential(atom&,atom&);
		friend void updatelist(std::vector<atom>&,double);
        friend void updatetensor(std::vector<atom>&);
        friend std::vector<double> totaltensor(std::vector<atom>& atomall);
        friend std::vector<double> cal_force(atom& one,atom& two);
        /*updateallposition include unpdate the force in this atom*/
        friend double updateallposition(std::vector<atom>& atomall,double lamda);
		friend std::vector<double> str_tensor(atom&,atom&);
		friend std::ostream& operator<<(std::ostream& os,atom& output);
		friend std::fstream& operator<<(std::fstream& fs,atom& output);
        friend double allpotential(std::vector<atom>& allatom);
	private:
		double x;
		double y;
		double radius;
		std::list<int> neighbor;
        std::vector<double> force;
        std::vector<double> stresstensor;
};

int count(std::vector<atom> all,atom& input,double r);
void print_radial_dis(double,double,std::vector<atom>&,std::string);
std::vector<double>& operator +=(std::vector<double>& one,std::vector<double>& two);
#endif
