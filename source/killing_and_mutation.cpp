#include <mpi.h>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <array>
#include <string>
#include <vector>
#include <algorithm>
#include <sys/time.h>
// #include <random>
#include "mt19937-64.h"
// #include <stdio.h>

#define Lx 100       // system size
#define Ly 200
#define Lz 200
#define timer 100000 // total calculation time
#define sample 1     // sample number for taking ensemble avarages
#define gTheta 8.0   // fixed
#define gK 50.0/8.0  // fixed

// ======important parameters=====
#define killrate 0.221    //C_k
#define mutantrate 0.05   //p_m
// ===============================

static std::vector<int> I,IH; //sorted cell index, sorted position index
static std::vector<double> T,KT; //sorted dividing time
static std::vector<int> S; //spatial profile
static double t; //external time

static int face[Ly*Lz]; //colony projection to yz plane
static int initface[Ly*Lz]; //initial projection
static int facechange[Ly*Lz]; //opinion change times

static double dataFACEb[timer+1]; //surface density in x
static double dataFACEs[timer+1]; //surface density in y 
static double dataM[timer+1]; //magnatization in whole
static double datam[timer+1]; //magnetization in cross section
static double datap[timer+1]; //persistence
static double datan[timer+1]; //average opinion change times
static double datasc[timer/1000+1][Ly/2-1]; //spacial correlation function
static double dataac[timer+1]; //autocorrelation function
static double datand[timer/1000+1][1500]; //distribution of opinion change times "n"

static double data[timer+1]; //exporting array
static double data2[timer+1];
static double data3[timer+1];
static double data4[timer+1];
static double data5[timer/1000+1][Ly/2-1];
static double data6[timer+1];
static double data7[timer+1];
static double data8[timer/1000+1][1500];


int randomt(int x){ //random number
	struct timeval tv;
	gettimeofday(&tv, NULL);
	init_genrand64(tv.tv_sec + tv.tv_usec);
    return((int)x*genrand64_real2());
}

double tau(double n){// make division time (make gamma random number)
   if(n == 0){return timer+1;}
   int int_kappa;
   double frac_kappa;
   struct timeval tv;
   gettimeofday(&tv, NULL);
   init_genrand64(tv.tv_sec + tv.tv_usec);
    
   int_kappa  = (int)(gK);
   frac_kappa = (gK) - (double)int_kappa;
    
   double u,uu;
   double b,p,x_frac,x_int;
   int i;
    
   /*integer part*/
   x_int=0;
   for(i=0;i<int_kappa;i++){
       x_int+=-std::log(genrand64_real3()); // add expnential random number with mean 1
   }
    
   /*fractional part*/
   if( std::fabs(frac_kappa) < 0.01 ) x_frac=0;
 
   else{
       b=(exp(1.0)+frac_kappa)/exp(1.0);
       while(1){
        
           u=genrand64_real3();
           p=b*u;
            
           uu=genrand64_real3();
            
           if(p<=1.0){
               x_frac=std::pow(p,1.0/frac_kappa);
               if(uu<=std::exp(-x_frac)) break;
           }
            
           else{
               x_frac=-std::log((b-p)/frac_kappa);
               if(uu<=std::pow(x_frac,frac_kappa-1.0)) break;
           }
        
       }
   }
    
   return (x_int+x_frac)*(gTheta/n);
}

double tau_killold(double mu){// make T6SS killing activation time (make gaussian random number)
 	struct timeval tv;
 	gettimeofday(&tv, NULL);
 	init_genrand64(tv.tv_sec + tv.tv_usec);

    double z=sqrt( -2.0*log(genrand64_real3()) ) * sin( 2.0*M_PI*genrand64_real3() );
    return mu + 0.5*mu*z;
}

double tau_kill(double mu){
	struct timeval tv;
 	gettimeofday(&tv, NULL);
 	init_genrand64(tv.tv_sec + tv.tv_usec);

    double z = -log(genrand64_real3())*mu;
    return z;
}

void initial(double m, double killrate2){// make random initial condition
	T.resize(Lx*Ly*Lz,0);
	S.resize(Lx*Ly*Lz,0);
	I.resize(Lx*Ly*Lz,0);
	IH.resize(Lx*Ly*Lz,0);
	KT.resize(Lx*Ly*Lz,0);
	double thre = (m+1)/2*100;

	if(killrate2 != 0){
		for(int i=0; i<Ly*Lx*Lz; i++){
			I[i] = i;
			KT[i] = tau_kill(gK*gTheta/killrate2);	
		}
		std::make_heap(KT.begin(), KT.end(), std::greater<double>());		
	}
	else{
		for(int i=0; i<Ly*Lx*Lz; i++){
			I[i] = i;
			KT[i] = timer+1;
		}
	}

	std::random_shuffle(I.begin(),I.end());


	for(int i=0; i<Ly*Lx*Lz; i++){
		T[i] = tau(1);
		IH[I[i]] = i;
		S[i] = 1;
		int val = randomt(100);
		if(val >= thre){S[i]= 0;}
	}// make initial condition
	for(int i=0; i<Ly*Lz; i++){
		face[i] = 0;
	}
	std::make_heap(T.begin(), T.end(), std::greater<double>());// make random initial condition
}

void initial_column(double killrate2){
	T.resize(Lx*Ly*Lz,0);
	S.resize(Lx*Ly*Lz,0);
	I.resize(Lx*Ly*Lz,0);
	IH.resize(Lx*Ly*Lz,0);
	KT.resize(Lx*Ly*Lz,0);
	
	if(killrate2 != 0){
		for(int i=0; i<Ly*Lx*Lz; i++){
			I[i] = i;
			KT[i] = tau_kill(gK*gTheta/killrate2);	
		}
		std::make_heap(KT.begin(), KT.end(), std::greater<double>());		
	}
	else{
		for(int i=0; i<Ly*Lx*Lz; i++){
			I[i] = i;
			KT[i] = timer+1;
		}
	}
	std::random_shuffle(I.begin(),I.end());
	for(int zz=0; zz<Lz; zz++){
		for(int yy=0; yy<Ly; yy++){
			int ss;
			int sq = std::pow((yy-Ly/2),2)+std::pow((zz-Lz/2),2);
			if(sq >= std::pow(25,2)){
				ss = 0;
			}
			else{
				ss = 1;
			}
			for(int xx=0; xx<Lx; xx++){
				int ind = zz*Ly*Lx+yy*Lx+xx;
				T[ind] = tau(1);
				IH[I[ind]] = ind;
				S[ind] = ss;
			}
		}
	}
	for(int i=0; i<Ly*Lz; i++){
		face[i] = 0;
	}
	std::make_heap(T.begin(), T.end(), std::greater<double>());// make random initial condition
}

namespace myheap { //edited by shitaro san
	namespace impl {
		template<typename RandomAccessIterator, typename Distance, typename Tp, typename Tpp, typename Compare> void push_heap(RandomAccessIterator first, Distance holeIndex, Distance topIndex, Tp value, Tpp valueI, Compare& comp);
		template<typename RandomAccessIterator, typename Distance, typename Tp, typename Tpp, typename Compare> void adjust_heap(RandomAccessIterator first, Distance holeIndex, Distance len, Tp value, Tpp valueI ,Compare comp);
		template<typename RandomAccessIterator, typename Distance, typename Tp, typename Tpp, typename Compare> void push_heap2(RandomAccessIterator first, Distance holeIndex, Distance topIndex, Tp value, Tpp valueI, Compare& comp);
		template<typename RandomAccessIterator, typename Distance, typename Tp, typename Tpp, typename Compare> void adjust_heap2(RandomAccessIterator first, Distance holeIndex, Distance len, Tp value, Tpp valueI ,Compare comp);


		template<typename RandomAccessIterator, typename Distance, typename Tp, typename Tpp, typename Compare>
		void push_heap(RandomAccessIterator first, Distance holeIndex, Distance topIndex, Tp value, Tpp valueI, Compare& comp) {
			Distance parent = (holeIndex - 1) / 2;
			//Tpp indexa = I[std::distance(first, first + holeIndex)];//??

			//while (holeIndex > topIndex && comp(*(first + parent), value)) {
			while (holeIndex > 0 && comp(*(first + parent), value)) {
				*(first + holeIndex) = std::move(*(first + parent));

				//I[std::distance(first, first + holeIndex)] = I[std::distance(first,first + parent)];//
				I[std::distance(first, first + holeIndex)] = I[std::distance(first,first + parent)];//
				IH[I[std::distance(first, first + holeIndex)]] = std::distance(first, first + holeIndex);//

				holeIndex = parent;
				parent = (holeIndex - 1) / 2;
			}
			*(first + holeIndex) = std::move(value);
			I[std::distance(first, first + holeIndex)] = valueI;//
			IH[valueI] = std::distance(first, first + holeIndex);//

		}
		template<typename RandomAccessIterator, typename Distance, typename Tp, typename Tpp, typename Compare>
		void push_heap2(RandomAccessIterator first, Distance holeIndex, Distance topIndex, Tp value, Tpp valueI, Compare& comp) {
			Distance parent = (holeIndex - 1) / 2;
			//Tpp indexa = I[std::distance(first, first + holeIndex)];//??

			//while (holeIndex > topIndex && comp(*(first + parent), value)) {
			while (holeIndex > 0 && comp(*(first + parent), value)) {
				*(first + holeIndex) = std::move(*(first + parent));
				holeIndex = parent;
				parent = (holeIndex - 1) / 2;
			}
			*(first + holeIndex) = std::move(value);
		}

		template<typename RandomAccessIterator, typename Distance, typename Tp, typename Tpp, typename Compare>
		void adjust_heap(RandomAccessIterator first, Distance holeIndex, Distance len, Tp value, Tpp valueI, Compare comp) {
			const Distance topIndex = holeIndex;
			Distance secondChild = holeIndex;

			while (secondChild < (len - 1) / 2) {
				secondChild = 2 * (secondChild + 1);
				if ( comp(*(first + secondChild), *(first + (secondChild - 1)) )) {
					--secondChild;
				}
				*(first + holeIndex) = std::move(*(first + secondChild));

				//I[std::distance(first, first + holeIndex)] = I[std::distance(first,first + secondChild)];//
				I[std::distance(first, first + holeIndex)] = I[std::distance(first,first + secondChild)];//
				IH[I[std::distance(first, first + holeIndex)]] = std::distance(first, first + holeIndex);//

				holeIndex = secondChild;
			}

			if ((len & 1) == 0 && secondChild == (len - 2) / 2) {// if l is odd..
				secondChild = 2 * (secondChild + 1);
				*(first + holeIndex) = std::move(*(first + (secondChild - 1)));
				// *(first + holeIndex) = std::move(*(first + (secondChild)));
				
				//I[std::distance(first, first + holeIndex)] = I[std::distance(first,first + (secondChild-1))];//
				I[std::distance(first, first + holeIndex)] = I[std::distance(first,first + secondChild-1)];//
				// I[std::distance(first, first + holeIndex)] = I[std::distance(first,first + secondChild)];//
				IH[I[std::distance(first, first + holeIndex)]] = std::distance(first, first + holeIndex);//

				holeIndex = secondChild - 1;
				// holeIndex = secondChild;
			}
			push_heap(first, holeIndex, topIndex, value, valueI, comp);
		}
		template<typename RandomAccessIterator, typename Distance, typename Tp, typename Tpp, typename Compare>
		void adjust_heap2(RandomAccessIterator first, Distance holeIndex, Distance len, Tp value, Tpp valueI, Compare comp) {
			const Distance topIndex = holeIndex;
			Distance secondChild = holeIndex;

			while (secondChild < (len - 1) / 2) {
				secondChild = 2 * (secondChild + 1);
				if ( comp(*(first + secondChild), *(first + (secondChild - 1)) )) {
					--secondChild;
				}
				*(first + holeIndex) = std::move(*(first + secondChild));

				holeIndex = secondChild;
			}

			if ((len & 1) == 0 && secondChild == (len - 2) / 2) {// if l is odd..
				secondChild = 2 * (secondChild + 1);
				*(first + holeIndex) = std::move(*(first + (secondChild - 1)));
				// *(first + holeIndex) = std::move(*(first + (secondChild)));
				
				holeIndex = secondChild - 1;
				// holeIndex = secondChild;
			}
			push_heap2(first, holeIndex, topIndex, value, valueI, comp);
		}
	} // namespace impl

	template<typename RandomAccessIterator, typename Compare>
	void insert_heap(RandomAccessIterator first, RandomAccessIterator last, typename std::iterator_traits<RandomAccessIterator>::difference_type insert_index, Compare comp) {
		using ValueType = typename std::iterator_traits<RandomAccessIterator>::value_type;
		using DistanceType = typename std::iterator_traits<RandomAccessIterator>::difference_type;

		if (last - first < 2) {
			return;
		}

		const DistanceType len = last - first;
		DistanceType parent = insert_index;

		ValueType value = std::move(*(first + parent));
		int valueI = I[std::distance(first, first + parent)];//

		impl::adjust_heap(first, parent, len, value, valueI, comp);
		if (parent == 0) {
			return;
		}
	}
	template<typename RandomAccessIterator, typename Compare>
	void insert_heap2(RandomAccessIterator first, RandomAccessIterator last, typename std::iterator_traits<RandomAccessIterator>::difference_type insert_index, Compare comp) {
		using ValueType = typename std::iterator_traits<RandomAccessIterator>::value_type;
		using DistanceType = typename std::iterator_traits<RandomAccessIterator>::difference_type;

		if (last - first < 2) {
			return;
		}

		const DistanceType len = last - first;
		DistanceType parent = insert_index;

		ValueType value = std::move(*(first + parent));
		int valueI = I[std::distance(first, first + parent)];//
		impl::adjust_heap2(first, parent, len, value, valueI, comp);
		if (parent == 0) {
			return;
		}
	}
} // namespace myheap


// **************************************TIME EVOLUTION*******************************//
void changeL(int minl, int flag){//devide in same row then push left direction
	int X = minl%Lx;
	int r,stat;

	int old = IH[minl-X];//old is the iterator of vanished cell in open boundary (or killed cell)

	stat = 1;
	for(int xx=X; xx>0; xx--){//find dead cell
		if(T[IH[minl-X+xx-1]] == timer+1){
			stat = xx;
			old = IH[minl-X+xx-1];
			break;
		}
	} 

	for(int xx=stat; xx<X; xx++){// slide states
		r = xx+minl-X;
		S[r-1] = S[r];
		//S[r] = S[r+1];
		I[IH[r]] -= 1;
		IH[r-1] = IH[r];
	}

	if(X != 0){// push out the ejected cell
		if(flag == 0){
			S[minl - 1] = S[minl];
		}
		else{
			S[minl - 1] = std::abs(1-S[minl]);
		}	
		I[old] = minl - 1;
		IH[minl - 1] = old;
		T[old] = tau(1) + t;
		myheap::insert_heap(T.begin(), T.end(), old, std::greater<double>());
	}

	auto a = 0;
	T[0] = tau(1) + t;
	myheap::insert_heap(T.begin(), T.end(), a, std::greater<double>()); //tau() + t
}

void changeR(int minl, int flag){//devide in same row then push right direction
	int X = minl%Lx;
	int r,stat;

	int old = IH[minl-X+Lx-1];
	
	stat = Lx-2;
	for(int xx=X; xx<Lx-1; xx++){//find dead cell
		if(T[IH[minl-X+xx+1]] == timer+1){
			stat = xx;
			old = IH[minl-X+xx+1];
			break;
		}
	} 

	for(int xx=stat; xx>X; xx--){
		r = xx+minl-X;
		S[r+1] = S[r];
		//S[r] = S[r-1];
		I[IH[r]] += 1;
		IH[r+1] = IH[r];
	}

	if(X != Lx-1){
		if(flag == 0){
			S[minl + 1] = S[minl];
		}
		else{
			S[minl + 1] = std::abs(1-S[minl]);
		}	
		I[old] = minl+1;
		IH[minl+1] = old;
		T[old] = tau(1) + t;
		myheap::insert_heap(T.begin(), T.end(), old, std::greater<double>());
	}

	auto a = 0;
	T[0] = tau(1) + t;
	myheap::insert_heap(T.begin(), T.end(), a, std::greater<double>());
}

void changeSL(int minl, int olabel, int flag){//invade from neighbor row then push left direction
	int X = minl%Lx;
	int r,stat;

	int old = IH[minl-X];//old is the iterator of vanished cell in open boundary

	stat = 1;
	for(int xx=X+1; xx>0; xx--){//find dead cell
		if(T[IH[minl-X+xx-1]] == timer+1){
			stat = xx;
			old = IH[minl-X+xx-1];
			break;
		}
	} 	

	for(int xx=stat; xx<X+1; xx++){
		r = xx+minl-X;
		S[r-1] = S[r];
		//S[r] = S[r+1];
		I[IH[r]] -= 1;
		IH[r-1] = IH[r];
	}

	if(flag == 0){
		S[minl] = S[olabel];
	}
	else{
		S[minl] = std::abs(1-S[olabel]);
	}
	// I[IH[minl]] -= 1;
	// IH[minl - 1] = IH[minl];

	I[old] = minl;
	IH[minl] = old;
	T[old] = tau(1) + t;
	myheap::insert_heap(T.begin(), T.end(), old, std::greater<double>());// renew daughter

	auto a = 0;
	T[0] = tau(1) + t;
	myheap::insert_heap(T.begin(), T.end(), a, std::greater<double>());// renew original
}

void changeSR(int minl, int olabel, int flag){//invade from neighbor row then push right direction
	int X = minl%Lx;
	int r,stat;

	int old = IH[minl-X+Lx-1];
	
	stat = Lx-2;
	for(int xx=X-1; xx<Lx-1; xx++){//find dead cell
		if(T[IH[minl-X+xx+1]] == timer+1){
			stat = xx;
			old = IH[minl-X+xx+1];
			break;
		}
	} 	

	for(int xx=stat; xx>X-1; xx--){
		r = xx+minl-X;
		S[r+1] = S[r];
		//S[r] = S[r-1];
		I[IH[r]] += 1;
		IH[r+1] = IH[r];
	}

	if(flag == 0){
		S[minl] = S[olabel];
	}
	else{
		S[minl] = std::abs(1-S[olabel]);
	}
	// I[IH[minl]] += 1;
	// IH[minl + 1] = IH[minl];

	I[old] = minl;
	IH[minl] = old;
	T[old] = tau(1) + t;
	myheap::insert_heap(T.begin(), T.end(), old, std::greater<double>()); //renew daughter cell from vanished element

	auto a = 0;
	T[0] = tau(1) + t;
	myheap::insert_heap(T.begin(), T.end(), a, std::greater<double>()); // renew original
}

void kill2(int sellabel, double killrate2){//T6SS killing event
	int X = sellabel%Lx;
	int Y = (sellabel%(Ly*Lx)-X)/Lx;
	int Z = (sellabel-sellabel%(Ly*Lx))/(Ly*Lx);

	KT[0] = t+tau_kill(gK*gTheta/killrate2);
	myheap::insert_heap2(KT.begin(), KT.end(), 0, std::greater<double>());

	int val = randomt(6);

	switch(val){
		case 0:
			if(S[sellabel] != S[Z*Ly*Lx + Y*Lx+ (X+1)%Lx]){
				T[IH[Z*Ly*Lx + Y*Lx+ (X+1)%Lx]] = timer+1;
				myheap::insert_heap(T.begin(), T.end(), IH[Z*Ly*Lx + Y*Lx+ (X+1)%Lx], std::greater<double>());
				break;
			}
			break;
		case 1:
			if(S[sellabel] != S[Z*Ly*Lx + Y*Lx+ (X+Lx-1)%Lx]){
				T[IH[Z*Ly*Lx + Y*Lx+ (X+Lx-1)%Lx]] = timer+1;
				myheap::insert_heap(T.begin(), T.end(), IH[Z*Ly*Lx + Y*Lx+ (X+Lx-1)%Lx], std::greater<double>());
				break;
			}
			break;
		case 2:
			if(S[sellabel] != S[Z*Ly*Lx + (Y+1)%Ly*Lx+ X]){
				T[IH[Z*Ly*Lx + (Y+1)%Ly*Lx+ X]] = timer+1;
				myheap::insert_heap(T.begin(), T.end(), IH[Z*Ly*Lx + (Y+1)%Ly*Lx+ X], std::greater<double>());
				break;
			}
			break;
		case 3:
			if(S[sellabel] != S[Z*Ly*Lx + (Y+Ly-1)%Ly*Lx+ X]){
				T[IH[Z*Ly*Lx + (Y+Ly-1)%Ly*Lx+ X]] = timer+1;
				myheap::insert_heap(T.begin(), T.end(), IH[Z*Ly*Lx + (Y+Ly-1)%Ly*Lx+ X], std::greater<double>());
				break;
			}
			break;
		case 4:
			if(S[sellabel] != S[(Z+1)%Lz*Ly*Lx + Y*Lx+ X]){
				T[IH[(Z+1)%Lz*Ly*Lx + Y*Lx+ X]] = timer+1;
				myheap::insert_heap(T.begin(), T.end(), IH[(Z+1)%Lz*Ly*Lx + Y*Lx+ X], std::greater<double>());
				break;
			}
			break;
		case 5:
			if(S[sellabel] != S[(Z+Lz-1)%Lz*Ly*Lx + Y*Lx+ X]){
				T[IH[(Z+Lz-1)%Lz*Ly*Lx + Y*Lx+ X]] = timer+1;
				myheap::insert_heap(T.begin(), T.end(), IH[(Z+Lz-1)%Lz*Ly*Lx + Y*Lx+ X], std::greater<double>());
				break;
			}
			break;
	}

}

void neighbor(int sellabel){// cell division event
	double val,val2,val3;

	val = randomt(1000);	//randomt(100)
	val3 = randomt(1000)/1000.0;

	if(1 <= val+1 && val+1 <= 250){ //left x-- (1 <= val+1 && val+1 <= yoko)
		if(0 <= val3 && val3 < mutantrate){
			changeL(sellabel, 1);
		}
		else{
			changeL(sellabel, 0);
		}
	}

	else if(251 <= val+1 && val+1 <= 500){ //right x++ (yoko+1 <= val+1 && val+1 <= yoko*2)
		if(0 <= val3 && val3 < mutantrate){
			changeR(sellabel, 1);
		}
		else{
			changeR(sellabel, 0);
		}
	}

	else if(501 <= val+1 && val+1 <= 625){ //foreward y-- (yoko*2+1 <= val+1 && val+1 <= yoko*2+tate/2.0)
		if(sellabel%(Ly*Lx) >= Lx){
			if(0 <= val3 && val3 < mutantrate){
				val2 = randomt(2);
				if(val2 == 0){
					changeSR(sellabel - Lx, sellabel, 1);
				}
				else{
					changeSL(sellabel - Lx, sellabel, 1);
				}
			}
			else{
				val2 = randomt(2);
				if(val2 == 0){
					changeSR(sellabel - Lx, sellabel, 0);
				}
				else{
					changeSL(sellabel - Lx, sellabel, 0);
				}			
			}
		}
		else{
			if(0 <= val3 && val3 < mutantrate){
				val2 = randomt(2);
				if(val2 == 0){
					changeSR(sellabel + Lx*(Ly-1), sellabel, 1);
				}
				else{
					changeSL(sellabel + Lx*(Ly-1), sellabel, 1);
				}
			}
			else{
				val2 = randomt(2);
				if(val2 == 0){
					changeSR(sellabel + Lx*(Ly-1), sellabel, 0);
				}
				else{
					changeSL(sellabel + Lx*(Ly-1), sellabel, 0);
				}			
			}
		}
	}

	else if(626 <= val+1 && val+1 <= 750){ //backward y++
		if(sellabel%(Ly*Lx) < Lx*(Ly-1)){
			if(0 <= val3 && val3 < mutantrate){
				val2 = randomt(2);
				if(val2 == 0){
					changeSR(sellabel + Lx, sellabel, 1);
				}
				else{
					changeSL(sellabel + Lx, sellabel, 1);
				}
			}
			else{
				val2 = randomt(2);
				if(val2 == 0){
					changeSR(sellabel + Lx, sellabel, 0);
				}
				else{
					changeSL(sellabel + Lx, sellabel, 0);
				}			
			}
		}
		else{
			if(0 <= val3 && val3 < mutantrate){
				val2 = randomt(2);
				if(val2 == 0){
					changeSR(sellabel - Lx*(Ly-1), sellabel, 1);
				}
				else{
					changeSL(sellabel - Lx*(Ly-1), sellabel, 1);
				}
			}
			else{
				val2 = randomt(2);
				if(val2 == 0){
					changeSR(sellabel - Lx*(Ly-1), sellabel, 0);
				}
				else{
					changeSL(sellabel - Lx*(Ly-1), sellabel, 0);
				}			
			}
		}	
	} 

	else if(751 <= val+1 && val+1 <= 875){ //upper z--
		if(sellabel >= Lx*Ly){
			if(0 <= val3 && val3 < mutantrate){
				val2 = randomt(2);
				if(val2 == 0){
					changeSR(sellabel - Ly*Lx, sellabel, 1);
				}
				else{
					changeSL(sellabel - Ly*Lx, sellabel, 1);
				}
			}
			else{
				val2 = randomt(2);
				if(val2 == 0){
					changeSR(sellabel - Ly*Lx, sellabel, 0);
				}
				else{
					changeSL(sellabel - Ly*Lx, sellabel, 0);
				}			
			}
		}
		else{
			if(0 <= val3 && val3 < mutantrate){
				val2 = randomt(2);
				if(val2 == 0){
					changeSR(sellabel + Lx*Ly*(Lz-1), sellabel, 1);
				}
				else{
					changeSL(sellabel + Lx*Ly*(Lz-1), sellabel, 1);
				}
			}
			else{
				val2 = randomt(2);
				if(val2 == 0){
					changeSR(sellabel + Lx*Ly*(Lz-1), sellabel, 0);
				}
				else{
					changeSL(sellabel + Lx*Ly*(Lz-1), sellabel, 0);
				}			
			}
		}
	}

	else{ //downer z++
		if(sellabel < Lx*Ly*(Lz-1)){
			if(0 <= val3 && val3 < mutantrate){
				val2 = randomt(2);
				if(val2 == 0){
					changeSR(sellabel + Ly*Lx, sellabel, 1);
				}
				else{
					changeSL(sellabel + Ly*Lx, sellabel, 1);
				}
			}
			else{
				val2 = randomt(2);
				if(val2 == 0){
					changeSR(sellabel + Ly*Lx, sellabel, 0);
				}
				else{
					changeSL(sellabel + Ly*Lx, sellabel, 0);
				}			
			}
		}
		else{
			if(0 <= val3 && val3 < mutantrate){
				val2 = randomt(2);
				if(val2 == 0){
					changeSR(sellabel - Lx*Ly*(Lz-1), sellabel, 1);
				}
				else{
					changeSL(sellabel - Lx*Ly*(Lz-1), sellabel, 1);
				}
			}
			else{
				val2 = randomt(2);
				if(val2 == 0){
					changeSR(sellabel - Lx*Ly*(Lz-1), sellabel, 0);
				}
				else{
					changeSL(sellabel - Lx*Ly*(Lz-1), sellabel, 0);
				}			
			}
		}
	}
	//time evolution function
}
// ************************************************************************************/

// **************************************DATA EXPORT*******************************//
void make_b_density(int checker){//measure surface density horizontal to x dir. in whole system
	double surf = 0;
	for(int xx = 0; xx<Lx; xx++){//yz plane.
		for(int zz = 0; zz<Lz; zz++){
			for(int yy = 0; yy<Ly-1; yy++){
				surf += std::abs(S[zz*Ly*Lx+(yy+1)*Lx+xx]-S[zz*Ly*Lx+yy*Lx+xx]);
			}
		}
		for(int yy = 0; yy<Ly; yy++){
			for(int zz = 0; zz<Lz-1; zz++){
				surf += std::abs(S[(zz+1)*Ly*Lx+yy*Lx+xx]-S[zz*Ly*Lx+yy*Lx+xx]);
			}
		}
	}
	dataFACEs[checker-1] += surf/(Lx*(2.0*Lx*Ly-(Ly+Lx)));
	surf = 0;
	for(int yy = 0; yy<Ly; yy++){//x dir.
		for(int zz = 0; zz<Lz; zz++){
			for(int xx = 0; xx<Lx-1; xx++){
				surf += std::abs(S[zz*Ly*Lx+yy*Lx+xx+1]-S[zz*Ly*Lx+yy*Lx+xx]);
			}
		}
	}
	checker = 1;
	dataFACEb[checker-1] += surf/(Ly*(Lx-1)*Lz);//make sueface
}

void make_b_projection(int checker){
	char fnm[40];
	FILE *fp;
	sprintf(fnm,"./mov20/p%d.csv",checker);
	fp = fopen(fnm,"w");
	for(int zz=0; zz<Lz; zz++){
		for(int yy=0; yy<Ly; yy++){
			double good = 0;
			for(int xx=0; xx<Lx; xx++){
				good += S[zz*Ly*Lx+yy*Lx+xx];
			}
			good = good/Lx;
			if(good > 0.5){
				fprintf(fp,"1,");
			}
			else{
				fprintf(fp,"0,");
			}
			// fprintf(fp, "%d,", S[zz*Ly*Lx+yy*Lx+150]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
}

void make_s_projection(int checker){
	char fnm[40],flm[40],frm[40];
	FILE *ft, *fl, *fr;

	sprintf(fnm,"./mov2/t%d.csv",checker);
	ft = fopen(fnm,"w");
	for(int yy=0; yy<Ly; yy++){
		for(int xx=0; xx<Lx; xx++){
			fprintf(ft, "%d,", S[yy*Lx+xx]);
		}
		fprintf(ft, "\n");
	}
	fclose(ft);

	sprintf(flm,"./mov2/l%d.csv",checker);
	fl = fopen(flm,"w");
	for(int yy=0; yy<Ly; yy++){
		for(int xx=0; xx<Lx; xx++){
			fprintf(fl, "%d,", S[(yy*Ly+xx)*Lx]);
		}
		fprintf(fl, "\n");
	}
	fclose(fl);

	sprintf(frm,"./mov2/r%d.csv",checker);
	fr = fopen(frm,"w");
	for(int yy=0; yy<Ly; yy++){
		for(int xx=0; xx<Lx; xx++){
			fprintf(fr, "%d,", S[(yy*Ly+Ly-1)*Lx+xx]);
		}
		fprintf(fr, "\n");
	}
	fclose(fr);
}

void make_b_density_p(int checker){
	double surf = 0;
	int old;
	for(int yy=0; yy<Ly; yy++){
		for(int zz=0; zz<Lz; zz++){
			double good = 0;
			for(int xx=0; xx<Lx; xx++){
				if(T[IH[zz*Ly*Lx+yy*Lx+xx]] != timer+1){good += S[zz*Ly*Lx+yy*Lx+xx];}
			}
			good = good/Lx;
			if(good > 0.5){
				old=face[yy*Lz+zz];
				face[yy*Lz+zz] = 1;
				facechange[yy*Lz+zz] += std::abs(old-1);
			}
			else{
				old=face[yy*Lz+zz];
				face[yy*Lz+zz] = 0;
				facechange[yy*Lz+zz] += old;
			}
		}
	}
	for(int zz = 0; zz<Lz; zz++){
		for(int yy = 0; yy<Ly-1; yy++){
			surf += std::abs(face[(yy+1)*Lz+zz]-face[yy*Lz+zz]);
		}
	}
	for(int yy = 0; yy<Ly; yy++){
		for(int zz = 0; zz<Lz-1; zz++){
			surf += std::abs(face[yy*Lz+zz+1]-face[yy*Lz+zz]);
		}
	}
	dataFACEs[checker-1] += surf/(2.0*Lz*Ly-(Ly+Lz));
	surf = 0;
	// for(int yy = 0; yy<Ly; yy++){//x dir.
	// 	for(int zz = 0; zz<Lz; zz++){
	// 		for(int xx = 0; xx<Lx-1; xx++){
	// 			surf += std::abs(S[zz*Ly*Lx+yy*Lx+xx+1]-S[zz*Ly*Lx+yy*Lx+xx]);
	// 		}
	// 	}
	// }
	dataFACEb[checker-1] += surf/(Ly*(Lx-1)*Lz);
}

void make_m(int checker){
	datam[checker-1] = 0;
	for(int i=0; i<Ly*Lz; i++){
		datam[checker-1] += (face[i]-0.5)*2;
	}
	datam[checker-1] = datam[checker-1]/(Lz*Ly);
	dataM[checker-1] = 0;
	for(int i=0; i<Ly*Lz*Lx; i++){
		if(T[IH[i]] != timer+1){dataM[checker-1] += (S[i]-0.5)*2;}
	}
	dataM[checker-1] = dataM[checker-1]/(Lz*Ly*Lx);
}

void make_persistence(int checker){// measure spin flip number of times and persistence
	datap[checker-1] = 0;
	double avec = 0;
	for(int i=0; i<Ly*Lz; i++){
		if(facechange[i] == 0){
			datap[checker-1] += 1;
		}
		avec += facechange[i];
	}
	datap[checker-1] = datap[checker-1]/(Lz*Ly);
	datan[checker-1] = avec/(Lz*Ly);
}

void make_spacial_corr(int checker){// measure spatial correlation function
	double ss;
	for(int dist=1; dist<Ly/2; dist++){
		ss = 0;
		for(int yy=0; yy<Ly; yy++){
			for(int zz=0; zz<Lz; zz++){
				ss += (face[yy*Lz+zz]-0.5)*2*(face[yy*Lz+(zz+dist)%Lz]-0.5)*2;
				if(zz-dist>=0){
					ss += (face[yy*Lz+zz]-0.5)*2*(face[yy*Lz+zz-dist]-0.5)*2;
				}
				else{
					ss += (face[yy*Lz+zz]-0.5)*2*(face[(yy+1)*Lz+zz-dist]-0.5)*2;
				}
				ss += (face[yy*Lz+zz]-0.5)*2*(face[((yy+dist)%Ly)*Lz+zz]-0.5)*2;
				if(yy-dist>=0){
					ss += (face[yy*Lz+zz]-0.5)*2*(face[(yy-dist)*Lz+zz]-0.5)*2;
				}
				else{
					ss += (face[yy*Lz+zz]-0.5)*2*(face[(Ly+yy-dist)*Lz+zz]-0.5)*2;
				}
			}
		}
		ss = ss/(4*Ly*Lz);
		datasc[checker-1][dist-1] = ss-std::pow(dataM[checker-1],2);
	}
}

void make_auto_corr(int checker){// measure autocorrelation function
	double ss = 0;
	for(int i=0; i<Ly*Lz; i++){
		ss += (initface[i]-0.5)*2*(face[i]-0.5)*2;
	}
	dataac[checker-1] = ss/(Ly*Lz);
}

void make_n_distribution(int checker){// make flip number of times distribution
	for(int i=0; i<Ly*Lz; i++){
		if(facechange[i]>=1500){continue;}
		datand[checker-1][facechange[i]] += 1;
	}
}
// ********************************************************************************//

// // R U ready?? // //
int main(int argc, char *argv[]){
	struct timeval s, e;
	gettimeofday(&s, NULL);
	double killrate2;
	int myid,totalcal,threadcal,numprocs;

	//=========================================MPI parallel calc====================================================//	
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);

	for(int trial=0; trial<sample; trial++){
		int checker = 1;
		int checker2 = 1;

		killrate2 = killrate;

		t = 0;
		initial(0, killrate2); //make initial condition

		for(int j=0; j<timer+1; j++){
			dataFACEs[j] = 0;
		}

		// //*****************************MAKE SURFACE DENSITY(direct)**************************************************//
		// make_b_density(checker);

		//*****************************MAKE projection**************************************************//
		// make_b_projection(checker);

		//*****************************MAKE cube surface projection**************************************************//
		// make_s_projection(checker);

		// //NECCESARY!!*************MAKE SURFACE DENSITY (projection)************************************************//
		make_b_density_p(checker);

		// // //*****************************MAKE MAGNETIZATION (projection)************************************************//
		make_m(checker);

		// //*****************************MAKE PERSISTENCE Prob. (projection)************************************************//
		for(int i=0; i<Lz*Ly; i++){
			facechange[i] = 0;
		}
		make_persistence(checker);

		// //*****************************MAKE AUTO CORRELATION (projection)************************************************//
		for(int i=0; i<Ly*Lz; i++){
			initface[i] = face[i];
		}
		make_auto_corr(checker);

		// //*****************************MAKE SPACIAL CORRELATION (projection)************************************************//
		make_spacial_corr(checker2);
		
		// //*****************************MAKE N(OPINION CHANGE TIMES) DISTRIBUTION (projection)*******************************//
		make_n_distribution(checker2);

		while(t <timer){
			if(T[0]<KT[0]){
				t = T[0];
			}
			else{
				t = KT[0];	
			}
						

			while(t > checker){
				checker += 1;
			// 	// //*****************************MAKE SURFACE DENSITY(direct)**************************************************//
				// make_b_density(checker);

			// 	// *****************************MAKE projection**************************************************//
				// make_b_projection(checker);

			// 	// *****************************MAKE cube surface projection**************************************************//
			// 	// make_s_projection(checker);

			// *****************************MAKE SURFACE DENSITY (projection)************************************************//
				make_b_density_p(checker);

			// 	// //*****************************MAKE MAGNETIZATION (projection)************************************************//
				make_m(checker);

			// 	// //*****************************MAKE PERSISTENCE Prob. (projection)************************************************//
				make_persistence(checker);

			// 	// //*****************************MAKE AUTO CORRELATION (projection)************************************************//
				make_auto_corr(checker);
			}

			while(t > checker2*1000){
				checker2 += 1;
				//*****************************MAKE SPACIAL CORRELATION (projection)************************************************//
				make_spacial_corr(checker2);

				//*****************************MAKE N(OPINION CHANGE TIMES) DISTRIBUTION (projection)*******************************//
				make_n_distribution(checker2);

			// }
			if(T[0]<KT[0]){
				neighbor(I[0]);	
			}
			else{
				kill2(randomt(Lx*Ly*Lz),killrate2);
			}

			// printf("\rNow %f/%d completed. (%d/%d) .",t,timer,i+1,sample);
			// fflush(stdout);
		}//time

		threadcal = trial+1;
	}// trial

	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Reduce(&dataFACEb, &data, timer+1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&dataFACEs, &data2, timer+1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&dataM, &data3, timer+1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&datam, &data4, timer+1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&datap, &data4, timer+1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	for(int i=0; i<timer/1000+1; i++){
		MPI_Reduce(&datasc[i], &data5[i], Ly/2-1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);			
	}
	MPI_Reduce(&datan, &data6, timer+1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&dataac, &data7, timer+1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	for(int i=0; i<timer/1000+1; i++){
		MPI_Reduce(&datand[i], &data8[i], 500, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);			
	}

	MPI_Reduce(&threadcal, &totalcal, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	

	if(myid == 0){ // data expose

		for(int j=0; j<timer+1; j++){//make average of surface density
			data[j] = data[j]/numprocs;
			data2[j] = data2[j]/numprocs;
			data3[j] = data3[j]/numprocs;
			data4[j] = data4[j]/numprocs;
			data6[j] = data6[j]/numprocs;
			data7[j] = data7[j]/numprocs;
		}
		for(int i=0; i<timer/1000+1; i++){
			for(int j=0; j<Ly/2-1; j++){
				data5[i][j] = data5[i][j]/numprocs;
			}
			for(int j=0; j<500; j++){
				data8[i][j] = data8[i][j]/numprocs;
			}
		}

		gettimeofday(&e, NULL);
		//std::cout << "\nNow printing.." <<std::endl;


		std::ofstream ofs("./expt2256/2Dperiod_x.csv"); // x surface dendity
		for(int j=0; j<timer+1; j++){
			ofs << data[j] << std::flush;
			if(j != timer){ofs << "," << std::flush;}
		}
		ofs << std::endl;

		std::ofstream ofs2("./expt221/2Dperiod_quench_yz.csv"); // yz surface density
		for(int j=0; j<timer+1; j++){
			ofs2 << data2[j] << std::flush;
			if(j != timer){ofs2 << "," << std::flush;}
		}
		ofs2 << std::endl;

		std::ofstream ofs3("./expt223/2Dperiod_quench_Mag.csv"); // magnetization growth
		for(int j=0; j<timer+1; j++){
			ofs3 << data3[j] << std::flush;
			if(j != timer){ofs3 << "," << std::flush;}
		}
		ofs3 << std::endl;

		std::ofstream ofs4("./expt225/2Dperiod_p.csv"); // persistence
		for(int j=0; j<timer+1; j++){
			ofs4 << data4[j] << std::flush;
			if(j != timer){ofs4 << "," << std::flush;}
		}
		ofs4 << std::endl;

		std::ofstream ofs4("./expt223/2Dperiod_quench_crossmag.csv"); // persistence
		for(int j=0; j<timer+1; j++){
			ofs4 << data4[j] << std::flush;
			if(j != timer){ofs4 << "," << std::flush;}
		}
		ofs4 << std::endl;		

		std::ofstream ofs5("./expt221/2Dperiod_SpacialCorr.csv"); // spacial correlation function
		for(int i=0; i<timer/1000+1; i++){
			for(int j=0; j<Ly/2-1; j++){
				ofs5 << data5[i][j] << std::flush;
				if(j != Ly/2-2){ofs5 << "," << std::flush;}
			}
			ofs5 << std::endl;
		}

		std::ofstream ofs6("./expt225/2Dperiod_n.csv"); // averaged opinion change times
		for(int j=0; j<timer+1; j++){
			ofs6 << data6[j] << std::flush;
			if(j != timer){ofs6 << "," << std::flush;}
		}
		ofs6 << std::endl;

		std::ofstream ofs7("./expt22543/2Dperiod_AutoCorr.csv"); // autocorrelation funtion
		for(int j=0; j<timer+1; j++){
			ofs7 << data7[j] << std::flush;
			if(j != timer){ofs7 << "," << std::flush;}
		}
		ofs7 << std::endl;

		std::ofstream ofs8("./expt225/2Dperiod_n_distribution.csv"); // distribution of opinion chnage times
		for(int i=0; i<timer/1000+1; i++){
			for(int j=0; j<1500; j++){
				ofs8 << data8[i][j] << std::flush;
				if(j != 1499){ofs8 << "," << std::flush;}
			}
			ofs8 << std::endl;
		}

		std::cout << (e.tv_sec - s.tv_sec) <<std::endl;
		std::cout << totalcal << std::endl;
	}

	MPI_Finalize();
	return 0;
	}
}