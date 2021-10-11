//============================================================================
// Name        : Montecarlo.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description :
//============================================================================

#include <iostream>
#include <iomanip>
#include <vector>
#include <chrono>
#include <random>
#include <numeric>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <omp.h>

#include "Classic_walk.hpp"
#include "Folded_walk.hpp"

using namespace std;



const size_t System_size = 400; //4*250
const size_t T = System_size /2; // 20000

const double dt =  0.1;

const double prob_to_move = pow(dt,1);

const int N_of_samples = 1'000; //50000

int main(
		//int argc, const char * argv[]
									 ) {

	if (prob_to_move > 0.5 || prob_to_move < 0){
		cout << "Probability should be in a range [0:1]" << endl;
		exit(1);
	}
    ofstream myFile, sz_strm;

    stringstream ss;
    ss << dt;
    string dt_for_filename;
    ss >> dt_for_filename;

    ios_base::openmode mode;
	mode = std::ofstream::out; //Erase previous file (if present)
	sz_strm.open("Data/Sz.dat", mode);
	//sz_strm.open("Data/Sz_dt" + dt_for_filename +".dat", mode);

	myFile.open("Data/classical_jammed.txt", mode);
	sz_strm.precision(9);
	myFile.precision(9);

	sz_strm << "# x \t sigma_z \t t" << endl;

    ////////////////////////////////////////////////////////////////
	const int Evolution_steps = T/dt;
    vector<vector<double> > contAvg(Evolution_steps, vector<double>(System_size ,0)); // vector<T> vec (how many elements, initialize by {.} )
    vector<vector<double> > contDev(contAvg);

    random_device rd;  // Will be used to obtain a seed for the random number engine
    mt19937 generator(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> distReal(0.0, 1.0);

    vector<int> posLeftMovers;
    vector<int> posRightMovers;
    posLeftMovers.reserve(System_size );
    posRightMovers.reserve(System_size );
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    auto start = chrono::steady_clock::now();
//#pragma omp parallel for num_threads(omp_get_num_procs())
    for(int count=0; count<N_of_samples; count++){

    	size_t particle_position = System_size /2;
    	vector<int> config = initialState_Folded(particle_position , System_size );






    	for (size_t n_evol = 0; n_evol < Evolution_steps; ++n_evol){


    		evolveState_Classic(config, prob_to_move);

		// we start collecting data

			for (size_t j = 0; j < System_size ; j++) {
				//contDev[n_evol][j] = (contDev[n_evol][j] + (1.0 * config[j] - contAvg[n_evol][j]) * (1.0 * config[j] - contAvg[n_evol][j])
				//				/ (count + 1)) * count / (count + 1);
				contAvg[n_evol][j] = (config[j] + 1.0 * count * contAvg[n_evol][j])	/ (count + 1.0);
				//contAvg[n_evol][j] = (config[j] + 1 * contAvg[n_evol][j])	;
			}
		}


    	if(floor(100.0*count/N_of_samples) < floor(100.0*(count+1)/N_of_samples)){
    		cout << "Done: " << (count+1.0) / N_of_samples << "%" << endl;
    	}
    }

    auto end = chrono::steady_clock::now();

    // writting to a file
    int n_evol = 0;
    for (auto &vec : contAvg){
    	int i = -vec.size()/2;
		if (true || n_evol % 10 == 0 ) {
			sz_strm << "\"t=" << dt * n_evol << "\"\n\n" ;
			for (auto elem : vec) {
				sz_strm << i++ << "\t" << elem << "\t" << n_evol << "\n";
				//myFile << elem << "\t\t\t";
			}
			sz_strm << "\n\n";
		}
		n_evol++;
    }
    sz_strm << endl;
    //myFile << endl;

//    for(auto &vec: contDev){
//        for(auto elem : vec){
//        	myFile << elem << "\t\t\t";
//        }
//        myFile << endl;
//    }

    myFile.close();
    sz_strm. close();


    cout << "Computation time: " << "\t" << chrono::duration_cast<chrono::seconds>(end - start).count() << endl;



    return 0;
}








