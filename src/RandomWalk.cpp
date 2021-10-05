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
#include <fstream>
#include <sstream>
#include <string>
#include <omp.h>

using namespace std;



const int L = 144; //4*250
const int T = sqrt(L); // 20000

const double dt =  0.1;

const double prob_to_move = pow(dt, 1);

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
	sz_strm.open("Data/Sz_dt" + dt_for_filename +".dat", mode);
	myFile.open("Data/classical_jammed.txt", mode);
	sz_strm.precision(9);
	myFile.precision(9);

	sz_strm << "# x \t sigma_z \t t" << endl;

    ////////////////////////////////////////////////////////////////

	const int Evolution_steps = T/dt;
    vector<vector<double> > contAvg(Evolution_steps, vector<double>(L,0)); // vector<T> vec (how many elements, initialize by {.} )
    vector<vector<double> > contDev(contAvg);


    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    mt19937_64 generator(seed);
    uniform_real_distribution<double> distReal(0,1);

//    double prob = 0.5; // = 1/6 //
    vector<int> config(L,0);
    vector<int> posLeftMovers;
    vector<int> posRightMovers;
    posLeftMovers.reserve(L);
    posRightMovers.reserve(L);
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    auto start = chrono::steady_clock::now();
//#pragma omp parallel for num_threads(omp_get_num_procs())
    for(int count=0; count<N_of_samples; count++){

    	//initial state
    	for (size_t j = 0; j < L ; j++) {
    		config[j] = 0;
//    		if (j % 2 == 0 ){
//    			config[j] = 0;
//    		} else {
//    			config[j] = 1;
//    		}
    	}

    	// local perturbation
    	config[L/2] = 1;
//    	if(config[L/2] == 1) {
//    		config[L/2]   = 0;
//    	} else {
//    		config[L/2+1] = 0;
//    	}

    	//config[L/2+1] = 0;

//    	config[L/2+2] = 0;
//    	config[L/2+3] = 0;



    	for (size_t n_evol = 0; n_evol < Evolution_steps; ++n_evol){
    		//cout << n_evol<< endl;
      	   // saving the positions of movable triples
			for (size_t j = 1; j < L; j++) {
				if (config[j - 1] == 0 && config[j] == 1) {
					posRightMovers.push_back(j);
				} else if (config[j - 1] == 1 && config[j] == 0) {
					posLeftMovers.push_back(j);
				}
			}


      	    int size = posLeftMovers.size();
      	    uniform_int_distribution<int> distPosition(0,posRightMovers.size()+size-1);
      	    //int pos = distPosition(generator);
			uniform_real_distribution<double> distReal(0,1);
			double prob = distReal(generator);

			if (posLeftMovers.size() != 0 && posRightMovers.size() != 0){
				if (prob < prob_to_move) {
					config[posLeftMovers[0] - 1] = 0;
					config[posLeftMovers[0] + 0] = 1;
				} else if (prob >= 1. - prob_to_move) {
					config[posRightMovers[0] - 1] = 1;
					config[posRightMovers[0] + 0] = 0;
				}
			} else if (posLeftMovers.size() == 0){
				if (prob < prob_to_move/2.0) {
					config[posRightMovers[0] - 1] = 1;
					config[posRightMovers[0] + 0] = 0;
				}
			} else if (posRightMovers.size() == 0){
				if (prob < prob_to_move/2.0 ) {
					config[posLeftMovers[0] - 1] = 0;
					config[posLeftMovers[0] + 0] = 1;
				}
			} else {
				cout << "No particles to move " << endl;
				exit(1);
			}


     	  	posLeftMovers.clear();
      	   	posRightMovers.clear();

		// we start collecting data

			for (int j = 0; j < L; j++) {
				contDev[n_evol][j] = (contDev[n_evol][j] + (1.0 * config[j] - contAvg[n_evol][j]) * (1.0 * config[j] - contAvg[n_evol][j])
								/ (count + 1)) * count / (count + 1);
				contAvg[n_evol][j] = (config[j] + count * contAvg[n_evol][j])	/ (count + 1);
				//contAvg[n_evol][j] = (config[j] + 1 * contAvg[n_evol][j])	;

			}
		}

	if(count%100==0) cout << "Done: " << count/1000 + 1 << " / " << N_of_samples/1000 << endl;
    }
    auto end = chrono::steady_clock::now();

    // writting to a file
    int n_evol = 0;
    for (auto &vec : contAvg){
    	int i = -vec.size()/2;
		if (true || n_evol % 11 == 0 ) {
			sz_strm << "\"t=" << n_evol*dt << "\"\n" << endl;
			for (auto elem : vec) {
				sz_strm << i++ << "\t" << elem << "\t" << n_evol*dt << "\n";
				myFile << elem << "\t\t\t";
			}
			sz_strm << "\n" << endl;
		}
		n_evol++;


    }
    myFile << endl;

    for(auto &vec: contDev){
        for(auto elem : vec){
        	myFile << elem << "\t\t\t";
        }
        myFile << endl;
    }


    myFile.close();
    sz_strm. close();

    cout << "Elapsed time: " << "\t" << chrono::duration_cast<chrono::seconds>(end - start).count() << endl;



    return 0;
}








