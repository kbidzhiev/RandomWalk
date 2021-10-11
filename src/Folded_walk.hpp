#pragma once

#include <vector>

using namespace std;

//random_device rd;  // Will be used to obtain a seed for the random number engine
//mt19937 generator(rd()); // Standard mersenne_twister_engine seeded with rd()
//std::uniform_real_distribution<> distReal(0.0, 1.0);


vector<int> initialState_Folded(const int flip_position, const size_t System_size){
	vector<int> config(System_size);
	config.reserve(System_size);

	for (size_t i = 0; i < System_size; ++i){
		if (i %3 == 0) 	config.push_back(0);
		else config.push_back(1);
	}
	config[flip_position] = 0;
	return config;
}
