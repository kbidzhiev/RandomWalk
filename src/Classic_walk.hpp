#pragma once


#include <vector>



using namespace std;

random_device rd;  // Will be used to obtain a seed for the random number engine
mt19937 generator(rd()); // Standard mersenne_twister_engine seeded with rd()
std::uniform_real_distribution<> distReal(0.0, 1.0);



vector<int> initialState_Classic(const size_t flip_position, const size_t System_size){
	vector<int> config(System_size ,0);
	config[flip_position] = 1;
	return config;
}


void do_jump_Classic (vector<int>& config, size_t particle_position, int direction){
	swap(config[particle_position],config[particle_position + direction]);
};

void evolveState_Classic(vector<int>& config, double prob_to_move){

	double prob = distReal(generator);
	int jump_distance = 1;
	int System_size = config.size();
	int particle_position(System_size/2);

	int direction = 0;
	for (size_t i = 0; i < config.size(); ++i){
		if (config[i] == 1) particle_position = i;
	}


	if (particle_position > jump_distance
			&& particle_position < (System_size - jump_distance)) {
		if (prob < prob_to_move) {
			direction = - 1;
		} else if (prob > 1 - prob_to_move) {
			direction =  1;
		} else {
			// we don't move, just stay at the same position
		}
	} else if (particle_position <= jump_distance) {
		direction = + 1;
	} else if (particle_position >= System_size - jump_distance) {
		direction = - 1;
	}
	do_jump_Classic(config, particle_position, direction);
}
