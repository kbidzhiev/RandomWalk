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


void do_jump (vector<int>& config, int particle_position){
	size_t distance = 1;
	swap(config[particle_position],config[particle_position + distance]);
};

void evolveState_Classic(vector<int>& config, double prob_to_move){

	double prob = distReal(generator);
	int jump_distance = 1;
	int System_size = config.size();
	int particle_position (System_size/2);

	for (size_t i = 0; i < config.size(); ++i){
		if (config[i] == 1) particle_position = i;
	}

	if (particle_position > jump_distance
			&& particle_position < (System_size - jump_distance)) {
		if (prob < prob_to_move) {
			do_jump(config, -particle_position);
		} else if (prob > 1 - prob_to_move) {
			do_jump(config, particle_position);
		} else {
			// we don't move, just stay at the same position
		}
	} else if (particle_position <= jump_distance) {
		do_jump(config, particle_position);
	} else if (particle_position >= System_size - jump_distance) {
		do_jump(config, -particle_position);
	}
}
