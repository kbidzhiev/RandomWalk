#pragma once

#include <vector>
#include <string>

using namespace std;


enum TYPE{
	Classic,
	Folded
};

class State{
public:
	State (size_t System_size, size_t spin_flip_pos, TYPE type)
	: type(type)
	, config(vector<int> (System_size ,0)){
		 if (type == TYPE::Classic){
			 config[spin_flip_pos] = 1;

		 } else if(type == TYPE::Folded){
			 for (size_t i = 0; i < config.size(); ++i){

			 }
		 }
	}



private:
	vector<int> config;
	vector<size_t> left_movers;
	vector<size_t> right_movers;
	TYPE type;
};
