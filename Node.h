#pragma once
#include <iostream>
#include <fstream>
#include <string>
using namespace std;
class Node {
	float x;
	float y;
public:
	Node();
	Node(float, float);
	void set(float, float);
	float get_x();
	float get_y();
	~Node();
};
