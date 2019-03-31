#include "Node.h"
Node::Node()
{

}
Node::Node(float X, float Y)
{
	x = X;
	y = Y;
}
void Node::set(float X, float Y)
{
	x = X;
	y = Y;
}
float Node::get_x()
{
	return x;
}
float Node::get_y()
{
	return y;
}
Node::~Node()
{

}