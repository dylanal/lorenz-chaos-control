// include files( only used files from standard library
#include "stdafx.h"
#include <string>
#include <vector>
#include <cstdlib>
#include <iostream>
#include <cstdio>
#include <ctime>
#include <cmath>
#include <random>
#include <chrono>
#include <map>
#include <limits>
#include <sstream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <algorithm>
using namespace std;

// node class: includes value, children pointers, layer in tree and class constructors.
class node {
public:
	string value1;
	int value2;
	string value3;
	vector<node*> children;
	int layer;

	// initialisation constructors
	node(string data) {
		this->value1 = data;
		children.resize(2);
		for (int i = 0; i < children.size(); i++) {
			children[i] = new node();
		}
		if (data == "+" || data == "-" || data == "*" || data == "/") {
			this->children[0]->layer = this->layer++;
			this->children[1]->layer = this->layer++;
		}
		else if (data == "tan") {
			this->children[0]->layer = this->layer++;
		}
	}
	node() : children() {
		value1 = "";
		value3 = "";
		value2 = 0;
		layer = 1;
	}

	// assignment operator (copying one node into the other)
	node& operator=(node other)
	{

		value1 = other.value1;
		value2 = other.value2;
		value3 = other.value3;
		layer = other.layer;
		if (other.value1 == "tan") {
			children.resize(1);
			for (int i = 0; i < 1; i++) {
				this->children[i] = new node();
			}
		}
		else if (other.value1 == "+" || other.value1 == "-" || other.value1 == "*" || other.value1 == "/") {
			children.resize(2);
			for (int i = 0; i < 2; i++) {
				this->children[i] = new node();
			}
		}
		return *this;
	}
};

// tree class: includes layer/node map, cost value, genetic operation probability, root node pointer, 
class tree
{
public:
	map<int, vector<node *>> layers;
	float cost;
	double upfit;
	double downfit;
	int maxlayer;
	node* root;

	// member function to recursively copy value of node pointers of another tree without copying pointer itself
	void copy(node * Node1, node * Node2) {
		*Node1 = *Node2;


		if (Node2->value1 == "+" || Node2->value1 == "-" || Node2->value1 == "*" || Node2->value1 == "/") {
			copy(Node1->children[0], Node2->children[0]);
			copy(Node1->children[1], Node2->children[1]);
		}
		else if (Node2->value1 == "tan") {
			copy(Node1->children[0], Node2->children[0]);
		}
	}

	//tree constructor
	tree() {
		vector<node *> nodes;
		layers.insert(pair<int, vector<node *>>(1, nodes));
		vector<node *> nodes1;
		layers.insert(pair<int, vector<node *>>(2, nodes1));
		vector<node *> nodes2;
		layers.insert(pair<int, vector<node *>>(3, nodes2));
		vector<node *> nodes3;
		layers.insert(pair<int, vector<node *>>(4, nodes3));
		vector<node *> nodes4;
		layers.insert(pair<int, vector<node *>>(5, nodes4));
		cost = 0;
		upfit = 0;
		downfit = 0;
		maxlayer = 0;
		root = new node();
	}

	tree(string key, int max) {
		this->cost = 0;
		this->upfit = 0;
		this->downfit = 0;
		this->root = new node(key);
		this->root->layer = 1;
		this->maxlayer = max;
		vector<node *> nodes;
		layers.insert(pair<int, vector<node *>>(1, nodes));
		vector<node *> nodes1;
		layers.insert(pair<int, vector<node *>>(2, nodes1));
		vector<node *> nodes2;
		layers.insert(pair<int, vector<node *>>(3, nodes2));
		vector<node *> nodes3;
		layers.insert(pair<int, vector<node *>>(4, nodes3));
		vector<node *> nodes4;
		layers.insert(pair<int, vector<node *>>(5, nodes4));
		this->layers[1].push_back(root);
	}

	// tree copy constructor
	tree(tree &other)
	{
		vector<node *> nodes;
		layers.insert(pair<int, vector<node *>>(1, nodes));
		vector<node *> nodes1;
		layers.insert(pair<int, vector<node *>>(2, nodes1));
		vector<node *> nodes2;
		layers.insert(pair<int, vector<node *>>(3, nodes2));
		vector<node *> nodes3;
		layers.insert(pair<int, vector<node *>>(4, nodes3));
		vector<node *> nodes4;
		layers.insert(pair<int, vector<node *>>(5, nodes4));
		cost = other.cost;
		upfit = other.upfit;
		downfit = other.downfit;
		maxlayer = other.maxlayer;
		root = new node();
		this->copy(root, other.root);
	}
};

// function to add control law to lorenz system (used in lyapunov calculator part
double addtree(node* Node, double x1, double x2, double x3) {
	double output;
	if (Node->value1 == "" && Node->value3 == "") {
		output = Node->value2;
	}
	else if (Node->value1 == "" && Node->value2 == 0) {
		if (Node->value3 == "x1") {
			output = x1;
		}
		else if (Node->value3 == "x2") {
			output = x2;
		}
		else if (Node->value3 == "x3") {
			output = x3;
		}
	}
	else if (Node->value1 == "+") {
		output = addtree(Node->children[0], x1, x2, x2) + addtree(Node->children[1], x1, x2, x2);
	}
	else if (Node->value1 == "-") {
		output = addtree(Node->children[0], x1, x2, x2) - addtree(Node->children[1], x1, x2, x2);
	}
	else if (Node->value1 == "*") {
		output = addtree(Node->children[0], x1, x2, x2) * addtree(Node->children[1], x1, x2, x2);
	}
	else if (Node->value1 == "/") {
		output = addtree(Node->children[0], x1, x2, x2) / addtree(Node->children[1], x1, x2, x2);
	}
	else if (Node->value1 == "tan") {
		output = tan(addtree(Node->children[0], x1, x2, x2));
	}
	return output;
}


// ---------------------------------------------------------------------------------//
// parts of lyapunov calculator
const int ATTRACTOR = 0; // 0 - Lorenz 

double sig, rho, bet;

double N0 = 100; // Time after which Lyapunov exponent is evaluated 
double T = 100;  // Total simulation time 
double dT = 0.01;// Time Step (important for stability of runge kutta 4th order)

				 //x^2
double SQ(double x)
{
	return x*x;
}

// state equation 1
double fun_x1(double x1, double x2, double x3)
{
	
		return sig*(x2 - x1);
	

}

// state equation 2
double fun_x2(double x1, double x2, double x3)
{
	
		return rho*x1 - x1*x3 - x2;
	

}

// state equation 3 with control law added to it 
double fun_x3(double x1, double x2, double x3, node *root)
{
	
		return x1*x2 - bet*x3 + x3*x2 + addtree(root, x1, x2, x2);
	

}

// runge kutta 4th order
double rk4(double(*f)(double, double, double), double dt, double x1, double x2, double x3, int i)
{
	double	k1, k2, k3, k4;

	k1 = dt * f(x1, x2, x3);
	k2 = dt * f(x1 + 0.5*k1, x2 + 0.5*k1, x3 + 0.5*k1);
	k3 = dt * f(x1 + 0.5*k2, x2 + 0.5*k2, x3 + 0.5*k2);
	k4 = dt * f(x1 + k3, x2 + k3, x3 + k3);

	switch (i)
	{
	case 0:
		return x1 + (k1 + 2 * k2 + 2 * k3 + k4) / 6.0;
		break;
	case 1:
		return x2 + (k1 + 2 * k2 + 2 * k3 + k4) / 6.0;
		break;
	case 2:
		return x3 + (k1 + 2 * k2 + 2 * k3 + k4) / 6.0;
		break;
	}
}

// runge kutta 4th order for control law state equation
double rk41(double(*f)(double, double, double, node *), double dt, double x1, double x2, double x3, int i, node *root)
{
	double	k1, k2, k3, k4;

	k1 = dt * f(x1, x2, x3, root);
	k2 = dt * f(x1 + 0.5*k1, x2 + 0.5*k1, x3 + 0.5*k1, root);
	k3 = dt * f(x1 + 0.5*k2, x2 + 0.5*k2, x3 + 0.5*k2, root);
	k4 = dt * f(x1 + k3, x2 + k3, x3 + k3, root);


	switch (i)
	{
	case 0:

		return x1 + (k1 + 2 * k2 + 2 * k3 + k4) / 6.0;
		break;
	case 1:

		return x2 + (k1 + 2 * k2 + 2 * k3 + k4) / 6.0;
		break;
	case 2:

		return x3 + (k1 + 2 * k2 + 2 * k3 + k4) / 6.0;
		break;
	}
}


// utility classes for lyapunov calculator
const int dimension = 3;

struct nVector
{
	double r[dimension];
};

typedef vector <nVector> VnV;

class Utility
{
public:

	nVector AddVector(nVector, nVector);
	nVector SubVector(nVector, nVector);
	nVector ScalarMul(nVector, double);
	nVector NormalizeVector(nVector);
	double ScalarProduct(nVector, nVector);
	double ModVector(nVector);
	VnV DoOrthogonalization(VnV);
	VnV MatrixMultiply(VnV, VnV);

	VnV Transpose(VnV x)
	{
		VnV xT;

		xT = x;

		for (int i = 0; i< dimension; i++)
		{
			for (int j = 0; j< dimension; j++)
			{
				xT[j].r[i] = x[i].r[j];
			}
		}

		return xT;
	}

	void PrintVector(VnV x)
	{
		for (int i = 0; i< dimension; i++)
		{
			for (int j = 0; j< dimension; j++)
			{
				cout << x[i].r[j] << "    ";
			}
			cout << endl;
		}
	}

};




// lyapunov exponent calculator: outputs three exponents
vector<float> lyexp(node *root) {

	vector<float> out(3);
	

		sig = 10.0;
		rho = 28.0;
		bet = 8.0 / 3.0;
	



	Utility u;

	nVector dy;
	nVector dummy; // dynamical variable
	nVector lyapunov; // lyapunov exponent
	nVector sumly;
	VnV v;
	VnV v1;
	VnV unit; // for lyapunov exponent
	VnV J; // jacabobian

	for (int i = 0; i < dimension; i++) {
		dummy.r[i] = 0;
	}

	v.assign(dimension, dummy);
	v1.assign(dimension, dummy);
	J.assign(dimension, dummy);

	double x1;
	double x2;
	double x3;

	//initial condition

	dy.r[0] = 1;
	dy.r[1] = 1;
	dy.r[2] = 1;
	sumly.r[0] = 0;
	sumly.r[1] = 0;
	sumly.r[2] = 0;

	// initial orthogonal vector
	v[0].r[0] = 1;
	v[0].r[1] = 0;
	v[0].r[2] = 0;
	v[1].r[0] = 0;
	v[1].r[1] = 1;
	v[1].r[2] = 0;
	v[2].r[0] = 0;
	v[2].r[1] = 0;
	v[2].r[2] = 1;

	unit = v;

	for (int n = 1; n <= int(T / dT); n++)
	{

		x1 = rk4(fun_x1, dT, dy.r[0], dy.r[1], dy.r[2], 0);
		x2 = rk4(fun_x2, dT, dy.r[0], dy.r[1], dy.r[2], 1);
		x3 = rk41(fun_x3, dT, dy.r[0], dy.r[1], dy.r[2], 2, root);

		if (isnan(x1) || isnan(x2) || isnan(x3) || isinf(x1) || isinf(x2) || isinf(x3)) {

			out[0] = numeric_limits<float>::quiet_NaN();
			out[1] = numeric_limits<float>::quiet_NaN();
			out[2] = numeric_limits<float>::quiet_NaN();

			return out;
			goto end;
		}

		dy.r[0] = x1;
		dy.r[1] = x2;
		dy.r[2] = x3;

		if (n > int(N0))
		{

			if (ATTRACTOR == 0)
			{
				J[0].r[0] = -sig;
				J[0].r[1] = sig;
				J[0].r[2] = 0;
				J[1].r[0] = rho - x3;
				J[1].r[1] = -1;
				J[1].r[2] = -x1;
				J[2].r[0] = x2;
				J[2].r[1] = x1;
				J[2].r[2] = -bet;

			}
			else
			{
				J[0].r[0] = 0; J[0].r[1] = -1; J[0].r[2] = -1;
				J[1].r[0] = 1; J[1].r[1] = sig; J[1].r[2] = 0;
				J[2].r[0] = x3; J[2].r[1] = 0; J[2].r[2] = x1 - bet;
			}


			J[0] = u.AddVector(unit[0], u.ScalarMul(J[0], dT));
			J[1] = u.AddVector(unit[1], u.ScalarMul(J[1], dT));
			J[2] = u.AddVector(unit[2], u.ScalarMul(J[2], dT));

			v1 = u.DoOrthogonalization(u.MatrixMultiply(J, u.Transpose(v)));

			for (int k = 0; k< dimension; k++)
			{
				lyapunov.r[k] = log(u.ModVector(v1[k])) / dT;
				v[k] = u.NormalizeVector(v1[k]);
				sumly.r[k] += lyapunov.r[k];


			}
		}
	}

	out[0] = sumly.r[0] / (int(T / dT) - N0);
	out[1] = sumly.r[1] / (int(T / dT) - N0);
	out[2] = sumly.r[2] / (int(T / dT) - N0);
	return out;
end:;
}

//----------------------------------------------------------------------------------//



// function to interpret tree as string to output (needs to be given root node input)
string printtree(node* Node) {
	string output;
	if (Node->value1 == "" && Node->value3 == "") {
		stringstream ss;
		ss << Node->value2;
		string str = ss.str();
		output = str;
	}
	else if (Node->value1 == "" && Node->value2 == 0) {
		output = Node->value3;
	}
	else if (Node->value1 == "+" || Node->value1 == "-" || Node->value1 == "*" || Node->value1 == "/") {
		output = "(" + printtree(Node->children[0]) + ")" + " " + Node->value1 + " " + "(" + printtree(Node->children[1]) + ")";
	}
	else if (Node->value1 == "tan") {
		output = Node->value1 + "(" + printtree(Node->children[0]) + ")";
	}
	return output;
}

// recursive function to generate tree (if given root) or subtree (if given another node) of a certain depth
void recursion(tree * Tree, node* Node, int maxlayer, map<int, string> symbol, int symno, int constmin, int constmax) {
	// first part  for operations with two children
	if (Node->value1 == "+" || Node->value1 == "-" || Node->value1 == "*" || Node->value1 == "/") {
		// making sure tree ends at certain depth
		if (Node->layer == maxlayer - 1) {
			unsigned seed = chrono::system_clock::now().time_since_epoch().count();
			default_random_engine generator(seed);
			uniform_int_distribution<int> distribution(0, 1);
			int dis = distribution(generator);
			if (dis == 0) {
				unsigned seed = chrono::system_clock::now().time_since_epoch().count();
				default_random_engine generator(seed);
				uniform_int_distribution<int> distribution(constmin, constmax);
				int a = distribution(generator);
				int b = distribution(generator);
				Node->children.resize(2);
				for (int i = 0; i < Node->children.size(); i++) {
					Node->children[i] = new node();
				}
				Node->children[0]->value2 = a;
				Node->children[1]->value2 = b;
				Node->children[0]->value3 = "";
				Node->children[1]->value3 = "";
				Node->children[0]->value1 = "";
				Node->children[1]->value1 = "";
				Node->children[0]->layer = Node->layer + 1;
				Node->children[1]->layer = Node->layer + 1;
				Tree->layers[Node->children[0]->layer].push_back(Node->children[0]);
				Tree->layers[Node->children[1]->layer].push_back(Node->children[1]);
			}
			else if (dis == 1) {
				unsigned seed = chrono::system_clock::now().time_since_epoch().count();
				default_random_engine generator(seed);
				uniform_int_distribution<int> distribution(1, 3);
				int a = distribution(generator);
				int b = distribution(generator);
				Node->children.resize(2);
				for (int i = 0; i < Node->children.size(); i++) {
					Node->children[i] = new node();
				}
				switch (a) {
				case 1:
					Node->children[0]->value3 = "x1";
					break;
				case 2:
					Node->children[0]->value3 = "x2";
					break;
				case 3:
					Node->children[0]->value3 = "x3";
					break;
				}
				switch (b) {
				case 1:
					Node->children[1]->value3 = "x1";
					break;
				case 2:
					Node->children[1]->value3 = "x2";
					break;
				case 3:
					Node->children[1]->value3 = "x3";
					break;
				}
				Node->children[0]->value2 = 0;
				Node->children[1]->value2 = 0;
				Node->children[0]->value1 = "";
				Node->children[1]->value1 = "";
				Node->children[0]->layer = Node->layer + 1;
				Node->children[1]->layer = Node->layer + 1;
				Tree->layers[Node->children[0]->layer].push_back(Node->children[0]);
				Tree->layers[Node->children[1]->layer].push_back(Node->children[1]);

			}
		}
		// if not final layer add operation to next part of tree
		else {
			unsigned seed = chrono::system_clock::now().time_since_epoch().count();
			default_random_engine generator(seed);
			uniform_int_distribution<int> distribution(1, symno);
			int a = distribution(generator);
			int b = distribution(generator);
			Node->children.resize(2);
			for (int i = 0; i < Node->children.size(); i++) {
				Node->children[i] = new node();
			}
			Node->children[0]->value1 = symbol[a];
			Node->children[1]->value1 = symbol[b];
			Node->children[0]->value3 = "";
			Node->children[1]->value3 = "";
			Node->children[0]->value2 = 0;
			Node->children[1]->value2 = 0;
			Node->children[0]->layer = Node->layer + 1;
			Node->children[1]->layer = Node->layer + 1;
			Tree->layers[Node->children[0]->layer].push_back(Node->children[0]);
			Tree->layers[Node->children[1]->layer].push_back(Node->children[1]);
			for (int i = 0; i < Node->children.size(); i++) {
				recursion(Tree, Node->children[i], maxlayer, symbol, symno, constmin, constmax);
			}
		}
	}
	// second part for operations with one child
	if (Node->value1 == "tan") {
		if (Node->layer == maxlayer - 1) {
			unsigned seed = chrono::system_clock::now().time_since_epoch().count();
			default_random_engine generator(seed);
			uniform_int_distribution<int> distribution(0, 1);
			int dis = distribution(generator);
			if (dis == 0) {
				unsigned seed = chrono::system_clock::now().time_since_epoch().count();
				default_random_engine generator(seed);
				uniform_int_distribution<int> distribution(constmin, constmax);
				int a = distribution(generator);
				Node->children.resize(1);
				for (int i = 0; i < Node->children.size(); i++) {
					Node->children[i] = new node();
				}
				Node->children[0]->value3 = "";
				Node->children[0]->value2 = a;
				Node->children[0]->value1 = "";
				Node->children[0]->layer = Node->layer + 1;
				Tree->layers[Node->children[0]->layer].push_back(Node->children[0]);
			}
			else if (dis == 1) {
				unsigned seed = chrono::system_clock::now().time_since_epoch().count();
				default_random_engine generator(seed);
				uniform_int_distribution<int> distribution(1, 3);
				int a = distribution(generator);
				Node->children.resize(1);
				for (int i = 0; i < Node->children.size(); i++) {
					Node->children[i] = new node();
				}
				switch (a) {
				case 1:
					Node->children[0]->value3 = "x1";
					break;
				case 2:
					Node->children[0]->value3 = "x2";
					break;
				case 3:
					Node->children[0]->value3 = "x3";
					break;
				}
				Node->children[0]->value2 = 0;
				Node->children[0]->value1 = "";
				Node->children[0]->layer = Node->layer + 1;
				Tree->layers[Node->children[0]->layer].push_back(Node->children[0]);
			}
		}
		else {
			unsigned seed = chrono::system_clock::now().time_since_epoch().count();
			default_random_engine generator(seed);
			uniform_int_distribution<int> distribution(1, symno);
			int a = distribution(generator);
			Node->children.resize(1);
			for (int i = 0; i < Node->children.size(); i++) {
				Node->children[i] = new node();
			}
			Node->children[0]->value1 = symbol[a];
			Node->children[0]->value3 = "";
			Node->children[0]->value2 = 0;
			Node->children[0]->layer = Node->layer + 1;
			Tree->layers[Node->children[0]->layer].push_back(Node->children[0]);
			for (int i = 0; i < Node->children.size(); i++) {
				recursion(Tree, Node->children[i], maxlayer, symbol, symno, constmin, constmax);
			}
		}
	}
}

// function to random node in a given layer (used in genetic operation crossover) utilising tree class member variable layers map
node * subtree(tree *Tree, int layer) {

	int x = Tree->layers[layer].size();
	unsigned seed = chrono::system_clock::now().time_since_epoch().count();
	default_random_engine generator(seed);
	uniform_int_distribution<int> distribution(0, x - 1);
	int a = distribution(generator);
	if (Tree->layers[layer][a]->value1 == "") {
		return subtree(Tree, layer);
	}
	else {
		return Tree->layers[layer][a];
	}
}

// function to calculate probability of tree getting chosen to proceed to next generation (utilises softmax function)
void fitness(vector<tree *> treelaws) {
	float sum = 0;
	for (int i = 0; i < treelaws.size(); i++) {
		sum = sum + (-1*(treelaws[i]->cost));
	}
	treelaws[0]->downfit = 0;
	treelaws[0]->upfit = (-1*(treelaws[0]->cost)) / sum;
	for (int i = 1; i < treelaws.size(); i++) {
		treelaws[i]->downfit = treelaws[i - 1]->upfit;
		treelaws[i]->upfit = treelaws[i]->downfit + (-1*exp(treelaws[i]->cost)) / sum;
	}
}

// function to update tree class member variable layers after genetic operation performed
void updatelayer(tree *Tree, node *Node) {
	if (Node->layer == Tree->maxlayer) {
		Tree->layers[Node->layer].push_back(Node);
	}
	else {
		Tree->layers[Node->layer].push_back(Node);
		for (int i = 0; i < Node->children.size(); i++) {
			updatelayer(Tree, Node->children[i]);
		}
	}
}

// function which performs all genetic operations given input of trees and probabilities for each genetic operation to occur
vector<tree *> geneop(int gennum, double percent, vector<tree *> treelaws, double probrec, double probrep, double probmut, double probcross, map<int, string> symbol) {
	vector<tree *> trees;
	int y = treelaws.size();
	int x = ceil((percent / 100)*y);
	for (int i = 0; i < x; i++) {
		tree *Tree1p = new tree();
		*Tree1p = *treelaws[i];
		updatelayer(Tree1p, Tree1p->root);
		trees.push_back(Tree1p);
		
	}
	
	int sum = x;
	while (sum < gennum) {
		unsigned seed = chrono::system_clock::now().time_since_epoch().count();
		default_random_engine generator(seed);
		uniform_real_distribution<double> distribution(0.0, 1.0);
		double l = distribution(generator);
		if (l <= probrep) {
			unsigned seed = chrono::system_clock::now().time_since_epoch().count();
			default_random_engine generator(seed);
			uniform_real_distribution<double> distribution(0.0, 1.0);
			double a = distribution(generator);
			for (int i = max(x - 1, 0); i < y; i++) {
				if ((a >= treelaws[i]->downfit) && (a < treelaws[i]->upfit)) {


					tree *Tree1p = new tree();
					*Tree1p = *treelaws[i];
					updatelayer(Tree1p, Tree1p->root);
					trees.push_back(Tree1p);
					sum++;
					goto stop;
				}
			}
		}
		else if ((l > probrep) && (l <= (probrep + probrec))) {
			unsigned seed = chrono::system_clock::now().time_since_epoch().count();
			default_random_engine generator(seed);
			uniform_int_distribution<int> distribution(1, 5);
			int gh = distribution(generator);
			string key = symbol[gh];
			tree *treelaw = new tree(key, 5);
			recursion(treelaw, treelaw->root, 5, symbol, 5, -10, 10);
			trees.push_back(treelaw);
		}
		else if ((l > (probrep + probrec)) && (l <= (probrep + probrec + probcross))) {

			unsigned seed = chrono::system_clock::now().time_since_epoch().count();
			default_random_engine generator(seed);
			uniform_real_distribution<double> distribution(0.0, 1.0);
			double a = distribution(generator);
			for (int i = 0; i < y; i++) {


				if ((a >= treelaws[i]->downfit) && (a < treelaws[i]->upfit)) {

					unsigned seed = chrono::system_clock::now().time_since_epoch().count();
					default_random_engine generator(seed);
					uniform_real_distribution<double> distribution(0.0, 1.0);
					double b = distribution(generator);

					for (int j = 0; j < y; j++) {


						if ((b >= treelaws[j]->downfit) && (b < treelaws[j]->upfit)) {

							tree *Tree1p = new tree();
							*Tree1p = *treelaws[i];

							tree *Tree2p = new tree();
							*Tree2p = *treelaws[j];

							updatelayer(Tree1p, Tree1p->root);
							updatelayer(Tree2p, Tree2p->root);
							node *one = subtree(Tree1p, 2);
							node *two = subtree(Tree2p, 2);

							if (one->value1 == "tan") {

								node *three = one->children[0];

								if (two->value1 == "tan") {

									node *four = two->children[0];
									one->children[0] = four;
									two->children[0] = three;
									trees.push_back(Tree1p);
									trees.push_back(Tree2p);
									sum = sum + 2;

									goto stop;
								}
								else {
									unsigned seed = chrono::system_clock::now().time_since_epoch().count();
									default_random_engine generator(seed);
									uniform_int_distribution<int> distribution(0, 1);
									int c = distribution(generator);
									node *four = two->children[c];
									one->children[0] = four;
									two->children[c] = three;
									trees.push_back(Tree1p);
									trees.push_back(Tree2p);
									
									sum = sum + 2;

									goto stop;
								}
							}
							else {
								unsigned seed = chrono::system_clock::now().time_since_epoch().count();
								default_random_engine generator(seed);
								uniform_int_distribution<int> distribution(0, 1);
								int d = distribution(generator);
								node *three = one->children[d];

								if (two->value1 == "tan") {
									node *four = two->children[0];
									one->children[d] = four;
									two->children[0] = three;
									trees.push_back(Tree1p);
									trees.push_back(Tree2p);
									
									sum = sum + 2;

									goto stop;
								}
								else {
									
									unsigned seed = chrono::system_clock::now().time_since_epoch().count();
									default_random_engine generator(seed);
									uniform_int_distribution<int> distribution(0, 1);
									int e = distribution(generator);
									node *four = two->children[e];
									one->children[d] = four;
									two->children[e] = three;
									trees.push_back(Tree1p);
									trees.push_back(Tree2p);
									
									sum = sum + 2;
									goto stop;
								}
							}
						}
					}
				}
			}
		}
		else if (l > (probrep + probrec + probcross)) {

			unsigned seed = chrono::system_clock::now().time_since_epoch().count();
			default_random_engine generator(seed);
			uniform_real_distribution<double> distribution(0.0, 1.0);
			double a = distribution(generator);

			for (int i = 0; i < y; i++) {

				if ((a >= treelaws[i]->downfit) && (a < treelaws[i]->upfit)) {

					tree *Tree1p = new tree();
					*Tree1p = *treelaws[i];
					updatelayer(Tree1p, Tree1p->root);
					node *one = subtree(Tree1p, 2);
					recursion(Tree1p, one, 5, symbol, 5, -10, 10);
					trees.push_back(Tree1p);
					sum++;
					goto stop;
				}
			}
		}
	stop:;
	}

	treelaws.clear();
	return trees;
}


//function to generate first generation of population, symbols = '-,+,*,/,tan' 
vector<string> control_gen(int maxlayer, map<int, string> symbol, int symno, int constmin, int constmax, int treeno, vector<tree *> treelaws) {
	vector<string> controls;
	for (int i = 0; i < treeno; i++) {

		recursion(treelaws[i], treelaws[i]->root, maxlayer, symbol, symno, constmin, constmax);

		controls.push_back(printtree(treelaws[i]->root));
	}
	return controls;
}

// function by which to sort the costs of the trees
bool myfunction(tree *i, tree *j) { return (i->cost < j->cost); }

// function to delete tree class member variable layers before updating it
void deletelayer(tree *Tree, int maxlayer) {
	for (int i = 1; i < maxlayer + 1; i++) {
		Tree->layers[i].clear();
	}
}

//----------------------------------------------------------------------------------//


int main()
{

	// parameters of script to change
	int treeno = 500; // initial population size
	int genum = 200; // population size in following generations
	int generations = 5; // number of iterations to perform

	// map linking operation to number (for use in recursion to randomly select operation)
	map<int, string> oper;
	oper.insert(pair<int, string>(1, "+"));
	oper.insert(pair<int, string>(2, "-"));
	oper.insert(pair<int, string>(3, "*"));
	oper.insert(pair<int, string>(4, "/"));
	oper.insert(pair<int, string>(5, "tan"));

	// initialise initial tree population
	vector<tree *> treelaws;
	for (int i = 0; i < treeno; i++) {
		unsigned seed = chrono::system_clock::now().time_since_epoch().count();
		default_random_engine generator(seed);
		uniform_int_distribution<int> distribution(1, 5);
		int a = distribution(generator);
		string key = oper[a];
		tree *treelaw = new tree(key, 5);
		treelaws.push_back(treelaw);
	}
	vector<string> control_laws = control_gen(5, oper, 5, -10, 10, treeno, treelaws);
	cout << "gen1" << endl;
	cout << endl;


	// initial sorting and printing out of best tree
	for (int i = 0; i < treelaws.size(); i++) {
		vector<float> out = {};
		out.resize(3);
		out = lyexp(treelaws[i]->root);
		// deleting error prone control laws
		if (isnan(out[0]) || isnan(out[1]) || isnan(out[2]) || isinf(out[0]) || isinf(out[1]) || isinf(out[2])) {
			treelaws.erase(treelaws.begin() + i);
			i--;
		}
		else {
			treelaws[i]->cost = exp(*max_element(out.begin(), out.end()));
		}
	}
	sort(treelaws.begin(), treelaws.end(), myfunction);
	cout << "cost 1: " << treelaws[0]->cost << " " << log(treelaws[0]->cost) << " " << treelaws[1]->cost << endl;
	cout << printtree(treelaws[0]->root) << endl;
	fitness(treelaws);
	


	// for however many generation perform genetic operations on population and print out the best tree for each generation
	for (int gen = 0; gen < generations; gen++) {
		cout << endl;
		cout << "gen " << gen + 2 << endl;
		cout << endl;
		treelaws = geneop(genum, 10, treelaws, 0.4, 0, 0.4, 0.2, oper);
		

		for (int i = 0; i < treelaws.size(); i++) {
			
			deletelayer(treelaws[i], 5);

			updatelayer(treelaws[i], treelaws[i]->root);

			vector<float> out = {};
			out.resize(3);
			out = lyexp(treelaws[i]->root);
			// deleting error prone control laws
			if (isnan(out[0]) || isnan(out[1]) || isnan(out[2]) || isinf(out[0]) || isinf(out[1]) || isinf(out[2])) {
				treelaws.erase(treelaws.begin() + i);
				i--;
			}
			else {
				treelaws[i]->cost = exp(*max_element(out.begin(), out.end()));
			}

		}
		sort(treelaws.begin(), treelaws.end(), myfunction);
		cout << "cost " << gen + 2 << ": " << treelaws[0]->cost << " " << log(treelaws[0]->cost) << endl;
		cout << printtree(treelaws[0]->root) << endl;
		fitness(treelaws);
	}
	cout << endl;


	// print out best control law overall
	cout << "The best control law is: " << printtree(treelaws[0]->root) << "with best cost : " << treelaws[0]->cost << endl;

	return 0;
}

//----------------------------------------------------------------------------------//
// code for utility class declarations above

VnV Utility::MatrixMultiply(VnV A, VnV B)
{
	nVector dummy;

	for (int i = 0; i < dimension; i++) {
		dummy.r[i] = 0;
	}

	VnV C;
	C.assign(dimension, dummy);

	for (int i = 0; i< dimension; i++)
	{
		for (int j = 0; j< dimension; j++)
		{
			for (int k = 0; k< dimension; k++)
			{
				C[i].r[j] += A[i].r[k] * B[k].r[j];
			}
		}
	}

	return Transpose(C);
}

double Utility::ModVector(nVector v1)
{
	double sum = 0;

	for (int i = 0; i< dimension; i++)
	{
		sum += v1.r[i] * v1.r[i];
	}

	return sqrt(sum);
}

double Utility::ScalarProduct(nVector v1, nVector v2)
{
	double sum = 0;

	for (int i = 0; i< dimension; i++)
	{
		sum += v1.r[i] * v2.r[i];
	}
	return sum;
}

nVector Utility::NormalizeVector(nVector v1)
{
	double sum = 0;
	double mul;

	nVector result;

	for (int i = 0; i< dimension; i++)
	{
		sum += v1.r[i] * v1.r[i];
	}

	for (int i = 0; i< dimension; i++)
	{
		result.r[i] = v1.r[i] / sqrt(sum);
	}

	return result;
}

nVector Utility::AddVector(nVector v1, nVector v2)
{
	nVector result;

	for (int i = 0; i< dimension; i++)
	{
		result.r[i] = v1.r[i] + v2.r[i];
	}

	return result;
}

nVector Utility::SubVector(nVector v1, nVector v2)
{
	nVector result;

	for (int i = 0; i< dimension; i++)
	{
		result.r[i] = v1.r[i] - v2.r[i];
	}

	return result;
}

nVector Utility::ScalarMul(nVector v1, double x)
{
	nVector result;

	for (int i = 0; i< dimension; i++)
	{
		result.r[i] = v1.r[i] * x;
	}

	return result;
}

VnV Utility::DoOrthogonalization(VnV in)
{
	VnV out, out_normalize;
	nVector sum, dummy;

	for (int i = 0; i < dimension; i++) {
		dummy.r[i] = 0;
	}

	out.assign(dimension, dummy);

	out[0] = in[0];

	for (int i = 1; i< in.size(); i++)
	{
		sum = dummy;
		for (int j = 0; j <= i - 1; j++)
		{
			sum = AddVector(sum, ScalarMul(out[j], ScalarProduct(out[j], in[i]) / SQ(ModVector(out[j]))));
		}

		out[i] = SubVector(in[i], sum);


	}

	return out;
}

//----------------------------------------------------------------------------------//


