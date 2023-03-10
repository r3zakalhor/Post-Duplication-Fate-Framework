// reconstruct gene tree.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <cinttypes>
#include <cstdlib>
#include <cerrno>
#include <sys/stat.h>
#include <list>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <string>
#include <cmath>
#include "probabilitiescalculation.h"

//#include <qgraphicsview>
//#include <qtextstream>
//#include <qprocess>
using namespace std;

// =================================================================
//                              Includes
// =================================================================
struct protein {
public:
	int nodeid, id, parent_id;
	long double m, h, w;
	int32_t start_pos_, length_, basal_level_, hamming_dist_;
};

struct Node {
public:
	protein pro;
	bool dup;
	vector<Node*> children;
	struct Node* left; //Reference to left child
	struct Node* right; //Reference to right child

	Node() { this->left = NULL; this->right = NULL; }

	bool IsLeaf() {
		//if (this->left == NULL && this->right == NULL)
		if (this->GetNbChildren() == 0)
			return true;
		else
			return false;
	}

	string GetLabel()
	{
		return to_string(pro.nodeid);
	}

	int GetNbChildren()
	{
		/*if (this->left != NULL && this->right != NULL)
			return 2;
		else
			return 1;*/
		return this->children.size();
	}

	Node* GetChild(int pos)
	{
		if (pos < 0 || pos >= this->GetNbChildren())
			throw "Bad child position";
		/*if (pos == 0)
			return this->left;
		else if (pos == 1)
			return this->right;*/
		else
			return this->children[pos];
	}
};

class GeneTree
{
public:
	//void init(QGraphicsScene* scene, QGraphicsView* view);
	void insert(int a);


	Node* construct_tree(fstream& proteins_file, long int gen) {
		vector<Node*> current_level, next_level;
		int numofnodes = 0;
		int time_gen = 0;
		Node* root = new Node();
		root->pro.nodeid = numofnodes;
		numofnodes++;
		root->dup = false;
		root->pro.id = 0;//stoi(pro_list[0][1]);
		root->pro.parent_id = 0;//stoi(pro_list[0][2]);
		/*root->pro.start_pos_ = stoi(pro_list[0][3]);
		root->pro.length_ = stoi(pro_list[0][4]);
		root->pro.m = stold(pro_list[0][10]);
		root->pro.w = stold(pro_list[0][11]);
		root->pro.h = stold(pro_list[0][12]);*/
		current_level.push_back(root);
		int numofpro = 0, i_numofpro = 0;
		//cout << current_level[0]->pro.id<< endl;
		//cout << initial->pro.id << endl;
		//current_level[0]->pro.id = 1;
		//cout << current_level[0]->pro.id << endl;
		//cout << initial->pro.id << endl;

		string line;
		getline(proteins_file, line);	//skip the header line
		while (getline(proteins_file, line) && time_gen <= gen)
		{
			stringstream str(line);

			vector<string> row;
			string word;
			while (getline(str, word, ','))
				row.push_back(word);

			if (row[0].size() > 1) {

				numofpro = stoi(row[8]);
				//cout << "num of pro: " << numofpro << endl;
				i_numofpro = 0;
				Node* next = new Node();
				next->pro.nodeid = numofnodes;
				numofnodes++;
				next->dup = false;
				next->pro.id = stoi(row[1]);
				next->pro.parent_id = stoi(row[2]);
				next->pro.start_pos_ = stoi(row[3]);
				next->pro.length_ = stoi(row[4]);
				next->pro.m = stold(row[10]);
				next->pro.w = stold(row[11]);
				next->pro.h = stold(row[12]);
				time_gen = stoi(row[15]);
				next_level.push_back(next);
				i_numofpro++;
				//cout << "i_numofpro: " << i_numofpro << endl;
			}
			else if (i_numofpro < numofpro) {

				//TODO: evil duplicated code
				Node* next = new Node;
				next->pro.nodeid = numofnodes;
				numofnodes++;
				next->dup = false;
				next->pro.id = stoi(row[1]);
				next->pro.parent_id = stoi(row[2]);
				next->pro.start_pos_ = stoi(row[3]);
				next->pro.length_ = stoi(row[4]);
				next->pro.m = stold(row[10]);
				next->pro.w = stold(row[11]);
				next->pro.h = stold(row[12]);
				time_gen = stoi(row[15]);
				next_level.push_back(next);
				i_numofpro++;
				//cout << "i_numofpro: " << i_numofpro << endl;
			}
			if (i_numofpro == numofpro) {

				for (int m = 0; m < next_level.size(); m++) {

					bool flag = false;
					for (int n = 0; n < current_level.size(); n++) {

						if (current_level[n]->pro.id == next_level[m]->pro.parent_id) {

							current_level[n]->children.push_back(next_level[m]);
							flag = true;

						}

						if (current_level[n]->GetNbChildren() > 1) {
							current_level[n]->dup = true;

						}
					}
					if (flag == false) {
						root->children.push_back(next_level[m]);
					}
				}
				current_level.clear();
				current_level = next_level;
				next_level.clear();
				i_numofpro = 0;
			}

		}
		cout << "last event at generation: " << gen << endl;
		root->dup = false;
		return root;
	}

	void preorderWalk();
	void postorderWalk();
	void inorderWalk();

	void deleteNode(int val);
	void deleteNode(Node* p);
	void show();

	int countNodes();
	int countLevels();
	int countLeftNodes();
	Node* findElem(int val);
protected:
private:
	int countNodes(Node* p);
	int countLevels(Node* p);
	int countLeftNodes(Node* p);

	void preorderWalk(Node* p);
	void postorderWalk(Node* p);
	void inorderWalk(Node* p);

	Node* findSuccessor(int val);

	//QByteArray _prepareGraph();
	//void _graphWalk(Node* p, QTextStream* stream);
	Node* findElem(int val, Node* p);

	Node* _root;
	//QGraphicsScene* _scene;
	//QGraphicsView* _view;
};

struct Trunk
{
	Trunk* prev;
	string str;

	Trunk(Trunk* prev, string str)
	{
		this->prev = prev;
		this->str = str;
	}
};

// Helper function to print branches of the binary tree
void showTrunks(Trunk* p)
{
	if (p == nullptr) {
		return;
	}

	showTrunks(p->prev);
	cout << p->str;
}

void printTree(Node* root, Trunk* prev, bool isLeft)
{
	if (root == nullptr) {
		return;
	}

	string prev_str = "    ";
	Trunk* trunk = new Trunk(prev, prev_str);

	printTree(root->right, trunk, true);

	if (!prev) {
		trunk->str = "---";
	}
	else if (isLeft)
	{
		trunk->str = ".---";
		prev_str = "   |";
	}
	else {
		trunk->str = "`---";
		prev->str = prev_str;
	}

	showTrunks(trunk);
	cout << " " << root->pro.nodeid << endl;

	if (prev) {
		prev->str = prev_str;
	}
	trunk->str = "   |";

	printTree(root->left, trunk, false);
}

void ReturnLeafNodes(Node* root, vector<protein>& leaflist)
{
	// if node is null, return
	if (!root)
		return;

	// if node is leaf node, print its data
	if (root->GetNbChildren() == 0)
	{
		//cout << root->data << " ";
		leaflist.push_back(root->pro);
		return;
	}

	// if left child exists, check for leaf
	// recursively
	if (root->GetNbChildren() == 1) {
		if (root->children[0])
			ReturnLeafNodes(root->children[0], leaflist);
	}

	// if right child exists, check for leaf
	// recursively
	if (root->GetNbChildren() == 2) {
		if (root->children[0])
			ReturnLeafNodes(root->children[0], leaflist);
		if (root->children[1])
			ReturnLeafNodes(root->children[1], leaflist);
	}
}

void WriteNodeChildren(string& str, Node* curNode, bool addBranchLengthToLabel, bool addInternalNodesLabel)
{
	if (curNode->IsLeaf())
	{

		str += curNode->GetLabel();

		//if (addBranchLengthToLabel && !curNode->IsRoot())
			//str += ":" + Util::ToString(curNode->GetBranchLength());
	}
	else if (curNode->GetNbChildren() > 1)
	{
		str += "(";
		for (int i = 0; i < curNode->GetNbChildren(); i++)
		{
			if (i != 0)
				str += ", ";

			Node* child = curNode->GetChild(i);
			WriteNodeChildren(str, child, addBranchLengthToLabel, addInternalNodesLabel);
		}

		str += ")";

		if (addInternalNodesLabel)
			str += curNode->GetLabel();

		//if (addBranchLengthToLabel && !curNode->IsRoot())
			//str += ":" + Util::ToString(curNode->GetBranchLength());
	}
	else if (curNode->GetNbChildren() == 1) {
		Node* child = curNode->GetChild(0);
		WriteNodeChildren(str, child, addBranchLengthToLabel, addInternalNodesLabel);
	}

}

string ToNewickString(Node* root, bool addBranchLengthToLabel, bool addInternalNodesLabel)
{
	string str;
	WriteNodeChildren(str, root, addBranchLengthToLabel, addInternalNodesLabel);
	str += ";";
	return str;
}


void ClassifyDuplicationFates(Node* root, ofstream& dups_list_file, ofstream& avg_dups_list_file) {
	vector<protein> left_leaflist;
	vector<protein> right_leaflist;
	if (!root)
		return;

	vector<long double> fates_probablities;
	long double sub_avg, cons_avg, new_avg, pseudo_avg, spec_avg, dblneo_avg;
	int numofleaves;

	if (root->dup == true) {
		if (root->GetNbChildren() == 1) {
			if (root->children[0] != nullptr) {
				ReturnLeafNodes(root->children[0], left_leaflist);
				if (left_leaflist.size() == 0)
					left_leaflist.push_back(root->children[0]->pro);
			}
		}
		if (root->GetNbChildren() == 2) {
			if (root->children[0] != nullptr) {
				ReturnLeafNodes(root->children[0], left_leaflist);
				if (left_leaflist.size() == 0)
					left_leaflist.push_back(root->children[0]->pro);
			}
			if (root->children[1] != nullptr) {
				ReturnLeafNodes(root->children[1], right_leaflist);
				if (right_leaflist.size() == 0)
					right_leaflist.push_back(root->children[1]->pro);
			}
		}
		long double g[3] = { root->pro.m, root->pro.h, root->pro.w };
		/*cout << endl << "node: " << root->pro.nodeid << endl;
		cout << "left_leaves: ";
		for (int i = 0; i < left_leaflist.size(); i++)
			cout << left_leaflist[i].nodeid << "  ";
		cout << "right_leaves: ";
		for (int i = 0; i < right_leaflist.size(); i++)
			cout << right_leaflist[i].nodeid << "  ";
		*/
		/*if (left_leaflist.size() == 1 && right_leaflist.size() == 1 && root->left == nullptr && root->right == nullptr) {
			// m, h, w
			long double a[3] = { leaflist[0].m, leaflist[0].h, leaflist[0].w };
			long double b[3] = { leaflist[1].m, leaflist[1].h, leaflist[1].w };
			cout << endl;
			fates_probablities = ProbabilitiesCalculation(g, a, b);
			dups_list_file << root->pro.nodeid << "," << g[0] << "," << g[1] << "," << g[2] << ","
				<< leaflist[0].nodeid << "," << a[0] << "," << a[1] << "," << a[2] << ","
				<< leaflist[1].nodeid << "," << b[0] << "," << b[1] << "," << b[2] << ","
				<< fates_probablities[0] << "," << fates_probablities[1] << "," << fates_probablities[2] << "," << fates_probablities[3] << "," << fates_probablities[4] << std::endl;
			//dups_list_file << endl;
			//dups_list_file << "21original node, g.m, g.h, g.w, first descendant id, a.m, a.h, a.w, second descendant id, b.m, b.h, b.w, P_subfunctionlization, P_conservation, P_newfunctionlizatoin, P_pseudogenization, P_specialization " << std::endl;

		}*/
		if (left_leaflist.size() > 0 && right_leaflist.size() > 0) {
			numofleaves = left_leaflist.size() * right_leaflist.size();
			sub_avg = 0; cons_avg = 0; new_avg = 0; pseudo_avg = 0; spec_avg = 0; dblneo_avg = 0;
			for (int i = 0; i < left_leaflist.size(); i++) {
				long double a[3] = { left_leaflist[i].m, left_leaflist[i].h, left_leaflist[i].w };
				for (int j = 0; j < right_leaflist.size(); j++) {
					long double b[3] = { right_leaflist[j].m, right_leaflist[j].h, right_leaflist[j].w };
					//cout << endl;
					fates_probablities = ProbabilitiesCalculation(g, a, b);
					/*if(isnan(fates_probablities[0])){
						cout << "P_subfunc=" << fates_probablities[0] << endl
								<< "P_cons=" << fates_probablities[1] << endl
								<< "P_neo=" << fates_probablities[2] << endl
								<< "P_pseudo=" << fates_probablities[3] << endl
								<< "P_spec=" << fates_probablities[4] << endl
								<< "ig_a=" << fates_probablities[5] << endl
								<< "ig_b=" << fates_probablities[6] << endl
								<< "ia_g=" << fates_probablities[7] << endl
								<< "ib_g=" << fates_probablities[8] << endl
								<< "ig_a_plus_b=" << fates_probablities[9] << endl
								<< "ia_plus_b_g=" << fates_probablities[10] << endl
								<< "pa=" << fates_probablities[11] << endl
								<< "pb=" << fates_probablities[12] << endl;
						cout << "done!" << endl;
					}*/
					sub_avg += fates_probablities[0];
					cons_avg += fates_probablities[1];
					new_avg += fates_probablities[2];
					pseudo_avg += fates_probablities[3];
					spec_avg += fates_probablities[4];
					dblneo_avg += fates_probablities[13];
					dups_list_file << root->pro.nodeid << "," << g[0] << "," << g[1] << "," << g[2] << ","
						<< left_leaflist[i].nodeid << "," << a[0] << "," << a[1] << "," << a[2] << ","
						<< right_leaflist[j].nodeid << "," << b[0] << "," << b[1] << "," << b[2] << ","
						<< fates_probablities[0] << "," << fates_probablities[1] << "," << fates_probablities[2] << "," << fates_probablities[3] << "," << fates_probablities[4] << "," << fates_probablities[13] << std::endl;
					//dups_list_file << endl;
				}
			}
			avg_dups_list_file << root->pro.nodeid << ","
				<< sub_avg / numofleaves << "," << cons_avg / numofleaves << "," << new_avg / numofleaves << "," << pseudo_avg / numofleaves << "," << spec_avg / numofleaves << "," << dblneo_avg / numofleaves << std::endl;

		}

	}

	for (int i = 0; i < root->GetNbChildren(); i++) {
		ClassifyDuplicationFates(root->children[i], dups_list_file, avg_dups_list_file);
		//ClassifyDuplicationFates(root->right, dups_list_file, avg_dups_list_file);
	}
}

void WritetoCSV(string proteinlist, long int gen) {

	vector<vector<string>> content;
	vector<string> row;
	string line, word;
	static std::ofstream dups_list_file("dups_fates_probablities.csv", std::ofstream::trunc);
	dups_list_file << "original node, g.m, g.h, g.w, first descendant id, a.m, a.h, a.w, second descendant id, b.m, b.h, b.w, P_subfunctionlization, P_conservation, P_newfunctionlizatoin, P_pseudogenization, P_specialization, P_dblneo " << std::endl;

	static std::ofstream avg_dups_list_file("dups_fates_avg_probablities.csv", std::ofstream::trunc);
	avg_dups_list_file << "original node, P_subfunctionlization, P_conservation, P_newfunctionlizatoin, P_pseudogenization, P_specialization, P_dblneo " << std::endl;

	//string infilename = "C:\\Users\\Manue\\Desktop\\tmp\\gene tree reconstuctor\\0.1 n n p.csv";
	string infilename = proteinlist;

	fstream file(infilename, ios::in);
	if (!file.is_open())
	{
		cout << "Could not open the file\n";
		return;
	}


	GeneTree gt;

	//ML: now the input csv is stramed and handled one line at a time to save memory
	Node* root = gt.construct_tree(file, gen);

	cout << "gene tree constructed" << endl;

	//cout << "left : " << root->left->pro.id <<  " right : \n" ;
	//printTree(root, nullptr, false);
	ClassifyDuplicationFates(root, dups_list_file, avg_dups_list_file);

	cout << "fates classified" << endl;

	auto it_root = root->children.begin();

	//ML: remove root children that only lead to one leaf
	while (it_root != root->children.end())
	{
		vector<protein> tmp;
		ReturnLeafNodes((*it_root), tmp);
		if (tmp.size() == 1)
		{
			it_root = root->children.erase(it_root);
		}
		else
		{
			++it_root;
		}
	}

	cout << "tree cleaned up" << endl;

	std::string newick_string;
	std::ofstream newick("newick.txt");
	newick_string = ToNewickString(root, true, true);
	newick << newick_string;
	return;
}
