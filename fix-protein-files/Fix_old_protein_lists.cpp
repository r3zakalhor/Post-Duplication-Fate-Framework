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
	int32_t start_pos_, length_, basal_level_, hamming_dist_, generation;
	string event, strand;
	bool ismapped;
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

	string GetID()
	{
		return to_string(pro.id);
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


class FixFile
{
public:
	//void init(QGraphicsScene* scene, QGraphicsView* view);
	void insert(int a);


	Node* fix_protein_lists(fstream& proteins_file, long int gen, string proteinfile) {

		// Convert fstream object to string
		std::ostringstream oss;
		oss << &proteins_file;
		std::string fileString = oss.str();
		cout << proteinfile << endl;
		// Create new filename by concatenating with "fixed_"
		std::string newFilename = "fixed_";
		newFilename = newFilename + proteinfile;
		cout << newFilename << endl;
		// Open the new file
		std::ofstream fixed_proteins_list;
		fixed_proteins_list.open(newFilename, std::ofstream::trunc);
		int in;
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
		std::streampos lastposition;
		string line;
		string event = "initial";
		int lastpos;
		getline(proteins_file, line);	//skip the header line
		fixed_proteins_list << line << std::endl;
		bool firstline = false;
		while (getline(proteins_file, line) && time_gen <= gen)
		{
			stringstream str(line);
			vector<string> row;
			string word;
			while (getline(str, word, ','))
				row.push_back(word);
			//cout << row[0] << endl;
			if (row[0].size() > 1 && row[0] != "INV" && row[0] != "TRANS") {
				//cout << row[0] << endl;
				lastposition = fixed_proteins_list.tellp();
				fixed_proteins_list << line << std::endl;
				firstline = true;
				//cout << line << endl;
				numofpro = stoi(row[8]);
				//cout << "event type " << row[0] << " at gen: " << row[15] << endl;
				//cout << "num of genes " << stoi(row[8]) << endl;

				Node* next = new Node();
				next->pro.nodeid = numofnodes;
				numofnodes++;
				next->dup = false;
				next->pro.event = row[0];
				event = row[0];
				next->pro.id = stoi(row[1]);
				next->pro.parent_id = stoi(row[2]);
				next->pro.start_pos_ = stoi(row[3]);
				next->pro.length_ = stoi(row[4]);
				next->pro.m = stold(row[10]);
				next->pro.w = stold(row[11]);
				next->pro.h = stold(row[12]);
				next->pro.ismapped = false;
				lastpos = next->pro.start_pos_;
				//next->pro.strand = row[15];
				next->pro.generation = stoi(row[15]);
				time_gen = stoi(row[15]);
				next_level.push_back(next);
				i_numofpro++;
				//cout << "i_numofpro: " << i_numofpro << endl;
			}
			else if (i_numofpro < numofpro && event != "INV" && row[0].size() <= 1 && event != "TRANS") {
				//cout << "2";

				//TODO: evil duplicated code
				Node* next = new Node;
				next->pro.nodeid = numofnodes;
				numofnodes++;
				next->dup = false;
				next->pro.event = event;
				next->pro.id = stoi(row[1]);
				next->pro.parent_id = stoi(row[2]);
				next->pro.start_pos_ = stoi(row[3]);
				next->pro.length_ = stoi(row[4]);
				next->pro.m = stold(row[10]);
				next->pro.w = stold(row[11]);
				next->pro.h = stold(row[12]);
				next->pro.ismapped = false;

				//next->pro.strand = row[15];
				next->pro.generation = stoi(row[15]);
				time_gen = stoi(row[15]);

				i_numofpro++;
				//cout << "i_numofpro: " << i_numofpro << endl;
				if (lastpos != next->pro.start_pos_) {
					next_level.push_back(next);
					lastposition = fixed_proteins_list.tellp();
					fixed_proteins_list << line << std::endl;
					firstline = false;

					// Move the write pointer to the beginning of the line you just wrote
				}
				else {
					for (int n = 0; n < current_level.size(); n++) {
						if (current_level[n]->pro.id == next->pro.parent_id) {
							fixed_proteins_list.seekp(lastposition);
							if (firstline) {
								row[0] = event;
								firstline = false;
							}
							for (size_t i = 0; i < row.size(); i++) {
								fixed_proteins_list << row[i];
								if (i < row.size() - 1) {
									fixed_proteins_list << ",";
								}
							}
							fixed_proteins_list << endl;
							next_level.pop_back();
							next_level.push_back(next);
							break;
						}
					}
				}
				lastpos = next->pro.start_pos_;
			}



			if (row[0].size() > 1 && (row[0] == "INV" || row[0] == "TRANS")) {
				//cout << row[0] << endl;
				//cout << "3";
				//fixed_proteins_list << line << std::endl;
				vector<Node*> findnodes;
				numofpro = stoi(row[8]);
				//cout << "event type " << row[0] << " at gen: " << row[15] << endl;
				//cout << "num of genes " << stoi(row[8]) << endl;
				in = 0;
				Node* next = new Node();
				next->pro.nodeid = numofnodes;
				numofnodes++;
				next->dup = false;
				next->pro.event = row[0];
				event = row[0];
				next->pro.id = stoi(row[1]);
				next->pro.parent_id = stoi(row[2]);
				next->pro.start_pos_ = stoi(row[3]);
				lastpos = stoi(row[3]);
				next->pro.length_ = stoi(row[4]);
				next->pro.m = stold(row[10]);
				next->pro.w = stold(row[11]);
				next->pro.h = stold(row[12]);
				next->pro.ismapped = false;
				//next->pro.strand = row[15];
				next->pro.generation = stoi(row[15]);
				time_gen = stoi(row[15]);
				next_level.push_back(next);
				i_numofpro++;
				for (int n = 0; n < current_level.size(); n++) {
					if (current_level[n]->pro.m == next->pro.m && current_level[n]->pro.w == next->pro.w && current_level[n]->pro.start_pos_ == next->pro.start_pos_) {
						findnodes.push_back(current_level[n]);
					}
					else if (current_level[n]->pro.m == next->pro.m && current_level[n]->pro.w == next->pro.w && current_level[n]->pro.length_ == next->pro.length_) {
						findnodes.push_back(current_level[n]);
					}
				}
				bool find = false;
				for (int n = 0; n < findnodes.size(); n++) {
					if (findnodes[n]->pro.id == next->pro.id) {
						row[2] = to_string(findnodes[n]->pro.id);
						find = true;
						findnodes[n]->pro.ismapped = true;
						break;
					}
				}
				if (find) {
					for (size_t i = 0; i < row.size(); i++) {
						fixed_proteins_list << row[i];
						if (i < row.size() - 1) {
							fixed_proteins_list << ",";
						}
					}
					fixed_proteins_list << endl;
				}
				else {
					if (!findnodes.empty()) {
						for (int n = 0; n < findnodes.size(); n++) {
							if (!findnodes[n]->pro.ismapped) {
								row[2] = to_string(findnodes[n]->pro.id);
								findnodes[n]->pro.ismapped = true;
								break;
							}
						}
						for (size_t i = 0; i < row.size(); i++) {
							fixed_proteins_list << row[i];
							if (i < row.size() - 1) {
								fixed_proteins_list << ",";
							}
						}
						fixed_proteins_list << endl;
					}
					else {
						bool match = false;
						int max = -1;
						for (int n = 0; n < current_level.size(); n++) {
							if (max < current_level[n]->pro.id)
								max = current_level[n]->pro.id;
							if (current_level[n]->pro.id == next->pro.parent_id) {
								match = true;
							}
						}
						if (!match)
							fixed_proteins_list << line << std::endl;
						else {
							in++;
							max = max + in;
							row[2] = to_string(max);
							for (size_t i = 0; i < row.size(); i++) {
								fixed_proteins_list << row[i];
								if (i < row.size() - 1) {
									fixed_proteins_list << ",";
								}
							}
							fixed_proteins_list << endl;
						}
					}
				}
				//cout << "i_numofpro: " << i_numofpro << endl;
			}
			else if (i_numofpro < numofpro && (event == "INV" || event == "TRANS") && row[0].size() <= 1) {
				//cout << "4";
				//fixed_proteins_list << line << std::endl;
				vector<Node*> findnodes;
				//TODO: evil duplicated code
				Node* next = new Node;
				next->pro.nodeid = numofnodes;
				numofnodes++;
				next->dup = false;
				next->pro.event = event;
				next->pro.id = stoi(row[1]);
				next->pro.parent_id = stoi(row[2]);
				next->pro.start_pos_ = stoi(row[3]);
				next->pro.length_ = stoi(row[4]);
				next->pro.m = stold(row[10]);
				next->pro.w = stold(row[11]);
				next->pro.h = stold(row[12]);
				next->pro.ismapped = false;
				//next->pro.strand = row[15];
				next->pro.generation = stoi(row[15]);
				time_gen = stoi(row[15]);

				i_numofpro++;
				//cout << "i_numofpro: " << i_numofpro << endl;
				if (lastpos != next->pro.start_pos_) {
					next_level.push_back(next);
					if (time_gen == 7611 && next->pro.id == 1038)
						cout << " 1) " << next->pro.id << ", " << next->pro.m << ", " << next->pro.h << ", " << next->pro.w << endl;
					for (int n = 0; n < current_level.size(); n++) {
						//if (time_gen == 10 && (current_level[n]->pro.id == 122 || current_level[n]->pro.id == 122) && next->pro.id == 19)
							//cout << " 2) " << current_level[n]->pro.id << ", " << current_level[n]->pro.m << ", " << current_level[n]->pro.h << ", " << current_level[n]->pro.w << ", " << current_level[n]->pro.ismapped << endl;
						if (current_level[n]->pro.m == next->pro.m && current_level[n]->pro.w == next->pro.w && current_level[n]->pro.start_pos_ == next->pro.start_pos_) {
							findnodes.push_back(current_level[n]);
						}
						else if (current_level[n]->pro.m == next->pro.m && current_level[n]->pro.w == next->pro.w && current_level[n]->pro.length_ == next->pro.length_) {
							findnodes.push_back(current_level[n]);
						}
					}
					if (next->pro.id == 1038 && time_gen == 7611)
						cout << "size: " << findnodes.size() << endl;
					bool find = false;
					for (int n = 0; n < findnodes.size(); n++) {
						if (next->pro.id == 1038 && time_gen == 7611)
							cout << n << ": " << to_string(findnodes[n]->pro.id) << " " << findnodes[n]->pro.ismapped << endl;
						if (findnodes[n]->pro.id == next->pro.id && !findnodes[n]->pro.ismapped) {
							row[2] = to_string(findnodes[n]->pro.id);
							find = true;
							//if (findnodes[n]->pro.id == 1038)
								//cout << next->pro.id << " at gen " << time_gen << endl;
							findnodes[n]->pro.ismapped = true;
							break;
						}
					}
					if (find) {
						for (size_t i = 0; i < row.size(); i++) {
							fixed_proteins_list << row[i];
							if (i < row.size() - 1) {
								fixed_proteins_list << ",";
							}
						}
						fixed_proteins_list << endl;
					}
					else {
						if (!findnodes.empty()) {
							for (int n = 0; n < findnodes.size(); n++) {
								if (!findnodes[n]->pro.ismapped) {
									//if (next->pro.id == 19 && time_gen == 10)
										//cout << "is not mapped: " << to_string(findnodes[n]->pro.id) << endl;
									row[2] = to_string(findnodes[n]->pro.id);
									//if (findnodes[n]->pro.id == 243)
										//cout << next->pro.id << " at gen " << time_gen << endl;
									findnodes[n]->pro.ismapped = true;
									break;
								}
							}
							for (size_t i = 0; i < row.size(); i++) {
								fixed_proteins_list << row[i];
								if (i < row.size() - 1) {
									fixed_proteins_list << ",";
								}
							}
							fixed_proteins_list << endl;
						}
						else {
							bool match = false;
							int max = -1;
							for (int n = 0; n < current_level.size(); n++) {
								if (max < current_level[n]->pro.id)
									max = current_level[n]->pro.id;
								if (current_level[n]->pro.id == next->pro.parent_id) {
									match = true;
								}
							}
							if (!match)
								fixed_proteins_list << line << std::endl;
							else {
								in++;
								max = max + in;
								row[2] = to_string(max);
								for (size_t i = 0; i < row.size(); i++) {
									fixed_proteins_list << row[i];
									if (i < row.size() - 1) {
										fixed_proteins_list << ",";
									}
								}
								fixed_proteins_list << endl;
							}
						}
					}
				}
				lastpos = next->pro.start_pos_;
			}
			if (i_numofpro == numofpro) {
				current_level.clear();
				current_level = next_level;
				next_level.clear();
				//cout << current_level.size() << endl;
				i_numofpro = 0;
				event = "initial";
			}

		}
		cout << "last event at generation: " << gen << endl;
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


int main(int argc, char* argv[]) {

	//cout << "You have entered " << argc
		//<< " arguments:" << "\n";
	string s;
	string version;
	string proteinfile;
	long int gen;
	vector<long double> fates_variables;
	fates_variables.resize(14);
	long double g[3], a[3], b[3];
	proteinfile = "proteins_list_after_events.csv";
	gen = 58;
	string infilename = proteinfile;

	fstream file(infilename, ios::in);
	if (!file.is_open())
	{
		cout << "Could not open the file\n";
		return 0;
	}
	FixFile ft;
	Node* root = ft.fix_protein_lists(file, gen, proteinfile);
	//cout << "root child: " << root->GetNbChildren() << endl;
	//printTree(root, "tree_structure.txt", "tree_newick.txt");
	cout << "file is fixed!" << endl;


	return 0;
}
