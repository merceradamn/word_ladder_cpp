/**************************************************************************************************
Name:		Adam Mercer
Project:	Word Ladder Problem
Descript:	The Word Ladder Problem: Given a start word and an end word, find a 
			shortest sequence of words (shortest word ladder) from the start to
			the end so that:
			1. Each word is a legal word from a given dictionary
			2. The next word is obtained from the previous word by substituting
				exactly one letter.
			
			If there is no such sequence, then say so. Note: Assume all words are
			in lowercase.
			Example: Given a start word zero and an end word five, here is the
			ladder:

			zero => hero => here => hire => fire => five
			We assume that all the words are legal.

Notes:		Part A is done. Part B is done.

**************************************************************************************************/

#include <iostream>			// std::cout
#include <fstream>			// std::ifstream
#include <iomanip>			// Formatting library
#include <string>			// String types, character traits and a set of converting functions
#include <cmath>			// Set of functions to compute common math operations & transformations
#include <cstdio>			// C-style I/O
#include <chrono>			// Deals with time
#include <vector>			// Vector library; expandable array container
#include <queue>			// Queue library
#include <map>				// Sorted map library
#include <unordered_map>	// Bucket insert map library
#include <stack>

using namespace std;	// So I don't have to type std:: on everything
typedef chrono::high_resolution_clock Clock;

// Vertex object for use in holding the vertex table
struct Vertex {
	// 0 = unknown, 1 = frontier, 2 = visited;
	int status;											// Holds the visit status
														// Should be -1 if this is starting vertex
	int distance;										// Distance to node from 
	string parent;										// String of parent vertex
};

// Vertex object used for an adjacency list
struct VertexEdge {
	string vert;										// String name of the edge
	int weight;											// Weight of the edge
};

// Used to initialize a vertex to base values
void initVertex(Vertex &v) {
	v.status = 0;										// "0" = unknown
	v.distance = 9999;									// "" = INF distance
	v.parent = "|";										// Name of the parent vertex
}

void initEdge(VertexEdge v) {
	v.vert = "NULL";
	v.weight = -1;
}

// Set the edge variables
void setEdge(VertexEdge &v, string &s, int& w) {
	v.vert = s;
	v.weight = w;
}

bool oneCharOff(const string & word1, const string & word2);
int oneCharOffSaver(const string & word1, const string & word2);
void filePrint(fstream& f);
map<string, vector<string>> computeAdjacentWords(const vector<string> & words);
unordered_map<string, string> findChain(const map<string, vector<string>> & adjacentWords,
	const string & first, const string & second);
vector<string> getChainFromPreviousMap(const unordered_map<string, string> & previous, 
	const string & second);
int computeWeightedPath(map<string, vector<VertexEdge>> &wwm, map<string, vector<Vertex>> & vt,
	const string & first, const string & second);
stack<string> getChainFromWeighted(map<string, vector<Vertex>> & vt, const string & first, const string & second);

int main() {

/*************************************************************************************************/
// VARIABLES AND OTHER DECLARATIONS
/*************************************************************************************************/
// Part A variables
	vector<string> words;								// Holds a vector with words right length
	map<string, vector<string>> wordMap;				// Holds map of words and adjacent words
	unordered_map<string, string> path;					// Creates a map to hold our path
	vector<string> res;									// Holds a vector that stores our ladder
	bool done = false;

// Part B variables
	map<string, vector<string>>::iterator WM;			// Iterator to our word map
	string w1, w2;										// Holds starting word and comparing
	char orig, repl;									// Original and Replacing char
	int LWeight, pos;									// Holds weight and char position
	int temp = 0;										// Temp variable
	int numbers[26][26];								// Holds letter weight grid
	char number;										// Used in letterweights section

// Testing
	map<string, vector<VertexEdge>> weightedWordMap;	// Weighted word map
	map<string, vector<Vertex>> wordTable;				// Table for holding vertex data
	map<string, vector<VertexEdge>>::iterator WWM;		// Iterator for weighted word map
	Vertex vertTemp;
	VertexEdge word;

// "Global" variables
	string top, bot;									// Strings to hold word pair
	string readIn;										// String to hold compared words
	unsigned __int64 times[4] = {0};					// Time holder
	auto start = Clock::now(), end = Clock::now();		// Holds times for calculating duration

// File-handling Section
	ifstream fileIn;									// Read in our dictionary
	string dict = "dictionary.txt";						// "dictionary.txt"
	fileIn.open(dict);									// Open the dictionary
	if (!fileIn.is_open())
		cout << "Dictionary didn't open." << endl;

	ifstream weights;									// Weight graph file
	string weight = "letterWeights.txt";				// "letterWeights.txt
	weights.open(weight);								// Open the weight file
	if (!weights.is_open())
		cout << "Weights didn't open." << endl;

/*************************************************************************************************/
// Main Program Section
/*************************************************************************************************/
//	LOAD THE DICTIONARY AND THE WORD MAP (ADJACENCY LIST)
	start = Clock::now();								// Start clock for time duration
	bool tVal = 0, bVal = 0;							// Values to check if words are good
	while (!fileIn.eof()) {								// While end-of-file indicator is false
		fileIn >> readIn;								// Read in the word
		words.push_back(readIn);						// Push the word onto the vector
	}
	fileIn.close();										// Close the file
	end = Clock::now();									// Stop clock and print duration
	times[0] = chrono::duration_cast<chrono::milliseconds>(end - start).count();
	cout << "Loaded the dictionary in " << times[0] << "ms." << endl;

	start = Clock::now();								// Start clock for time duration
	wordMap = computeAdjacentWords(words);				// Build adjacency list, store in wordMap
	end = Clock::now();									// Stop clock and print duration
	times[1] = chrono::duration_cast<chrono::milliseconds>(end - start).count();
	cout << "Loaded the word map in " << times[1] << "ms." << endl;

// Gets the raw data from the weights chart and stores in a two D array for processing
	for (int i = 0; i < 26; i++) {
		for (int j = 0; j < 26; j++) {
			weights >> number;
			if (number == '-') {
				numbers[i][j] = -1;
			}
			else
				numbers[i][j] = number - 48;
		}
	}
	weights.close();

//	WORD PAIR RETRIEVAL AND EVALUALTION
	while (!done) {
		// Gets word pair, checks for valid, then adds all equal length words to queue to read
		/* Get the word pair and store it */
		cout << "Start> ";	cin >> top;						// Get the start word
		cout << "End> ";	cin >> bot;						// Get the end word

		if (top.length() != bot.length()) {					// Check word lengths for match
			cout << "Length mismatch." << endl;
			continue;
		}

		/* Checks if either word doesn't exist in the dictionary. Returns if word is missing */
		fileIn.open(dict);
		tVal = false; bVal = false;
		while (!fileIn.eof()) {								// While end-of-file indicator is false
			fileIn >> readIn;
			if (top == readIn) {							// Set top valid flag if match
				tVal = true; 
			}				
			if (bot == readIn) {							// Set bot valid flag if match
				bVal = true; 
			}				

		}
		fileIn.close();										// Close the file we have our words

		if (tVal == 0 || bVal == 0) {
			if (tVal == 0)
				cout << "Start word doesn't exist." << endl;

			if(bVal == 0)
				cout << "End word doesn't exist." << endl;

			continue;
		}

// LADDER CALCULATION AND DELIVERY
		start = Clock::now();								// Start clock for time duration
		path = findChain(wordMap, top, bot);				// Calls the function to fill our path
		res = getChainFromPreviousMap(path, bot);			// Gets the chain if it exists
		end = Clock::now();									// Stop clock, print duration

		if (res.size() > 1) {								// If we have a solvable ladder
			cout << "Found the unweighted ladder in " << res.size() << " steps." << endl;
			times[2] = chrono::duration_cast<chrono::milliseconds>(end - start).count();
			cout << "Found in " << times[2] << "ms." << endl;
			for (unsigned int j = 0; j < res.size(); j++) {
				cout << res[j] << endl;						// Print the word
			}
			cout << endl;
		}
		else {
			cout << "No ladder is possible." << endl;
			continue;
		}

/*************************************************************************************************/
// Weighted Ladder Code
/*************************************************************************************************/
// Read all the words from wordMap and calculate the weights and adding them
		WM = wordMap.begin();								// Get an iterator for wordMap
		initEdge(word);										// Initialize dummy vertex edge
		int x = 0;
		while (WM != wordMap.end()) {						// While we have words to copy
			if (top.length() == WM->first.length()) {		// Compare top with each word
				w1 = WM->first;								// Get first word
				wordTable[WM->first].push_back(vertTemp);	// Push a vertex to the table
				initVertex(wordTable[WM->first][0]);		// Set data in vertex to default

				// Get every vector for the specific key
				for (unsigned int i = 0; i < WM->second.size(); i++) {
					// Calculate the edge weight
					w2 = WM->second[i];						// Get second word
					pos = oneCharOffSaver(w1, w2);			// Compare the words and return pos
					// If pos returned greater than -1 we have a good pair
					if (pos > -1) {
						orig = w1.at(pos);					// Get the original character
						repl = w2.at(pos);					// Get the replacing character

						orig -= 97; repl -= 97;				// Offset for ASCII values

						// If we have a letter within our range
						if ((orig >= 0 && orig <= 25) && (repl >= 0 && repl <= 25)) {
							LWeight = numbers[orig][repl];	// Get the weight for the swap
						}
						else {
							LWeight = 10;					// Set the weight to 0
						}
					}

					// Add the word and weight to the new map
					weightedWordMap[WM->first].push_back(word);
					setEdge(weightedWordMap[WM->first][i], WM->second[i], LWeight);
				}
			}
			WM++;											// Move the iterator
		}

		// Call Dijkstra's
		start = Clock::now();								// Start clock for time duration
		temp = computeWeightedPath(weightedWordMap, wordTable, top, bot);
		end = Clock::now();									// Stop clock, print duration

		if (temp == 0) {
			times[3] = chrono::duration_cast<chrono::milliseconds>(end - start).count();
			cout << "Found the weighted ladder in " << times[3] << "ms. " << endl;

			// Build the path and print it
			getChainFromWeighted(wordTable, top, bot);
			cout << endl;
		}
		// Clear the variables
		path.clear();
		res.clear();
		weightedWordMap.clear();
		wordTable.clear();
	}

	return 0;
}

// Returns true if word1 and word2 are the same length
// and differ by only one character
bool oneCharOff(const string & word1, const string & word2) {
	if (word1.length() != word2.length())
		return false;

	int diffs = 0;

	for (unsigned int i = 0; i < word1.length(); ++i)
		if (word1[i] != word2[i])
			if (++diffs > 1)
				return false;

	return diffs == 1;
}

// Returns the position in the word
int oneCharOffSaver(const string & word1, const string & word2) {
	if (word1.length() != word2.length())
		return -1;

	int pos = -1;
	int diffs = 0;

	for (unsigned int i = 0; i < word1.length(); ++i)
		if (word1[i] != word2[i])
			if (++diffs > 1)
				return -1;
			else
				pos = i;

	return pos;
}

// Assumes a file is open in read mode and prints the file
// Prints each string on a separate line
void filePrint(fstream& f) {
	string data;
	if (f.is_open()) {
		while (!f.eof()) {
			f >> data;
			cout << data << endl;
		}
	}
	else
		cout << "Can't print file. Not open." << endl;

	return;
}

// Computes a map in which the keys are words and values are vectors of words
// that differ in only one character from the corresponding key.
// Uses an efficient algorithm that is O(N log N) with a map
map<string, vector<string>> computeAdjacentWords(const vector<string> & words) {
	map<string, vector<string>> adjWords;
	map<int, vector<string>> wordsByLength;

	// Groups elements from words by length in wordsByLength
	// Auto detects what type variable we need and then iterates
	// over the entire words vector
	for (auto & str : words)
		wordsByLength[str.length()].push_back(str);

	// Work on each group separately
	// For each group containing words of length groupNum
	for (auto & entry : wordsByLength) {
		const vector<string> & groupsWords = entry.second;
		int groupNum = entry.first;

		// Work on each position in each group
		// For each position (ranging from 0 to len-1)
		for (int i = 0; i < groupNum; ++i) {
			// Remove one character in specified position, computing representative.
			// Words with same representatives are adjacent; so populate a map ...

			// Make an empty map<string,vector<string>> repsToWords for each word
			map<string, vector<string>> repToWord;

			// Obtain word's representative by removing position p
			for (auto & str : groupsWords) {
				string rep = str;
				rep.erase(i, 1);
				repToWord[rep].push_back(str);
			}

			// and then look for map values with more than one string
			// Use cliques in repsToWords to update adjWords map
			for (auto & entry : repToWord) {
				const vector<string> & clique = entry.second;
				if (clique.size() >= 2)
					for (unsigned int p = 0; p < clique.size(); ++p)
						for (unsigned int q = p + 1; q < clique.size(); ++q) {
							adjWords[clique[p]].push_back(clique[q]);
							adjWords[clique[q]].push_back(clique[p]);
						}
			}
		}
	}
	return adjWords;
}

// Runs the shortest path calculation from the adjacency map, returning a vector
// that contains the sequence of word changes to get from first to second.
unordered_map<string, string> findChain(const map<string, vector<string>> & adjacentWords,
	const string & first, const string & second) {
	unordered_map<string, string> previousWord;
	queue<string> q;

	q.push(first);

	while (!q.empty()) {
		string current = q.front(); q.pop();
		auto itr = adjacentWords.find(current);

		vector<string> adj = itr->second;
		for (string & str : adj)
			if (previousWord[str] == "") {
				previousWord[str] = current;
				q.push(str);
			}
	}
	previousWord[first] = "";

	return previousWord;
}

// After the shortest path calculation has run, computes the vector that
// contains the sequence of words changes to get from first to second.
vector<string> getChainFromPreviousMap(const unordered_map<string, string> & previous,
	const string & second) {

	vector<string> result;
	auto & prev = const_cast<unordered_map<string, string> &>(previous);

	for (string current = second; current != ""; current = prev[current])
		result.push_back(current);

	reverse(begin(result), end(result));
	return result;
}

// Pass the weighted word map and the vert table and word pair
// Calculate the table for path calculation
int computeWeightedPath(map<string, vector<VertexEdge>> &wwm, map<string, vector<Vertex>> & vt,
	const string & first, const string & second) {

	// Iterators for the problem
	map<string, vector<VertexEdge>>::iterator wwmIT = wwm.begin();
	map<string, vector<Vertex>>::iterator vtIT = vt.begin();
	map<string, vector<Vertex>>::iterator vtCheck = vtIT;

	stack<string> path;

	// Set all the data in the table to defaults (unknown/inf cost/-1 parent)
	for (unsigned int i = 0; i < vt.size(); i++) {
		for (unsigned int j = 0; j < vtIT->second.size(); j++) {
			initVertex(vtIT->second[j]);
		}
		vtIT++;
	}

	// Get to the starting vertex and set the stage for dijkstra's
	vtIT = vt.find(first);
	wwmIT = wwm.find(first);
	if (vtIT == vt.end() || wwmIT == wwm.end()) {
		cout << "Start word not found. Returning..." << endl;
		return -1;
	}

	// Set the parameters for our starting vertex
	vtIT->second[0].distance = 0;		// Distance to this one is zero

	string lowest;						// holds name of the current string with the lowest weight
	string current = vtIT->first;		// Set current to what word we're starting on
	int weight = 0;
	bool fail = 0;

	do{
		// Move to the current vertex in weightmap and vert table
		vtIT = vt.find(current);
		wwmIT = wwm.find(current);

		vtIT->second[0].status = 1;		// Set the vertex to known since we're here

		// Check each edge for current vertex for weight and status
		for (unsigned int i = 0; i < wwmIT->second.size(); i++) {
			vtCheck = vt.find(wwmIT->second[i].vert);
			if (vtCheck->second[0].status != 1) {
				weight = vtIT->second[0].distance + wwmIT->second[i].weight;
				
				// If weight is less than the weight at the vertex neighbor then change it
				if (weight < vtCheck->second[0].distance) {
					vtCheck->second[0].distance = weight;	// Change weight since found cheaper
					vtCheck->second[0].parent = vtIT->first;// Change parent to current vertex
				}
			}
		}
		
		// Go through the vert table and store the lowest weight
		// When we finish reading the table move to the lowest weight
		vtIT = vt.begin();									// Set iterator to begin of vert table
		weight = 999;
		while (vtIT != vt.end()) {							// Look for unknown vert w/ lowest cost
			if (vtIT->second[0].distance < weight && vtIT->second[0].status != 1) {
				// Set current to this vertex and compare the rest
				current = vtIT->first;
				weight = vtIT->second[0].distance;
			}
			vtIT++;
		}

		// If we didn't find a lower weight we are done and should return
		if (weight == 999) {
			fail = 1;
		}

	} while (fail == 0);	// If we can't find a frontier node then we're done
	
	return 0;
}

stack<string> getChainFromWeighted(map<string, vector<Vertex>> & vt, const string & first ,const string & second) {
	// Start from the end word and ride the parents up to the start word
	stack<string> path;
	int pathCost = 0;
	map<string, vector<Vertex>>::iterator vtIT;
	map<string, vector<Vertex>>::iterator vtStart = vt.find(first);
	vtIT = vt.find(second);

	// Check to see if we can find the second word in the 
	if (vtIT == vt.end()) {
		cout << "End word doesn't exist." << endl;
		return path;
	}

	// Get the path from the end word to the start and push onto a stack
	pathCost = vtIT->second[0].distance;
	while (vtIT->second[0].parent != "|") {
		path.push(vtIT->first);
		vtIT = vt.find(vtIT->second[0].parent);
		if (vtIT == vt.end()) {
			cout << "Error getting next parent." << endl;
			return path;
		}
	}

	if (vtIT->second[0].parent == "|") {
		path.push(vtIT->first);
	}

	// Print the path in correct order
	cout << "Found the shortest path (Dijkstra's). Path costs " << pathCost << ".\n";
	while (!path.empty()) {
		cout << path.top() << endl;
		path.pop();
	}

	return path;
}