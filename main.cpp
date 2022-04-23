#include <iostream>
#include <vector>
#include "RNA.h"
using namespace std;

int main() {
  string RNAsequence;
  cout << "Enter the RNA Sequence: \n";
  cin >> RNAsequence;
  RNA rna;
  int n = RNAsequence.size();
  rna.n = n;

  vector<vector<int>> opt = vector<vector<int>>(n + 1, vector<int>(n + 1, 0));  
  vector<vector<int>> split(n + 1, vector<int> (n + 1, -1)); 

  rna.find_and_count_pairs(opt, RNAsequence, split);
  cout << "\nTotal number of pairs - " << opt[1][n] << "\n\n";
  rna.print_pairs(split, RNAsequence, 1, n);
  return 0;
}
