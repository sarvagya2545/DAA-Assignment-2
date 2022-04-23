#include <iostream>
#include <vector>
using namespace std;

int n;

bool match(string s, int x, int y) {
  if (s[x - 1] == 'A' && s[y - 1] == 'U' ||
      s[x - 1] == 'G' && s[y - 1] == 'C' ||
      s[x - 1] == 'U' && s[y - 1] == 'A' ||
      s[x - 1] == 'C' && s[y - 1] == 'G') {
    return true;
  }

  return false;
}

void find_and_count_pairs(vector<vector<int>> &opt, string RNAsequence, vector<vector<int>> &split) {
  // k = window size
  for (int k = 5; k <= n - 1; k++) { 
    for (int i = 1; i <= n - k; i++) {
      int j = i + k;
      int temp = -1000000000;

      // for all t in [i , j - 5], if match (t, j) exists, let temp =
      // max(1 + opt[i][t - 1] + opt[t + 1][j - 1]) for all values of t
      int t_val = -1;
      for (int t = i; t <= j - 5; t++) {
        if (match(RNAsequence, t, j)) {
          int incl = 1 + opt[i][t - 1] + opt[t + 1][j - 1];
          if(temp < incl) {
            temp = incl;
            t_val = t;
          }
        }
      } 

      if(opt[i][j - 1] > temp) {
        opt[i][j] = opt[i][j - 1];
        // do not join
        split[i][j] = -1;
      } else {
        // join {t, j} splitting the RNA
        split[i][j] = t_val;
        opt[i][j] = temp;
      }
      // opt[i][j] = max(opt[i][j - 1], temp);
    }
  }
}

int pair_num = 1;

void print_pairs(vector<vector<int>> &split, string &RNAseq, int x = 1, int y = n) {
  if(abs(x - y) < 5) return;
  
  // pair<int,int> p = split[x][y];
  if(split[x][y] == -1) {
    print_pairs(split, RNAseq, x, y - 1);
  } else {
    int t = split[x][y];

    // print current pair
    cout << "pair number " << pair_num++ << " - ";
    cout << "(" << RNAseq[t - 1] << "," << RNAseq[y - 1] << ")\t";
    cout << "(" << t << "," << y << ")\n";
    
    print_pairs(split, RNAseq, x, t - 1);
    print_pairs(split, RNAseq, t + 1, y - 1);
  }
}

int main() {
  string RNAsequence;
  cin >> RNAsequence;
  n = RNAsequence.size();

  vector<vector<int>> opt = vector<vector<int>>(n + 1, vector<int>(n + 1, 0));  
  vector<vector<int>> split(n + 1, vector<int> (n + 1, -1)); 

  find_and_count_pairs(opt, RNAsequence, split);
  cout << "\nTotal number of pairs - " << opt[1][n] << "\n\n";
  print_pairs(split, RNAsequence, 1, n);
  return 0;
}
