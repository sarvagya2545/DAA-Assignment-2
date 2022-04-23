#include <iostream>
#include <vector>
using namespace std;

class RNA
{
public:
    int n;
    int pair_num = 1;
    RNA()
    {
    }

    /// @brief Function to tell whether or not 2 bases can form a bond in the RNA or not
    /// @param s the String representation of RNA
    /// @param x index value of base 1
    /// @param y index value of base 2
    /// @returns boolean value telling if the bases at index x & y can form a bond or not
    bool match(string s, int x, int y)
    {
        if (((s[x - 1] == 'A') && (s[y - 1] == 'U')) ||
            ((s[x - 1] == 'G') && (s[y - 1] == 'C')) ||
            ((s[x - 1] == 'U') && (s[y - 1] == 'A')) ||
            ((s[x - 1] == 'C') && (s[y - 1] == 'G')))
        {
            return true;
        }

        return false;
    }

    /// @brief Function to find and count the max number of bonds formed in the RNA
    /// @param opt DP Matrix where opt[i][j] = max number of bonds formed when the strand is RNA[i..j]
    /// @param RNAsequence String representation of the RNA
    /// @param split an empty matrix to be filled such that split[i][j] = t if (t,j) bond is formed, else split[i][j] = -1 if no bond from j is formed
    void find_and_count_pairs(vector<vector<int>> &opt, string &RNAsequence, vector<vector<int>> &split)
    {
        // k = window size
        for (int k = 5; k <= n - 1; k++)
        {
            for (int i = 1; i <= n - k; i++)
            {
                int j = i + k;
                int temp = -1000000000;

                // for all t in [i , j - 5], if match (t, j) exists, let temp =
                // max(1 + opt[i][t - 1] + opt[t + 1][j - 1]) for all values of t
                int t_val = -1;
                for (int t = i; t <= j - 5; t++)
                {
                    if (match(RNAsequence, t, j))
                    {
                        int incl = 1 + opt[i][t - 1] + opt[t + 1][j - 1];
                        if (temp < incl)
                        {
                            temp = incl;
                            t_val = t;
                        }
                    }
                }

                if (opt[i][j - 1] > temp)
                {
                    opt[i][j] = opt[i][j - 1];
                    // do not join
                    split[i][j] = -1;
                }
                else
                {
                    // join {t, j} splitting the RNA
                    split[i][j] = t_val;
                    opt[i][j] = temp;
                }
            }
        }
    }

    /// @brief Function to print all the bonds form between bases in a RNA strand
    /// @param split A grid where for 1 <= i < j <= n, split[i][j] = t or -1, t is the position where (t,j) is a bond
    /// @param RNAseq string representation for the input RNA
    /// @param x left index of the part of string considered
    /// @param y right index of the part of string considered
    void print_pairs(vector<vector<int>> &split, string &RNAseq, int x, int y)
    {
        if (abs(x - y) < 5)
            return;

        if (split[x][y] == -1)
        {
            print_pairs(split, RNAseq, x, y - 1);
        }
        else
        {
            int t = split[x][y];

            // print current pair
            cout << "pair number " << pair_num++ << " - ";
            cout << "(" << RNAseq[t - 1] << "," << RNAseq[y - 1] << ")\t";
            cout << "(" << t << "," << y << ")\n";

            print_pairs(split, RNAseq, x, t - 1);
            print_pairs(split, RNAseq, t + 1, y - 1);
        }
    }
};