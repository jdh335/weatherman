#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <limits>
#include <vector>
#include <map>
#include <tuple>

using namespace std;

#define MATCH -3
#define INDEL 5
#define SUB 1

struct AlignResult {
    double align_cost;
    string seq1;
    string seq2;
};

struct ScoreResult {
    double score;
    int row;
    int col;
};

class GeneSequencing {
public:
    string string1, string2;

    AlignResult align(string seq1, string seq2, bool banded, int align_length) {

        ScoreResult** score_matrix = new ScoreResult*[align_length+1];
        for (int i = 0; i <= align_length; ++i) {
            score_matrix[i] = new ScoreResult[align_length+1];
        }

        string1 = "-" + seq1.substr(0, align_length);
        string2 = "-" + seq2.substr(0, align_length);

        for (int i = 0; i <= align_length; i++) {
            for (int j = 0; j <= align_length; j++) {
                double score = numeric_limits<double>::infinity();
                ScoreResult result = {score, 0, 0};

                double left = (i==0) ? numeric_limits<double>::infinity() : (score_matrix[i-1][j]).score + INDEL;
                double top = (j==0) ? numeric_limits<double>::infinity() : (score_matrix[i][j-1]).score + INDEL;
                double diagonal = (i==0 || j==0) ? numeric_limits<double>::infinity(): (score_matrix[i-1][j-1]).score + (string1[i] == string2[j]) ?  MATCH : SUB;

                if (left < score) {
                    score = left;
                    result = {left, i-1, j};
                }
                if (top < score) {
                    score = top;
                    result = {top, i, j-1};
                }
                if (diagonal < score) {
                    score = diagonal;
                    result = {diagonal, i-1, j-1};
                }
                
                score_matrix[i][j] = result;
            }
        }
        // editDistance();
        ScoreResult final = score_matrix[10][10];
        return {0,0,0};
    }

    // void editDistance() {
    //     for (int i = 0; i <= string1.length(); i++) {
    //         for (int j = 0; j <= string2.length(); j++) {
    //             score_matrix[i][j] = computeScore(i, j, string1[i], string2[j]);
    //         }
    //     }
    // }

    // ScoreResult computeScore(int i, int j, char char1, char char2) {
    //     ScoreResult result;
    //     double score = numeric_limits<double>::infinity();

    //     double left = (score_matrix[i-1][j]).score + INDEL;
    //     double top = (score_matrix[i][j-1]).score + INDEL;
    //     double diagonal = (score_matrix[i-1][j-1]).score + (char1 == char2) ?  MATCH : SUB;

    //     if (left < score) {
    //         score = left;
    //         result = {left, i-1, j};
    //     }
    //     if (top < score) {
    //         score = top;
    //         result = {top, i, j-1};
    //     }
    //     if (diagonal < score) {
    //         score = diagonal;
    //         result = {diagonal, i-1, j-1};
    //     }
    //     return result;
    // }
    
};


// Main driver function
int main() {
    string path = "./GenomeSequences/";
    GeneSequencing gs;
    int align_length = 10;
    string seq[8];
    string line[8];
    ifstream file1(path+"BCov-ENT.txt");
    while (getline(file1, line[0])) seq[0] += line[0];
    ifstream file2(path+"BCov-LUN.txt");
    while (getline(file2, line[1])) seq[1] += line[1];
    ifstream file3(path+"BCov-Mebus.txt");
    while (getline(file3, line[2])) seq[2] += line[2];
    ifstream file4(path+"BCov-Quebec.txt");
    while (getline(file4, line[3])) seq[3] += line[3];
    ifstream file5(path+"Mouse_Hepatitis.txt");
    while (getline(file5, line[4])) seq[4] += line[4];
    ifstream file6(path+"Murine_hepatitis1.txt");
    while (getline(file6, line[5])) seq[5] += line[5];
    ifstream file7(path+"Murine_hepatitis2.txt");
    while (getline(file7, line[6])) seq[6] += line[6];
    ifstream file8(path+"Murine_hepatitis3.txt");
    while (getline(file8, line[7])) seq[7] += line[7];


    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 8; j++) {
            AlignResult result = gs.align(seq[i], seq[j], false, align_length);
            cout << "Sequence " << i << " and Sequence " << j << " Result: " << result.align_cost << endl;
        }
    }


    return 0;
}