#ifndef LK_MATRIX
#define LK_MATRIX

#include <vector>

using namespace std;

class LKSolver {
  public:
    int size;
    LKSolver(vector<pair<double, double> > &coords, vector<int> &ids);
    vector<int> getCurrentTour();
    double getCurrentTourDistance();
    void optimizeTour();
    vector<pair<double, double> > getTour();
    void printTour();
    void printTourIds(bool showCoords=false);


  private:
    vector<int> tour;
    vector<vector<int> > edgeFlags;
    vector<pair<double, double> > coords;
    vector<int> ids;
    void joinLocations(int i, int j);
    vector<vector<double> > edgeDistances;
    void LKMove(int tourStart);
    void reverse(int start, int end);
    bool isTour();
};

#endif
