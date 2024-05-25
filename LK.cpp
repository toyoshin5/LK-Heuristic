#include "LKSolver.h"
#include <vector>
#include <cmath>
#include <set>
#include <iostream>
#include <cassert>
#include <cstdlib>
#include <ctime>
#include <sys/time.h>
#include <stdio.h>
#include <unistd.h>

using namespace std;

pair<int, int> make_sorted_pair(int x, int y) {
  if (x < y) {
    return make_pair(x, y);
  } else {
    return make_pair(y, x);
  }
}

LKSolver::LKSolver(vector<pair<double, double> > &cds, vector<int> &idss) {
  this->coords = cds;
  this->ids = idss;
  size = ids.size();

  // 経路を初期化
  tour = vector<int>(size, 0);

  // ランダムな経路で初期化 
  for (int i = 0; i < size; i++) {
    tour[i] = (i + 1) % size;
  }

  // coordsをランダムにシャッフル
  srand(time(NULL));
  for (int i = 0; i < size; i++) {
    int j = rand() % size;
    pair<double, double> temp = coords[i];
    coords[i] = coords[j];
    coords[j] = temp;
    int tempId = ids[i];
    ids[i] = ids[j];
    ids[j] = tempId;
  }
  // 各都市間の距離を計算しておく
  edgeDistances = vector<vector<double> > (size, vector<double> (size, 0));

  double edgeDistance;
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {

      // 距離を計算
      edgeDistance = sqrt(pow((coords[i].first - coords[j].first), 2) + 
          pow((coords[i].second - coords[j].second), 2));
      edgeDistances[i][j] = edgeDistance;
      edgeDistances[j][i] = edgeDistance;
    }
  }

}

// 経路の総距離を取得
double LKSolver::getCurrentTourDistance() {
  int currentIndex = 0;
  double distance = 0;
  for (int i = 0; i < size; i++) {
    //cout << edgeDistances[i][tour[i]] << "; ";
    distance += edgeDistances[i][tour[i]];
  }
  //cout << endl;
  return distance;
}

void LKSolver::LKMove(int tourStart) {
  // 変数の初期化
  set<pair<int,int> > broken_set, joined_set; // 壊されたエッジと結合されたエッジのペアを追跡するためのセット
  vector<int> tour_opt = tour; //最適化されたツアーを一時的に保存するためのベクター
  double g_opt = 0; //最適なツアーの長さ
  double g = 0; // 現在のツアーの長さ
  double g_local; // 局所的なツアーの長さ
  int lastNextV = tourStart; // 前回つなげた都市
  int fromV;// 現在の都市
  int nextV;// 次に訪れる都市
  int nextFromV;
  int lastPossibleNextV;
  pair<int, int> broken_edge;
  double y_opt_length;
  double broken_edge_length;
  double g_opt_local;
  // 現在の都市をtourStartの次として設定
  fromV = tour[lastNextV];
  long initialTourDistance = getCurrentTourDistance();

  do {
    nextV = -1;

    // 現在のエッジを一つ前から切り離す
    broken_edge = make_sorted_pair(lastNextV, fromV);
    broken_edge_length = edgeDistances[broken_edge.first][broken_edge.second];
    // 切り離すエッジがすでに結合セットに含まれている場合、ループを終了
    if (joined_set.count(broken_edge) > 0) break;

    // 現在の都市からつながっている都市を順番に見ていく

    for (int possibleNextV = tour[fromV]; nextV == -1 && possibleNextV != tourStart; possibleNextV = tour[possibleNextV]) {
      //繋ぎ変えた結果得られる利得を計算(大きいほどよい)
      g_local = broken_edge_length - edgeDistances[fromV][possibleNextV];
      // この都市をつなげるべきかどうかを判断
      if (!(
        broken_set.count(make_sorted_pair(fromV, possibleNextV)) == 0 &&
        g + g_local > 0 &&
        joined_set.count(make_sorted_pair(lastPossibleNextV, possibleNextV)) == 0 &&
        tour[possibleNextV] != 0 && 
        possibleNextV != tour[fromV] 
      )) {
        // この都市をつなげるべきではない場合
        lastPossibleNextV = possibleNextV;
        continue;
      }
      // この都市をつなげるべき場合
      nextV = possibleNextV;
    }
    
    if (nextV != -1) {
      // 現在のエッジを一つ前から切り離す
      broken_set.insert(broken_edge);
      // 現在の都市と次の都市のエッジを結合する
      joined_set.insert(make_sorted_pair(fromV, nextV));
      y_opt_length = edgeDistances[fromV][tourStart];
      g_opt_local = g + (broken_edge_length - y_opt_length);
      // エッジの交換が現在の最適解を改善するかどうか
      if (g_opt_local > g_opt) {
        g_opt = g_opt_local;
        tour_opt = tour;
        tour_opt[tourStart] = fromV;
      }
      g += broken_edge_length - edgeDistances[fromV][nextV];
      // 現在の都市から次の都市の手前までのツアーを反転
      reverse(fromV, lastPossibleNextV);
      nextFromV = lastPossibleNextV;
      tour[fromV] = nextV;
      lastNextV = nextV;
      fromV = nextFromV;

    }

  } while (nextV != -1);

  tour = tour_opt;
  long distanceAfter = getCurrentTourDistance();
  assert(distanceAfter <= initialTourDistance);

  assert(isTour());

}

void LKSolver::optimizeTour() {
  // we need to test for convergence and also the difference from last time
  int diff;
  double old_distance = getCurrentTourDistance();
  double new_distance = 0;
  for (int j = 0; j < 100; j++) {
    for (int i = 0; i < size; i++) {
      LKMove(i);
    }
    new_distance = getCurrentTourDistance();
    diff = old_distance - new_distance;
    if (j != 0) {
      assert(diff >= 0);
      if (diff == 0) {
        cout << j << "イテレーションで収束しました" << endl;
        break;
      }
    };
    old_distance = new_distance;
  }
}

// 経路の指定区間の反転
void LKSolver::reverse(int start, int end) {
  int current = start;
  int next = tour[start];
  int nextNext;
  do {
    //cout << "reversing" << endl;
    // この反復後にどこに行く必要があるかを先読みする
    nextNext = tour[next];
    tour[next] = current;
    current = next;
    next = nextNext;
  } while (current != end); // 終了条件は、最後に反転した都市が終了都市になるまで
}

// 巡回しているかどうかを確認
bool LKSolver::isTour() {
  int count = 1;
  int start = tour[0];
  while (start != 0) {
    start = tour[start];
    count++;
  }
  return (count == size);
}

void LKSolver::printTour() {
  int current = 0;
  do {
    cout << current << " ; ";
    current = tour[current];
  } while (current != 0);
  cout << endl;
}

void LKSolver::printTourIds(bool showCoords) {
  int current = 0;
  do {
    cout << ids[current];
    if (showCoords) {
      cout << "," << coords[current].first << "," << coords[current].second;
    }
    cout << endl;
    current = tour[current];
  } while (current != 0);
}

// Main

vector<int> id;
vector<pair<double,double> > coord;

int main(){
  int n;
  double x, y;
  while(cin >> n){
    id.push_back(n);
    cin >> x >> y;
    coord.push_back(make_pair(x, y));
  }

  LKSolver mat(coord, id);

  mat.optimizeTour();
  cout << mat.getCurrentTourDistance() << endl;
  mat.printTourIds(true);
}
