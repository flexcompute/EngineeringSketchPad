#include <string>
using namespace std;

#include "rcm.hpp"
#include "genrcmi.hpp"

//#define DEBUG
#ifdef DEBUG
#include <iostream>
#include <vector>
#endif

extern "C"
void genrcmi(int node_num, int adj_num, int *adj_row, int *adj, int *perm)
{
#ifdef DEBUG
    int bandwidth = adj_bandwidth(node_num, adj_num, adj_row, adj);
    cout << endl << "  ADJ bandwidth = " << bandwidth << endl;
#endif

  genrcm(node_num, adj_num, adj_row, adj, perm);

#ifdef DEBUG
  vector<int> perm_inv(node_num);
  perm_inverse3(node_num, perm, perm_inv.data());


  bandwidth = adj_perm_bandwidth(node_num, adj_num, adj_row, adj,
                                 perm, perm_inv.data());

  cout << "  Permuted ADJ bandwidth = " << bandwidth << endl;
#endif
}
