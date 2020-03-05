/* Miscellaneous functions for "linkcomm".
 *
 * Author: Alex T. Kalinka (alex.t.kalinka@gmail.com)
 *
 * Will be loaded as a shared object into R (dyn.load('name.so')) and called from within R.
 * Based on the Link communities derived from the algorithm in:
 * Ahn et al. (2010). Link communities reveal multiscale complexity in networks. Nature 466:761-765.
 *
 */

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <iterator>
#include <algorithm>
#include <set>
#include <map>

#include <R.h>

extern "C" {

using namespace std;

void cutTreeAt(int *ma, int *mb, double *heights, double *cutat, int *csize, int *numM)

	{

	unsigned int i, j;
	float prog;
	vector<int> mergeA;
	vector<int> mergeB;
	vector<int> drop;
	vector< vector<int> > clusters;
	vector<int> current;

	copy(ma, ma + *numM, back_inserter(mergeA));
	copy(mb, mb + *numM, back_inserter(mergeB));

	remove("linkcomm_metaclusters.txt");

	ofstream outfile;
	outfile.open("linkcomm_metaclusters.txt", ios::out);
	if(! outfile.is_open()){
			Rprintf("\nERROR: can't open linkcomm_metaclusters.txt!\n"); return;
			}


	// Loop through merges from lowest to highest.
	// At each merge extract cluster elements until we reach the cut point.
	for(i = 0; i < mergeA.size(); i++){

		prog = (i+0.0)/(mergeA.size()-1)*100;

		Rprintf("\r   Extracting clusters... %3.2f%%",prog);

		R_FlushConsole();
		R_ProcessEvents();
		
		if(mergeA.at(i) < 0 && mergeB.at(i) < 0){
			clusters.push_back(vector<int>());
			clusters[i].push_back(-1*mergeA.at(i));
			clusters[i].push_back(-1*mergeB.at(i));
			sort(clusters[i].begin(),clusters[i].end());
			current.push_back(i); // Update current cluster subs.
		}else if(mergeA.at(i) > 0 && mergeB.at(i) < 0){
			clusters.push_back(vector<int>());
			clusters[i] = clusters[mergeA.at(i)-1];
			clusters[i].push_back(-1*mergeB.at(i));
			sort(clusters[i].begin(),clusters[i].end());
			current.push_back(i);
			drop.push_back(mergeA.at(i) - 1);
		}else if(mergeA.at(i) < 0 && mergeB.at(i) > 0){
			clusters.push_back(vector<int>());
			clusters[i] = clusters[mergeB.at(i)-1];
			clusters[i].push_back(-1*mergeA.at(i));
			sort(clusters[i].begin(),clusters[i].end());
			current.push_back(i);
			drop.push_back(mergeB.at(i) - 1);
		}else{
			clusters.push_back(vector<int>());
			clusters[i].insert(clusters[i].end(), clusters[mergeA.at(i)-1].begin(), clusters[mergeA.at(i)-1].end());
			clusters[i].insert(clusters[i].end(), clusters[mergeB.at(i)-1].begin(), clusters[mergeB.at(i)-1].end());
			sort(clusters[i].begin(),clusters[i].end());
			current.push_back(i);
			drop.push_back(mergeA.at(i) - 1);
			drop.push_back(mergeB.at(i) - 1);
			}


		}

	// Remove nested clusters.
	sort(drop.begin(),drop.end());
	for(i = 0; i < drop.size(); i++){
		current.erase(current.begin() + drop.at(i) - i);
		}

	// Write meta-clusters to disk.
	for(i = 0; i < current.size(); i++){
		for(j = 0; j < clusters[current.at(i)].size(); j++){
			outfile << clusters[current.at(i)].at(j) << " ";
			}
		outfile << endl;
		}
	
	*csize = current.size();

	outfile.close();


	}


void getJaccards(int *nodes, int *clusters, int *clusids, unsigned int *numNodes, double *dissvec, bool *verbose)
	
	{

	unsigned int i, j = 0;
	int runn = 0;
	float prog;
	set<int> tempNodes;
	map<int, set<int> > clusMap;
	set<int> common;
	set<int> total;
	set<int>::iterator sit;

	// Construct map of clusids and their nodes as integers.
	for(i = 0; i < *numNodes; i++){
		
		if(*verbose){
			prog = (i+0.0)/(*numNodes)*100;

			Rprintf("\r   Calculating cluster similarities 1/2... %3.2f%%",prog);

			R_FlushConsole();
			R_ProcessEvents();
			}

		if(clusters[i] == clusids[j]){
			tempNodes.insert(nodes[i]);
		}else{
			clusMap.insert( pair<int,set<int> >(clusids[j],tempNodes) );
			tempNodes.clear();
			tempNodes.insert(nodes[i]);
			j++;
			}
		if(i == *numNodes-1){
			clusMap.insert( pair<int,set<int> >(clusids[j],tempNodes) );
			tempNodes.clear();
			}
		}

	// Loop through clusters and calculate Jaccards.
	for(i = 0; i < clusMap.size()-1; i++){
		
		if(*verbose){
			prog = (i+0.0)/(clusMap.size()-2)*100;

			Rprintf("\r   Calculating cluster similarities 2/2... %3.2f%%",prog);

			R_FlushConsole();
			R_ProcessEvents();
			}

		for(j = i+1; j < clusMap.size(); j++){
			
			set_intersection(clusMap[clusids[i]].begin(),clusMap[clusids[i]].end(),clusMap[clusids[j]].begin(),clusMap[clusids[j]].end(), inserter(common, common.begin()));

			set_union(clusMap[clusids[i]].begin(),clusMap[clusids[i]].end(),clusMap[clusids[j]].begin(),clusMap[clusids[j]].end(), inserter(total, total.begin()));

			// Jaccard coefficient.
			dissvec[j - 1 + runn] = 1.0 - (common.size()+0.0)/(total.size()+0.0);

			common.clear();
			total.clear();

			}

		runn = runn + clusMap.size() - 2 - i;

		
		}


	}




void getNumClusters(int *unn, int *nodes, int *counts, unsigned int *numnodes, unsigned int *nrows, bool *verbose)
	{

	unsigned int i, j;
	int nd;
	float prog;

	// Loop through nodes and count their community membership.
	for(i = 0; i < *numnodes; i++){

		if(*verbose){
			prog = (i+0.0)/(*numnodes-1.0)*100;

			Rprintf("\r   Finishing up...4/4... %3.2f%%",prog);

			R_FlushConsole();
			R_ProcessEvents();
			}   

		nd = unn[i];

		for(j = 0; j < *nrows; j++){

			if(nodes[j] == nd){
				counts[i]++;
				}
			}

		}

	}







	}












