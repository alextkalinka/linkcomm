/* Function for calculating link densities and the partition density for clusters generated within "getLinkCommunities".
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

void getLinkDensities(int *ma, int *mb, int *ea, int *eb, unsigned int *numedg, int *clusnums, int *numcl, double *pdens, double *heights, double *pdmax, int *csize, bool *removetrivial, bool *bipartite, int *bip, bool *verbose)

	{

	unsigned int i, j, k, ne, one = 0, csum = clusnums[0];
	int p = 0, nn, count = 0, rm = 0, n1, n2, which1, which2, nn0, nn1;
	float prog;
	double ldens, maxp, best = 0.0, denom;
	vector<int> mergeA;
	vector<int> mergeB;
	vector<int> edgeA;
	vector<int> edgeB;
	vector<int> drop;
	vector<int> dropTemp;
	vector<int> rem;
	map<int, vector<int> > clusters;
	vector<int> current;
	vector<int> currentTemp;
	vector<int> bestC;
	vector<int> cdel; // Indices of clusters that can be deleted to free up memory.
	vector<int> cdelTemp;
	vector<int>::iterator bestIt;
	set<int> nodes;
	set<int> nodes0; // For bipartite networks.
	set<int> nodes1;

	copy(ma, ma + *numedg-1, back_inserter(mergeA));
	copy(mb, mb + *numedg-1, back_inserter(mergeB));

	copy(ea, ea + *numedg, back_inserter(edgeA));
	copy(eb, eb + *numedg, back_inserter(edgeB));

	remove("linkcomm_clusters.txt");

	ofstream outfile;
	outfile.open("linkcomm_clusters.txt", ios::out);
	if(! outfile.is_open()){
			Rprintf("\nERROR: can't open linkcomm_clusters.txt!\n"); return;
			}


	// Loop through merges from lowest to highest.
	// At each merge extract cluster elements.
	for(i = 0; i < *numedg-1; i++){
		
		if(*verbose){
			prog = (i+0.0)/(*numedg-2)*100;

			Rprintf("\r   Calculating link densities... %3.2f%%",prog);

			R_FlushConsole();
			R_ProcessEvents();
			}

		if(mergeA.at(i) < 0 && mergeB.at(i) < 0){
			clusters[i].push_back(-1*mergeA.at(i));
			clusters[i].push_back(-1*mergeB.at(i));
			sort(clusters[i].begin(),clusters[i].end());
			current.push_back(i); // Update current cluster subs.
			one = 0;
		}else if(mergeA.at(i) > 0 && mergeB.at(i) < 0){
			clusters[i] = clusters[mergeA.at(i)-1];
			clusters[i].push_back(-1*mergeB.at(i));
			cdel.push_back(mergeA.at(i)-1); // Index of merged cluster so it can be deleted when possible.
			sort(clusters[i].begin(),clusters[i].end());
			current.push_back(i);
			drop.push_back(mergeA.at(i) - 1);
			one = 1;
		}else if(mergeA.at(i) < 0 && mergeB.at(i) > 0){
			clusters[i] = clusters[mergeB.at(i)-1];
			clusters[i].push_back(-1*mergeA.at(i));
			cdel.push_back(mergeB.at(i)-1);
			sort(clusters[i].begin(),clusters[i].end());
			current.push_back(i);
			drop.push_back(mergeB.at(i) - 1);
			one = 1;
		}else{
			clusters[i].insert(clusters[i].end(), clusters[mergeA.at(i)-1].begin(), clusters[mergeA.at(i)-1].end());
			clusters[i].insert(clusters[i].end(), clusters[mergeB.at(i)-1].begin(), clusters[mergeB.at(i)-1].end());
			cdel.push_back(mergeA.at(i)-1);
			cdel.push_back(mergeB.at(i)-1);
			sort(clusters[i].begin(),clusters[i].end());
			current.push_back(i);
			drop.push_back(mergeA.at(i) - 1);
			drop.push_back(mergeB.at(i) - 1);
			one = 2;
			}

		// Clear out memory.
		if(one != 0){
			cdelTemp = cdel;
			rm = 0;
			for(j = 1; j <= one; j++){
				bestIt = find(bestC.begin(), bestC.end(), cdel.at( cdel.size()-j ));
				if(bestIt == bestC.end()){
					clusters.erase( cdel.at( cdel.size()-j ) );
					cdelTemp.erase( cdelTemp.end() - j + rm );
					rm++;
					}
				}
			cdel = cdelTemp;
			cdelTemp.clear();
			}


		if(i+1 == csum){ // Reached end of this height, calculate link densities for the current clusters.
			currentTemp = current;
			sort(drop.begin(),drop.end());
			for(j = 0; j < drop.size(); j++){
				bestIt = find(currentTemp.begin(), currentTemp.end(), drop.at(j));
				currentTemp.erase( bestIt );
				}
			drop.clear();

			ldens = 0;
			for(j = 0; j < currentTemp.size(); j++){
				// Number of edges.
				ne = clusters[currentTemp.at(j)].size();
				// Number of nodes.
				for(k = 0; k < ne; k++){
					n1 = edgeA.at(clusters[currentTemp.at(j)].at(k)-1);
					n2 = edgeB.at(clusters[currentTemp.at(j)].at(k)-1);
					nodes.insert(n1);
					nodes.insert(n2);
					if(*bipartite){
						which1 = bip[(n1-1)];
						which2 = bip[(n2-1)];
						if(which1 == 0){
							nodes0.insert(n1);
						}else{
							nodes1.insert(n1);
							}
						if(which2 == 0){
							nodes0.insert(n2);
						}else{
							nodes1.insert(n2);
							}
						}
					}
				nn = nodes.size();
				if(!(*bipartite)){
					ldens = ldens + (double(ne)*(double(ne)-double(nn)+1.0))/((double(nn)-2.0)*(double(nn)-1.0));
				}else{
					nn0 = nodes0.size();
					nn1 = nodes1.size();
					denom = (double(nn0)*double(nn1)*2.0)-((double(nn)-1.0)*2.0);
					if(!(denom==0)){
						ldens = ldens + (double(ne)*(double(ne)-double(nn)+1.0))/denom;
						}
					}
				nodes.clear();
				nodes0.clear();
				nodes1.clear();
				}

			pdens[p] = (2.0/ *numedg)*ldens;
			maxp = pdens[p];

			if(maxp > best){
				best = maxp;
				*pdmax = heights[p];
				bestC = currentTemp;
				// Clusters at this height are more optimal so we can free memory used by previously merged clusters.
				for(j = 0; j < cdel.size(); j++){
					clusters.erase( cdel.at(j) );
					}
				cdel.clear();
				}

			// Clusters at this height are not optimal so we can free memory from previously merged clusters
			// that do not belong to the currently optimal cluster set.
			cdelTemp = cdel;
			rm = 0;
			for(j = 0; j < cdel.size(); j++){
				bestIt = find(bestC.begin(), bestC.end(), cdel.at(j));
				if(bestIt == bestC.end()){
					clusters.erase( cdel.at(j) );
					cdelTemp.erase( cdelTemp.begin() + j - rm );
					rm++;
					}
				}
			cdel = cdelTemp;
			cdelTemp.clear();

			current = currentTemp;
			currentTemp.clear();

			if(p < (*numcl-1)){
				p++;
				csum = csum + clusnums[p];
				}
			
			}

		}

	// Write optimal clusters to disk
	for(i = 0; i < bestC.size(); i++){
		if(*removetrivial && clusters[bestC.at(i)].size() == 2){count++; continue;} // Delete trivial clusters of size 2.
		for(j = 0; j < clusters[bestC.at(i)].size(); j++){
			outfile << clusters[bestC.at(i)].at(j) << " ";
			}
		outfile << endl;
		}

	*csize = bestC.size() - count;
	
	outfile.close();

	}

	}







