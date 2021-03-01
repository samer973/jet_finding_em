#include "hexbin.h"
#include "jet_header.h"

/* Software to emulate the hardware 2-layer jet-finding algorithm (integers). *Layers 1 and 2*           
 *
 * 2019
 *
 * Authors: Bennett Greenberg, Felipe Paucar-Velasquez, Yuri Gershtein, Sam Leigh
 * Rutgers, the State University of New Jersey
 *  Revolutionary for 250 years
 */

//Holds data from tracks, converted from their integer versions.
int nzbins = 6;
int main(int argc, char ** argv){
	nzbins = 6;
	int eventstart = 0;            
	int eventend = 999999; //999999;         
	if(argc == 2){
		nzbins = atoi(argv[1]);                                                 
	}
	else if(argc == 3){
		eventstart = atoi(argv[1]);
		eventend = atoi(argv[2]);
	}
	else if(argc == 4){
		nzbins = atoi(argv[1]);
		eventstart = atoi(argv[2]);
		eventend = atoi(argv[3]);
	}
	else if(argc > 4){
		cout << "Unable to understand arguments. Please include only nzbins, or eventstart eventend, or nzbins eventstart eventend." << endl;
		exit(0);
	}
	cout << "Running with " << nzbins << " zbins. " << endl;
	string fname = "phi";
     //number of bits the firmware will use for pT
	const int pTbits = 9;
     //the three parts of the incoming track
	int pTinverse; //as of now, it is actually just pT
	int eta;
	int pT;
	string data_in;
	string bin_data;
	string filename;
	//Open input file, read data to tracks array
	//number of tracks in each phi bin is 24 for now.
	ifstream in_tracks[nphibins*2/3];
	int ntrks[nphibins*2/3];
	for(int i = 0; i < nphibins/3; ++i){
		string s;
		stringstream out;
		out << i;
		s = out.str();
		filename = fname + s + "_p.dat";
		in_tracks[i*2+1].open(filename.c_str());
		ntrks[i*2+1] = 0;
	}
	for(int i = 0; i < nphibins/3; ++i){
		string s;
		stringstream out;
		out << i;
		s = out.str();
		filename = fname + s + "_n.dat";
		in_tracks[i*2].open(filename.c_str());
		ntrks[i*2] = 0;
	}
	ofstream out_clusts;
	string outname = "int_em_out.txt";
	out_clusts.open(outname.c_str());
	int zint;
	int etaint;
	int pTint;
	int phi;
	int eta_start;
	int eta_reg;
	int zb1;
	int zb2;
	int rel_phi_bin;
	int phi_bin;
	bool z_big;
	bool z_vbig;
	string data;
	string bin_z;
	string bin_eta;
	string bin_pT;
	string bin_phi;
	string bin_nt;
	string bin_nx;
	int ntracks;
	int maxntrk = 0;
	int last_maxntrk;
	int nevents = 0;
	int nreal_events = 0;
	for(nevents = 0; nevents <= eventend; ++nevents){
	    ntracks = 0;
	    last_maxntrk = maxntrk;
	    maxntrk = 0;
	    vector<track_data> tracks(0);
         //read from one phi slice at a time
	    for(int pslice = 0; pslice < nphibins*2/3; ++pslice){
		getline(in_tracks[pslice], data_in);
		if(data_in == "") {
			    goto data_read;
		}
	//Get all the lines of 0s of the last event (consequence of all phisectors having same number of tracks for each event).
		for(int tk=ntrks[pslice]; tk < last_maxntrk; ++tk){
			getline(in_tracks[pslice], data_in);
		}
		ntrks[pslice] = 0;
		while(true) {
			//if data is 0, event has ended. 
			if(data_in == "0x000000000000000000000000"){                        
				if(nevents < eventstart){
					cout << "ok";
					ntracks = 0;
					break;   
				}
				break;                                                       
			}//end data is 0; event has ended
			if(data_in == "") { //if end of file reached, we're all done reading in data.
				    goto data_read;
			}
			ntrks[pslice]++;
			bin_data = hex_to_bin(data_in, 96); 
			pTinverse = bin_to_int(bin_data.substr(82, 14));
			eta = bin_to_int(bin_data.substr(53, 16));
			track_data trkd;
			if(bin_data.substr(0, 1) == "1"){
				trkd.xbit = true;
			}
			else {
				trkd.xbit = false;
			}
			pT = pTinverse;
			if(pT > 511) pT = 511;
			trkd.pT = pT;
			cout << "eta str: " + bin_data.substr(53, 16) << endl;
			cout << "eta(0): " + bin_data.substr(53, 1) << endl;
			cout << "eta num: ";
			cout << eta;
			cout << endl;
			if(bin_to_int(bin_data.substr(53, 1)) == 1) {
				eta_start = bin_to_int(bin_data.substr(54, 2))*3;
			}
			else if(bin_to_int(bin_data.substr(53, 1)) == 0) {
				eta_start = (4+bin_to_int(bin_data.substr(54, 2)))*3;
			}
			if(bin_to_int(bin_data.substr(56, 13)) < 2730) {
				eta_reg = 0;
			}
			else if(bin_to_int(bin_data.substr(56, 13)) < 5460) {
				eta_reg = 1;			
			}
			else {
				eta_reg = 2;
			}
			cout << "eta bin: ";
			cout << (eta_start+eta_reg);
			cout << endl;
			trkd.eta = (eta_start+eta_reg); //-1.0 * log(sqrt(tf*tf + 1.0) - tf);
			if(bin_to_int(bin_data.substr(42, 11)) > 683) {
				z_big = true;
			}
			else {
				z_big = false;
			}
			if(bin_to_int(bin_data.substr(42, 11)) > 1365) {
				z_vbig = true;
			}
			else {
				z_vbig = false;
			}
			if(bin_to_int(bin_data.substr(41, 1)) == 1) {
				if(z_vbig) {
					zb1 = 2;
				}
				else {
					zb1 = 0;
				}
				if(z_big) {
					zb2 = 15;
				}
				else {
					zb2 = 1;
				}
			}
			else {
				if(z_big) {
					zb1 = 4;
				}
				else {
					zb1 = 2;
				}
				if(z_vbig) {
					zb2 = 15;
				}
				else {
					zb2 = 3;
				}
			}
			trkd.z1 = zb1;
			trkd.z2 = zb2;
			if(bin_to_int(bin_data.substr(69, 1)) == 1) {
				if(2047-bin_to_int(bin_data.substr(70, 11)) > 682) {
					rel_phi_bin = 0;
				}
				else {
					rel_phi_bin = 1;
				}
			}
			else {
				if(bin_to_int(bin_data.substr(70, 11)) < 682) {
					rel_phi_bin = 1;
				}
				else {
					rel_phi_bin = 2;
				}
			}
			if(pslice % 2 == 0) {
				phi_bin = 3*(pslice/2) + rel_phi_bin;
			}
			else {
				phi_bin = 3*(pslice-1)/2 + rel_phi_bin;
			}
			cout << "zb1: ";
			cout << zb1;
			cout << endl;
			cout << "zb2: ";
			cout << zb2;
			cout << endl;
			cout << "phi bin: ";
			cout << phi_bin;
			cout << endl;
			trkd.phi = phi_bin; 
			trkd.bincount = 0;
			++ntracks;                                       
			tracks.push_back(trkd);
			getline(in_tracks[pslice], data_in);
		}
		if(maxntrk < ntrks[pslice]) maxntrk = ntrks[pslice];
	}
			
	//Call L2cluster
	//
	//Then print to output file all clusters from the most energetic zbin 
data_read:
		if(ntracks == 0){ continue;}
		cout << "****EVENT " << nreal_events << " ****" << endl;
		nreal_events++;
		maxzbin mzb = L2_cluster(tracks, nzbins, ntracks); //left off here
		if(mzb.isEmpty == true) { 
			continue;
		}
	for(int kk = 0; kk < nzbins-1; ++kk){	
	        if(kk != mzb.znum) continue;
		vector<etaphibin> clusters = all_zbins[kk].clusters;
		for(int k = 0; k < all_zbins[kk].nclust; ++k){
			//Convert z, eta, and pT to their integer equivalents.
			//Measure from the middle of each bin to avoid errors.
			//zbins will go from 0 to 4.
			zint = all_zbins[kk].znum;
			etaint = clusters[k].eta;
			phi = clusters[k].phi;
			if(clusters[k].pTtot < pow(2, pTbits)){
				pTint = clusters[k].pTtot;
			}
			else {
				pTint = (1 << pTbits) - 1;
			}
			int nt = clusters[k].numtracks;
			int nx = clusters[k].nx_tracks;
			//Convert integers to binary strings.
			bin_z = int_to_bin(zint, 4);
			bin_eta = int_to_bin(etaint, 5);
			bin_pT = int_to_bin(pTint, pTbits);
			bin_phi = int_to_bin(phi, 5);
			bin_nt = int_to_bin(nt, 5);
			bin_nx = int_to_bin(nx, 4);
			//Concatenate into one binary string, convert to hex, write to output file.
			data = bin_nt + bin_nx + bin_z + bin_eta + bin_phi + bin_pT;
			data = bin_to_hex(data);
			if(bin_nt != "00000" && (pTint < 50 || (bin_nt != "00001" && pTint < 100) || (bin_nt != "00001" && bin_nt != "00010"))) {
				out_clusts << data << endl;	
			}
		}
         } //for each zbin
		if(mzb.ht == 0) cout << "WARNING: HT = 0 (Event " << nevents << ")" << endl;
		out_clusts << "00000000" << endl;
    }
	for(int it = 0; it < nphibins*2/3; ++it) {
        	in_tracks[it].close();
	}
	out_clusts.close();
	cout << "Data written to " << outname << endl;
	return 0;
}

vector<maxzbin> all_zbins(nzbins);
//input array of track_data, output zbin of maximum ht.
maxzbin L2_cluster(vector<track_data> tracks, int nzbins, int ntracks){
	maxzbin mzb = all_zbins[0];
    //returns NULL if there are no tracks for this event.
        if(ntracks == 0){
	      mzb.isEmpty = true;
	      return mzb;
	}
	else {
	     mzb.isEmpty = false;
	}
	//Create grid of phibins! 
	etaphibin epbins[nphibins][netabins];
	for(int i = 0; i < nphibins; ++i){
            for(int j = 0; j < netabins; ++j){
		epbins[i][j].phi = i;
		epbins[i][j].eta = j;
	     }//for each phibin
	 } //for each etabin (finished creating epbins)

	 //Last zbin won't be used (goes beyond maximum z)
	for(int zbin = 0; zbin < nzbins-1; ++zbin){
	
	      //First initialize pT, numtracks, nx_tracks, used to 0 (or false)
	        for(int i = 0; i < nphibins; ++i){
			for(int j = 0; j < netabins; ++j){
				epbins[i][j].pTtot = 0;
				epbins[i][j].used = false;
				epbins[i][j].numtracks = 0;
				epbins[i][j].nx_tracks = 0;
			}//for each phibin
		} //for each phibin

	      //Fill in etaphibins grid with pT from each track.
		for(int k = 0; k < ntracks; ++k) {
			for(int i = 0; i < nphibins; ++i){
				for(int j = 0; j < netabins; ++j){
					if((tracks[k].z1 == zbin || tracks[k].z2 == zbin) &&
					   (tracks[k].eta == epbins[i][j].eta) && 
					   (tracks[k].phi == epbins[i][j].phi) &&
					   (tracks[k].pT > 0 && tracks[k].bincount != 2)){
						++tracks[k].bincount;
						if((epbins[i][j].pTtot + tracks[k].pT) > 511) {
							epbins[i][j].pTtot = 511;
						}
						else {
							epbins[i][j].pTtot += tracks[k].pT;
						}
						++epbins[i][j].numtracks;
						if(tracks[k].xbit) {
							++epbins[i][j].nx_tracks;
						}
						cout << "track " << k << " zbin " << zbin << " phibin " << i << " etabin " << j << " pT " << tracks[k].pT << " xbit " << tracks[k].xbit << endl; 
					} //if right bin
				} //for each phibin: j loop
			}//for each phibin: i loop
		} //for each track: k loop

    //Uncomment to print out pT of each eta and phi bin.
		for(int i = 0; i < nphibins; ++i)
			for(int j = 0; j < netabins; ++j)
				if(epbins[i][j].pTtot != 0) {
//				        cout << "zmin " << zmin << " zmax " << zmax << endl;
					cout << "zbin " << zbin << " epbins[" << i << "][" << j << "] pTtot: " << epbins[i][j].pTtot << endl;
                                      //  cout << epbins[i][j].phi << "\t" << epbins[i][j].pTtot << endl;					
				}
	

	  //First do clustering in Layer 1: maximum possible nclust for each eta slice would be a cluster in every other phibin.
		vector<etaphibin> L1clusters[nphibins];
                for(int phislice = 0; phislice < nphibins; ++phislice){
			L1clusters[phislice] = L1_cluster(epbins[phislice]);
			for(int ind = 0; L1clusters[phislice][ind].pTtot != 0; ++ind){
				L1clusters[phislice][ind].used = false;
			}
		}

	//Create clusters array to hold output cluster data for Layer2; can't have more clusters than tracks.
		vector<etaphibin> L2cluster;// = (etaphibin *)malloc(ntracks * sizeof(etaphibin));

	//Find eta-phibin with maxpT, make center of cluster, add neighbors if not already used.
		int hipT = 0;
		int nclust = 0;
		int phibin = 0;
		int imax;
	     //index of clusters array for each phislice.
		int index1;
		int E1 =0;
		int E0 =0;
		int E2 =0;
		int trx1, trx2;
		int xct1, xct2;
		int used1, used2, used3, used4;

			//Find eta-phibin with highest pT.
		for(phibin = 0; phibin < nphibins; ++phibin){
		    while(true){
			hipT = 0;
			for(index1 = 0; L1clusters[phibin][index1].pTtot > 0; ++index1){
				if(!L1clusters[phibin][index1].used && L1clusters[phibin][index1].pTtot >= hipT){
					hipT = L1clusters[phibin][index1].pTtot;
					imax = index1;
				}
			}//for each index within the phibin
		      //If highest pT is 0, all bins are used.
			if(hipT == 0){
				break;
			}
			E0 = hipT;   //E0 is pT of first phibin of the cluster.
			E1 = 0;
			E2 = 0;
			trx1 = 0;
			trx2 = 0;
			xct1 = 0;
			xct2 = 0;
			L2cluster.push_back(L1clusters[phibin][imax]);
			L1clusters[phibin][imax].used = true;
		//Add pT of upper neighbor.
		//E1 is pT of the middle phibin (should be highest pT)
			if(phibin != nphibins-1){
				used1 = -1;
				used2 = -1;
				for (index1 = 0; L1clusters[phibin+1][index1].pTtot != 0; ++index1){
					if(L1clusters[phibin+1][index1].used){
						continue;
					}
					if(abs(L1clusters[phibin+1][index1].eta - L1clusters[phibin][imax].eta) <= 1){
						E1 += L1clusters[phibin+1][index1].pTtot;
						trx1 += L1clusters[phibin+1][index1].numtracks;
						xct1 += L1clusters[phibin+1][index1].nx_tracks;
						if(used1 < 0)
							used1 = index1;
						else
							used2 = index1;
					}//if cluster is within one phibin
				} //for each cluster in above phibin
			//if E1 isn't higher, E0 and E1 are their own cluster.
				if(E1 < E0){
					L2cluster[nclust].pTtot += E1;   
					L2cluster[nclust].numtracks += trx1;
					L2cluster[nclust].nx_tracks += xct1;
					if(used1 >= 0)
						L1clusters[phibin+1][used1].used = true;
					if(used2 >= 0)
						L1clusters[phibin+1][used2].used = true;
					++nclust;
					continue;
				}
				
				if(phibin != nphibins-2){
                                      //E2 will be the pT of the third phibin (should be lower than E1).
					used3 = -1;
					used4 = -1;
					for (index1 = 0; L1clusters[phibin+2][index1].pTtot != 0; ++index1){
						if(L1clusters[phibin+2][index1].used){
							continue;
						}
						if(abs(L1clusters[phibin+2][index1].eta - L1clusters[phibin][imax].eta) <= 1){
							E2 += L1clusters[phibin+2][index1].pTtot;
							trx2 += L1clusters[phibin+2][index1].numtracks;
							xct2 += L1clusters[phibin+2][index1].nx_tracks;
							if(used3 < 0)
								used3 = index1;
							else
								used4 = index1;
						}
		
					}
				     //if indeed E2 < E1, add E1 and E2 to E0, they're all a cluster together.
				     //  otherwise, E0 is its own cluster.
					if(E2 < E1){
						L2cluster[nclust].pTtot += E1 + E2;
						L2cluster[nclust].numtracks += trx1 + trx2;
						L2cluster[nclust].nx_tracks += xct1 + xct2;
						L2cluster[nclust].phi = L1clusters[phibin+1][used1].phi;	
						if(used1 >= 0)
							L1clusters[phibin+1][used1].used = true;
						if(used2 >= 0)
							L1clusters[phibin+1][used2].used = true;
						if(used3 >= 0)
							L1clusters[phibin+2][used3].used = true;
						if(used4 >= 0)
							L1clusters[phibin+2][used4].used = true;
					}
					++nclust;
					continue;
				} // end Not nphibins-2
				else{
					L2cluster[nclust].pTtot += E1;
					L2cluster[nclust].numtracks += trx1;
					L2cluster[nclust].nx_tracks += xct1;
					L2cluster[nclust].phi = L1clusters[phibin+1][used1].phi;
					if(used1 >= 0)
						L1clusters[phibin+1][used1].used = true;
					if(used2 >= 0)
						L1clusters[phibin+1][used2].used = true;
					++nclust;
					continue;
				}
			}//End not last phibin(26)
			else { //if it is phibin 26
				L1clusters[phibin][imax].used = true;
				++nclust;
			}
		    }//while hipT not 0
		}//for each phibin
		//for(int db=0;db<nclust;++db)cout<<L2cluster[db].phi<<"\t"<<L2cluster[db].pTtot<<"\t"<<L2cluster[db].numtracks<<endl;	
	//Now merge clusters, if necessary
		for(int m = 0; m < nclust -1; ++m){
                     for(int n = m+1; n < nclust; ++n)
                        if(L2cluster[n].eta == L2cluster[m].eta && (abs(L2cluster[n].phi - L2cluster[m].phi) <= 1 || abs(L2cluster[n].phi - L2cluster[m].phi) >= 26)){
                                if(L2cluster[n].pTtot > L2cluster[m].pTtot){
                                        L2cluster[m].phi = L2cluster[n].phi;
                                }
                                L2cluster[m].pTtot += L2cluster[n].pTtot;
                                L2cluster[m].numtracks += L2cluster[n].numtracks;
				L2cluster[m].nx_tracks += L2cluster[n].nx_tracks;
                                for(int m1 = n; m1 < nclust-1; ++m1){
                                        L2cluster[m1] = L2cluster[m1+1];
                                }
                                nclust--;
                                m = -1;
                                break; 
                        }//end if clusters neighbor in eta
                }//end for (m) loop     
//		for(int db=0;db<nclust;++db)cout<<L2cluster[db].phi<<"\t"<<L2cluster[db].pTtot<<"\t"<<L2cluster[db].numtracks<<endl;	
          //sum up all pTs in this zbin to find ht.
		int ht = 0;
		for(int k = 0; k < nclust; ++k){
			//(pTint < 50 || (bin_nt != "00001" && pTint < 100) || (bin_nt != "00001" && bin_nt != "00010")
			if(L2cluster[k].pTtot < 50 || (L2cluster[k].pTtot < 100 && L2cluster[k].numtracks >= 2) || L2cluster[k].numtracks >= 3) {
				ht += L2cluster[k].pTtot;
			}
			if(L2cluster[k].pTtot > 511) L2cluster[k].pTtot = 511;
                }
		if(ht > 511) ht = 511;

	   //if ht is larger than previous max, this is the new vertex zbin.
		all_zbins[zbin].znum = zbin;
		all_zbins[zbin].nclust = nclust;
		all_zbins[zbin].clusters.clear();
		etaphibin allzclust;
		for(int k = 0; k < nclust; ++k){
			allzclust.phi = L2cluster[k].phi;                               
			allzclust.eta = L2cluster[k].eta;                             
			allzclust.pTtot = L2cluster[k].pTtot;
			allzclust.numtracks = L2cluster[k].numtracks;
			if(allzclust.numtracks > MAXNT) allzclust.numtracks = MAXNT;
			allzclust.nx_tracks = L2cluster[k].nx_tracks;
			if(allzclust.nx_tracks > MAXNX) allzclust.nx_tracks = MAXNX;
			all_zbins[zbin].clusters.push_back(allzclust);
		}
	//	for(int db=0;db<nclust;++db)cout<<all_zbins[zbin].clusters[db].phi<<"\t"<<all_zbins[zbin].clusters[db].pTtot<<endl;	
		all_zbins[zbin].ht = ht;
		if(ht >= mzb.ht){
			mzb = all_zbins[zbin];
		}
	       //Prepare for next zbin!
	     } //for each zbin
       return mzb;
}

vector<etaphibin> L1_cluster(etaphibin phislice[netabins]){

		vector<etaphibin> clusters(netabins/2);
	//Find eta-phibin with maxpT, make center of cluster, add neighbors if not already used.
		int my_pt, left_pt, right_pt, right2pt;
		int nclust = 0;

		for(int etabin = 0; etabin < netabins; ++etabin){
			//assign values for my pT and neighbors' pT
			if(phislice[etabin].used) continue;
			if(phislice[etabin].pTtot > 511) {
				my_pt = 511;
			}
			else {
				my_pt = phislice[etabin].pTtot;
			}
			//cout << my_pt << endl;
			if(etabin > 0 && !phislice[etabin-1].used) {
				left_pt = phislice[etabin-1].pTtot;
			} 
			else {
				left_pt = 0;
			}
			if(etabin < netabins - 1 && !phislice[etabin+1].used) {
				right_pt = phislice[etabin+1].pTtot;
				if(etabin < netabins - 2 && !phislice[etabin+2].used) {
					right2pt = phislice[etabin+2].pTtot;
				} 
				else {
					right2pt = 0;
				}
			} 
			else {
				right_pt = 0;
			}
		
		//if I'm not a cluster, move on.
			if(my_pt < left_pt || my_pt <= right_pt) {
			   //if unused pT in the left neighbor, spit it out as a cluster.
			        if(left_pt > 0) {
					cout << "my pt " << my_pt << " left_pt " << left_pt << " right_pt " << right_pt; 
					clusters[nclust] = phislice[etabin-1];
					//clusters.push_back(phislice[etabin-1]);
					phislice[etabin-1].used = true;
					++nclust;
				}
				continue;
			}

		//I guess I'm a cluster-- should I use my right neighbor?
		// Note: left neighbor will definitely be used because if it 
		//       didn't belong to me it would have been used already
//			clusters.push_back(phislice[etabin]);
			clusters[nclust] = phislice[etabin];
			phislice[etabin].used = true;
			if(left_pt > 0) {
				if((clusters[nclust].pTtot + left_pt) > 511) {
					clusters[nclust].pTtot = 511;
				}
				else {
					clusters[nclust].pTtot += left_pt;
					//cout << clusters[nclust].pTtot << endl;
					//cout << etabin << endl;
				}
				
				clusters[nclust].numtracks += phislice[etabin-1].numtracks;
				clusters[nclust].nx_tracks += phislice[etabin-1].nx_tracks;
			}
			if(my_pt >= right2pt && right_pt > 0) {
				if((clusters[nclust].pTtot + right_pt) > 511) {
					clusters[nclust].pTtot = 511;

				}
				else {
					clusters[nclust].pTtot += right_pt;
					//cout << clusters[nclust].pTtot << endl;
					//cout << etabin << endl;
				}
				//cout << clusters[nclust].pTtot << endl;
				clusters[nclust].numtracks += phislice[etabin+1].numtracks;
				clusters[nclust].nx_tracks += phislice[etabin+1].nx_tracks;
				phislice[etabin+1].used = true;
			}

			++nclust;
		} //for each etabin                       
	                         
	//Now merge clusters, if necessary
		for(int m = 0; m < nclust -1; ++m){
			if(abs(clusters[m+1].eta - clusters[m].eta) <= 1){
				if(clusters[m+1].pTtot > clusters[m].pTtot){
					clusters[m].eta = clusters[m+1].eta;
				}
				if((clusters[m].pTtot + clusters[m+1].pTtot) > 511) {
					clusters[m].pTtot = 511;
				}
				else {
					clusters[m].pTtot += clusters[m+1].pTtot;
				}
				//cout << clusters[m].pTtot << endl;
				clusters[m].numtracks += clusters[m+1].numtracks;
				clusters[m].nx_tracks += clusters[m+1].nx_tracks;
				for(int m1 = m+1; m1 < nclust-1; ++m1){
					clusters[m1] = clusters[m1+1];
				}
				nclust--;
				m = -1;
			}//end if clusters neighbor in eta
		}//end for (m) loop

	//zero out remaining unused clusters.
	for(unsigned int i = nclust; i < clusters.size(); ++i){
		clusters[i].pTtot = 0;
	}
	return clusters;
}
