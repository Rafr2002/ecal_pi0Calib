#include "TChain.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TF1.h"
#include "TH2F.h"
#include "TF1.h"
#include "TLegend.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <string>
#include "TCanvas.h"
// for ifstream
#include <fstream>
#include "TDecompSVD.h"
//include random number generator
#include <TRandom3.h>



using namespace std;

const double pi0_mass_pdg = 0.1349766;  // PDG pi0 mass in GeV

    const Double_t z_calo = 6; // position of calorimeter from the target in m
    const Double_t z_target = 0.09;    // position of target
    const Double_t z_origin = 0.0;
    const Double_t vertex_z = z_target - z_origin;  // position of vertex, where the pi0 is created, right in the middle of the target

    TVector3 vertex(0,0,z_target);


// 
// run the script with "root -l ecal_pi0calib.C"

// Find the highest and second-highest energy clusters
void find_highest_energy(const Double_t ecal_e[], int nclus, int &imax1, int &imax2, double &e1, double &e2) {
    for (int icl = 0; icl < nclus; ++icl) {
        if (ecal_e[icl] > e1) {
            // If this cluster has the highest energy so far:
            e2 = e1; imax2 = imax1;               // The old highest becomes the second-highest, The old highest index becomes the second-highest index
            e1 = ecal_e[icl]; imax1 = icl;        // Update highest energy and highest index
        } else if (ecal_e[icl] > e2) {
            // If this cluster is not the highest, but is higher than the second-highest:
            e2 = ecal_e[icl]; imax2 = icl;        // Update second-highest energy and the second-highest index
        }
    }
}
// Best pair of clusters within the time window
void best_pair_cuts(
    const Double_t ecal_e[], const Double_t clus_nblk[],
    const Double_t ecal_x[], const Double_t ecal_y[], 
    const Double_t clus_a_time[], 
    int nclus, int &best_icl, int &best_jcl, double &best_dt) {
    for (int icl = 0; icl < nclus; ++icl) {
        for (int jcl = icl + 1; jcl < nclus; ++jcl) {
            if ((ecal_e[icl] + ecal_e[jcl]) < 0.5) continue;            // sum of energies of the two clusters larger than 0.5 GeV
            if (ecal_e[icl] < 0.2 || ecal_e[jcl] < 0.2) continue;       // each cluster energy larger than 0.2 GeV
            if (clus_nblk[icl] < 2 || clus_nblk[jcl] < 2) continue;     // each cluster has at least 2 blocks
            Double_t deltaR = sqrt(pow(ecal_x[icl] - ecal_x[jcl], 2) + pow(ecal_y[icl] - ecal_y[jcl], 2));
            if (deltaR < 0.3) continue;                                // distance between clusters less than 0.3 m
            //pass_deltar++;
            if (clus_a_time[icl] < 80 || clus_a_time[icl] > 120) continue;  // cluster time in the range 80-120 ns
            if (clus_a_time[jcl] < 80 || clus_a_time[jcl] > 120) continue;  // cluster time in the range 80-120 ns
            double dt = clus_a_time[icl] - clus_a_time[jcl];
            if (fabs(dt) < fabs(best_dt) && fabs(dt) <= 50) { // time difference between clusters
                // Update the best pair of clusters
                best_dt = dt;
                best_icl = icl;
                best_jcl = jcl;
            }
        }
    }
}

// Best pair out of N clusters within the time window
void best_pair_cuts_NCompare(
    const Double_t ecal_e[], const Double_t clus_nblk[],
    const Double_t ecal_x[], const Double_t ecal_y[], 
    const Double_t clus_a_time[], 
    int nclus, int &best_icl, int &best_jcl, double &best_dt) {
    for (int icl = 0; icl < nclus; ++icl) {
        for (int jcl = icl + 1; jcl < nclus; ++jcl) {
            if ((ecal_e[icl] + ecal_e[jcl]) < 0.5) continue;            // sum of energies of the two clusters larger than 0.5 GeV
            if (ecal_e[icl] < 0.2 || ecal_e[jcl] < 0.2) continue;       // each cluster energy larger than 0.2 GeV
            if (clus_nblk[icl] < 2 || clus_nblk[jcl] < 2) continue;     // each cluster has at least 2 blocks
            Double_t deltaR = sqrt(pow(ecal_x[icl] - ecal_x[jcl], 2) + pow(ecal_y[icl] - ecal_y[jcl], 2));
            if (deltaR < 0.3) continue;                                // distance between clusters less than 0.3 m
            //pass_deltar++;
            if (clus_a_time[icl] < 80 || clus_a_time[icl] > 120) continue;  // cluster time in the range 80-120 ns
            if (clus_a_time[jcl] < 80 || clus_a_time[jcl] > 120) continue;  // cluster time in the range 80-120 ns
            double dt = clus_a_time[icl] - clus_a_time[jcl];
            if (fabs(dt) < fabs(best_dt) && fabs(dt) <= 50) { // time difference between clusters
                // Update the best pair of clusters
                best_dt = dt;
                best_icl = icl;
                best_jcl = jcl;
            }
        }
    }
}

// best pair by chi2 with m_pi0 and dt
void best_pair_chi2(
    const Double_t ecal_e[], const Double_t clus_nblk[],
    const Double_t ecal_x[], const Double_t ecal_y[], 
    const Double_t clus_a_time[], 
    int nclus, int &best_icl, int &best_jcl, double &best_dt, double &best_mass){

    // ######################
    best_icl = best_jcl = -1; best_dt = 1e9;  best_mass = 0.0;
    


    std::vector<int> cand; //candidate clusters 
    // std::vector<int> cand(8, 0); //candidate clusters 

    for (int icl = 0; icl < nclus; ++icl) {
        if (ecal_e[icl] < 0.2) continue; // each cluster energy larger than 0.2 GeV
        if ((int)clus_nblk[icl] < 2) continue;     // each cluster has at least 2 blocks
        if (clus_a_time[icl] < 80 || clus_a_time[icl] > 120) continue; // cluster time in the range 80-120 ns
        cand.push_back(icl);
    }

    if (cand.size() < 2) return;

    // ---- keep only top-N by energies ----
    int Nmax=8;
    std::sort(cand.begin(), cand.end(),
    [&](int a, int b){ return ecal_e[a] > ecal_e[b]; });
    if ((int)cand.size() > Nmax) {cand.resize(Nmax);};

     // ---- loop pairs among candidates ----
    for (size_t i = 0; i < cand.size(); ++i){
        const int icl = cand[i];
        for (size_t j = i+1; j < cand.size(); ++j){
            const int jcl = cand[j];

           double E1 = ecal_e[icl], E2 = ecal_e[jcl];
            if (E1 + E2 < 0.6) continue;  // sum of energies of the two clusters less than 0.6 GeV

           Double_t deltaR = sqrt(pow(ecal_x[icl] - ecal_x[jcl], 2) + pow(ecal_y[icl] - ecal_y[jcl], 2));
            if (deltaR < 0.3) continue;                                // distance between clusters less than 0.3 m

            // timing
            if (clus_a_time[icl] < 80 || clus_a_time[icl] > 120) continue;  // cluster time in the range 80-120 ns
            if (clus_a_time[jcl] < 80 || clus_a_time[jcl] > 120) continue;  // cluster time in the range 80-120 ns
            // if (fabs(dt) > 40) continue; // ns time difference between clusters
            // if (fabs(dt) < fabs(best_dt) && fabs(dt) <= 50)
            double dt = clus_a_time[icl] - clus_a_time[jcl];

            // directions and mass
            TVector3 dir1 = (TVector3(ecal_x[icl], ecal_y[icl], z_calo) - vertex).Unit();
            TVector3 dir2 = (TVector3(ecal_x[jcl], ecal_y[jcl], z_calo) - vertex).Unit();

            TLorentzVector ph1(dir1.X()*E1, dir1.Y()*E1, dir1.Z()*E1, E1);
            TLorentzVector ph2(dir2.X()*E2, dir2.Y()*E2, dir2.Z()*E2, E2);

            double opening_angle = dir1.Angle(dir2)* (180.0 / TMath::Pi()); 
            double M  = (ph1 + ph2).M();
            double theta_exp = TMath::ACos(1-(pi0_mass_pdg*pi0_mass_pdg) / (2*sqrt(E1*E2)));

            
            //Tweak later from data
            double sigmaM = 0.04;   // GeV mass resolution 
            double sigmaT = 2.0;    // ns timing resolution 

            double sigmaTheta = 1;     // degrees if wTheta>0


            double wM         = 1.0;      // weight for mass term
            double wT         = 1.0;      // weight for timing term
            double wTheta     = 0.0;      // weight for theta-consistency
            

            // chi2 terms
            double best_chi2 = 1e10;
            double chi2_m  = pow(((M - pi0_mass_pdg)/sigmaM), 2);
            double chi2_t  = (dt/sigmaT) * (dt/sigmaT);
            double chi2_th   = pow(((opening_angle - theta_exp)/sigmaTheta), 2);
            double chi2 = wM*chi2_m + wT*chi2_t + wTheta*chi2_th;

            if (chi2 < best_chi2){
                best_chi2 = chi2;
                best_icl  = icl;
                best_jcl  = jcl;
                best_dt   = dt;
                best_mass = M;
            }
        }
    }
}

void My_ecal_pi0calib(int run_start, int run_end) {
    //Moved to global 
    // const Double_t z_calo = 6; // position of calorimeter from the target in m
    // const Double_t z_target = 0.09;    // position of target
    // const Double_t z_origin = 0.0;
    // const Double_t vertex_z = z_target - z_origin;  // position of vertex, where the pi0 is created, right in the middle of the target

    const int nbclusmax=100;    // maximum number of clusters
    const int sizemax=1000;     // maximum size of the cluster (how many blocks in each cluster)

    //Moved to global******* 
    // const double pi0_mass_pdg = 0.1349766;  // PDG pi0 mass in GeV
    TChain *ch = new TChain("T");

    ifstream infile("runlist.txt");
    string run_num;
    set<int> UncalibratedCh;
    while (getline(infile, run_num)) {
      if(stof(run_num)>=run_start && stof(run_num)<=run_end){
	ch->Add(Form("/volatile/halla/sbs/jonesdc/output/reduced/rootfiles/gep5_fullreplay_nogems_%s_*.root", run_num.c_str()));
      }
    }

 
    // Initialize the pointers before setting the branch address
    // Define the  energy, x, y of each cluster
    Double_t ecal_e[nbclusmax]; 
    Double_t ecal_x[nbclusmax];
    Double_t ecal_y[nbclusmax];

    ch->SetBranchStatus("*", 0); // disable all Branches
    ch->SetBranchAddress("earm.ecal.clus.e", &ecal_e);
    ch->AddBranchToCache("earm.ecal.clus.e", kTRUE);
    ch->SetBranchAddress("earm.ecal.clus.x", &ecal_x);
    ch->AddBranchToCache("earm.ecal.clus.x", kTRUE);
    ch->SetBranchAddress("earm.ecal.clus.y", &ecal_y);
    ch->AddBranchToCache("earm.ecal.clus.y", kTRUE);

    // Define the number of clusters
    Double_t nclus = 0;
    ch->SetBranchAddress("earm.ecal.nclus", &nclus);
    ch->AddBranchToCache("earm.ecal.nclus", kTRUE);
    // Define the time of each cluster
    Double_t clus_a_time[1000];
    ch->SetBranchAddress("earm.ecal.clus.atimeblk", &clus_a_time);
    ch->AddBranchToCache("earm.ecal.clus.atimeblk", kTRUE);

    // Define the number of blocks in each cluster and the id of each block per cluster (max 100 blocks per cluster)
    Double_t clus_nblk[nbclusmax];       // number of blocks in the cluster
    Double_t clus_id[nbclusmax];         // block id
    Double_t clus_row[nbclusmax];        // block row in the calorimeter
    Double_t clus_col[nbclusmax];        // block column in the calorimeter
    //Double_t clus_eblk[nbclusmax];       // block energy with highest energy in the cluster
    //Double_t clus_idblk[nbclusmax];      // block id with highest energy in the cluster
    Double_t goodblock_e[1000];
    Double_t goodblock_id[10000];
    Double_t goodblock_col[10000];
    Double_t goodblock_row[10000];
    Double_t goodblock_cid[10000];
    //Double_t energy_blk[1000];

    ch->SetBranchAddress("earm.ecal.clus.nblk", &clus_nblk);
    ch->AddBranchToCache("earm.ecal.clus.nblk", kTRUE);
    ch->SetBranchAddress("earm.ecal.clus.id", &clus_id);
    ch->AddBranchToCache("earm.ecal.clus.id", kTRUE);
    ch->SetBranchAddress("earm.ecal.clus.row", &clus_row);
    ch->AddBranchToCache("earm.ecal.clus.row", kTRUE);
    ch->SetBranchAddress("earm.ecal.clus.col", &clus_col);
    ch->AddBranchToCache("earm.ecal.clus.col", kTRUE);
    //ch->SetBranchAddress("earm.ecal.clus.eblk", &clus_eblk);
    //ch->SetBranchAddress("earm.ecal.idblk", &clus_idblk);
    ch->SetBranchAddress("earm.ecal.goodblock.e", &goodblock_e);
    ch->AddBranchToCache("earm.ecal.goodblock.e", kTRUE);
    ch->SetBranchAddress("earm.ecal.goodblock.id", &goodblock_id);
    ch->AddBranchToCache("earm.ecal.goodblock.id", kTRUE);
    ch->SetBranchAddress("earm.ecal.goodblock.col", &goodblock_col);
    ch->AddBranchToCache("earm.ecal.goodblock.col", kTRUE);
    ch->SetBranchAddress("earm.ecal.goodblock.row", &goodblock_row);
    ch->AddBranchToCache("earm.ecal.goodblock.row", kTRUE);
    ch->SetBranchAddress("earm.ecal.goodblock.cid", &goodblock_cid);
    ch->AddBranchToCache("earm.ecal.goodblock.cid", kTRUE);


    Double_t ngoodADChits = 0;
    ch->SetBranchAddress("earm.ecal.ngoodADChits", &ngoodADChits);
    ch->AddBranchToCache("earm.ecal.ngoodADChits",kTRUE);
    TH1F *h_pi0_mass = new TH1F("h_pi0_mass", "Pi0 Invariant Mass;M_{#pi^{0}} [GeV];Events", 80, 0, 0.6);
    TH1F *h_pi0_mass_corr = new TH1F("h_pi0_mass_reco", "Reconstructed #pi^{0} Invariant Mass;M_{#pi^{0}} [GeV];Events", 80, 0, 0.6);

    //Pi0 energy spectra
    TH1F *h_Epi0      = new TH1F("h_Epi0",
    "Pi^{0} Energy;E_{#pi^{0}} [GeV];Events", 140, 0, 14);
    TH1F *h_Epi0_corr = new TH1F("h_Epi0_corr",
    "Reconstructed Pi^{0} Energy;E_{#pi^{0}} [GeV];Events", 140, 0, 14);


    Long64_t nEvents = ch->GetEntries();
    cout << "Number of events: " << nEvents << endl;


    // --- Discover geometry from data: collect all unique block IDs ---
    int minRow = INT_MAX, maxRow = INT_MIN;
    int minCol = INT_MAX, maxCol = INT_MIN;
    std::set<int> blockIDs;
    std::map<int, int> blockID_to_row, blockID_to_col;

    for (Long64_t i = 0; i < nEvents/100; ++i) {
        ch->GetEntry(i);
        if (nclus < 1) continue;
        for (unsigned int b = 0; b < (unsigned int)ngoodADChits; ++b) {
            int bid = static_cast<int>(goodblock_id[b]);
            int row = static_cast<int>(goodblock_row[b]);
            int col = static_cast<int>(goodblock_col[b]);
            if (bid >= 0) {
                blockIDs.insert(bid);
                blockID_to_row[bid] = row;
                blockID_to_col[bid] = col;
                if (row < minRow) minRow = row;
                if (row > maxRow) maxRow = row;
                if (col < minCol) minCol = col;
                if (col > maxCol) maxCol = col;
            }
        }
    }
    int nlin = maxRow - minRow + 1;
    int ncol = maxCol - minCol + 1;

    /*for (Long64_t i = 0; i < nEvents/100; ++i) {
        ch->GetEntry(i);
        if (nclus < 1) continue;
        // All good blocks (only active blocks)
        for (unsigned int b = 0; b < (unsigned int)ngoodADChits; ++b) {
            int bid = static_cast<int>(goodblock_id[b]);
            if (bid >= 0) blockIDs.insert(bid);
        }
    }*/

    // Build mapping: blockID <-> matrix index
    std::map<int, int> blockID_to_idx;
    std::vector<int> idx_to_blockID;
    int idx = 0;
    for (int bid : blockIDs) {
        blockID_to_idx[bid] = idx++;
        idx_to_blockID.push_back(bid);
    }
    int nblocks = idx_to_blockID.size();

    cout << "Number of unique blocks: " << nblocks << endl;

    // Matrix used for calibration
    TMatrixD A(nblocks, nblocks);
    TVectorD B(nblocks);
    A.Zero();
    B.Zero();

    TH1D *hBestDt = new TH1D("hBestDt", "Time Difference between Photon Clusters",500,-50,50);
    TH1D *hBestDE = new TH1D("hBestDE", "Energy Difference between Photon Clusters",200,-3,3);
    TH2D *hEvsE = new TH2D("hEvsE", "Photon 2 Energy_{Cluster} vs Photon 1 E_{Cluster}",200,0,7,200,0,7);
    TH1D *hBestDx = new TH1D("hBestDx", "Distance between Photon Clusters",50,0,5);
    TH1D *hAngle = new TH1D("hAngle", "Angle between Photon Clusters",100,0,30);
    TH1D *hAngleCut = new TH1D("hAngleCut", "Angle between Photon Clusters After Cut",200,-30,30);
    hAngleCut->SetLineColor(kRed);
  
    // Block occupancy for later use
    std::vector<int> occupancy(nblocks, 0);

    // Define the counters for the cuts
    Double_t pass_clus=0, pass_deltar=0, pass_time=0, pass_mass=0;
    // initialize the random number generator
    gRandom->SetSeed(0);
    //Double_t pi0_mass_pdg_random = pi0_mass_pdg;    
    std::vector<std::vector<std::vector<int>>> clusterBlocks(nEvents);

    for (Long64_t i = 0; i < nEvents; i++) {
        ch->GetEntry(i);
        if (i % 10000 == 0) {
            cout << "Processing event (correction) " << i << " / " << nEvents << "\r";
            cout.flush();
        }
        // Check that there are at least two clusters
        if (nclus < 2) continue;
        pass_clus++;

        // Find the two clusters with the best time difference and apply cuts
        // double best_dt = 1e9;
        // int best_icl = -1, best_jcl = -1;
        // best_pair_cuts(ecal_e, clus_nblk, ecal_x, ecal_y, clus_a_time, nclus, best_icl, best_jcl, best_dt);

        int best_icl=-1, best_jcl=-1;
        double best_dt=0.0, best_mass=0.0;
        // pick one pair per event by chi2
        best_pair_chi2(ecal_e, clus_nblk, ecal_x, ecal_y, clus_a_time, nclus, best_icl, best_jcl, best_dt, best_mass);

        pass_time++;

	if(best_icl >=0 && best_jcl >=0){
	  hBestDE->Fill(ecal_e[best_icl]-ecal_e[best_jcl]);
	  hBestDx->Fill(sqrt(pow(ecal_x[best_icl]-ecal_x[best_jcl],2)+pow(ecal_y[best_icl]-ecal_y[best_jcl],2)));
	}
	hBestDt->Fill(best_dt);
	if(fabs(best_dt) > 4) continue;
        // Find indices of the two clusters with the highest energies
        int imax1 = -1, imax2 = -1; // Index of highest and second-highest energy cluster
        double e1 = -1, e2 = -1;  // Energy of highest and second-highest energy cluster
        find_highest_energy(ecal_e, nclus, imax1, imax2, e1, e2);

        if (best_icl >= 0 && best_jcl >= 0) {
            if (!((best_icl == imax1 && best_jcl == imax2) || (best_icl == imax2 && best_jcl == imax1))) {
                continue; // skip if not the two highest-energy clusters
            }
            TVector3 pos1(ecal_x[best_icl], ecal_y[best_icl], z_calo);  // in m
            TVector3 pos2(ecal_x[best_jcl], ecal_y[best_jcl], z_calo);  // in m

            TVector3 vertex(0, 0, z_target);
            TVector3 dir1 = (pos1 - vertex).Unit();
            TVector3 dir2 = (pos2 - vertex).Unit();

            TLorentzVector ph1(dir1.X() * ecal_e[best_icl], dir1.Y() * ecal_e[best_icl], dir1.Z() * ecal_e[best_icl], ecal_e[best_icl]);
            TLorentzVector ph2(dir2.X() * ecal_e[best_jcl], dir2.Y() * ecal_e[best_jcl], dir2.Z() * ecal_e[best_jcl], ecal_e[best_jcl]);


            Double_t pi0_mass = (ph1 + ph2).M();

            Double_t opening_angle = dir1.Angle(dir2) * (180.0 / TMath::Pi());
	    hAngle->Fill(opening_angle);
            // if (opening_angle < 3.5 || opening_angle > 15) continue;  // The lower cut (e.g., < 6°) removes nearly collinear photon pairs → likely merged.
                                                                    // The upper cut (e.g., > 80°) removes highly unphysical, possibly misreconstructed pairs.
	    hAngleCut->Fill(opening_angle);
	    
            if (pi0_mass < 0.02 || pi0_mass > 0.6) continue;
            pass_mass++;
            h_pi0_mass->Fill(pi0_mass);

            double Epi0 = (ph1 + ph2).E(); // == ecal_e[best_icl] + ecal_e[best_jcl]
            h_Epi0->Fill(Epi0);


            // compute expected total π0 energy as using the scale factor
            //double pi0_mass_smeared = gRandom->Gaus(pi0_mass_pdg, 0.005); // 5 MeV width max
            double expected_E = (ecal_e[best_icl] + ecal_e[best_jcl]) * (pi0_mass_pdg / pi0_mass);
            //double expected_E = pi0_mass_pdg / pi0_mass;
            hEvsE->Fill(ecal_e[best_icl],ecal_e[best_jcl]);
            // **Loop over all blocks in cluster **:
            // zero the per-block energy accumulator
            static std::vector<double> energy(nblocks, 0.0);
            energy.assign(nblocks, 0.0);
    
            // Loop over all clusters and blocks per event
            for (int cl = 0; cl < nclus; ++cl) {
	      if(!(cl == best_icl || cl == best_jcl))continue;
                // Assign energy to each good block
                for (unsigned int b = 0; b < (unsigned int)ngoodADChits; ++b) {
                    int rawID = (int)goodblock_id[b];
                    int cid = (int)goodblock_cid[b];
                    //int iblock = num_btom[rawID];
                    int iblock = blockID_to_idx.count(rawID) ? blockID_to_idx[rawID] : -1;
                    if (iblock < 0 || iblock >= nblocks) continue;
                    if (cid != cl) continue; // Only assign energy to blocks in this cluster
                    energy[iblock] += goodblock_e[b];
                    if (energy[iblock] > 0.0) occupancy[iblock]++;
                }
            }
            // Now fill A and B:
            for (int j = 0; j < nblocks; ++j) {
                for (int k = 0; k < nblocks; ++k) {
                    A(j,k) += energy[j] * energy[k]/(expected_E*expected_E);
                    // A(j,k) += energy[j] * energy[k]/expected_E;
                }
                B(j) += energy[j]/expected_E;
                // B(j) += energy[j];
            } 
        }
    }
    cout<<"Events passing number of clusters cut: "<<pass_clus<<endl;
    //cout<<"Events passing deltaR cut: "<<pass_deltar<<endl;
    cout<<"Events passing time cut: "<<pass_time<<endl;
    cout<<"Events passing all cuts: "<<pass_mass<<endl;
  
    // ===================== OCCUPANCY PRUNING & REGULARIZATION =====================
    const int N_min = 30;
    for (int i = 0; i < nblocks; ++i) {
        if (occupancy[i] < N_min) {
            // zero out row and column i
            for (int j = 0; j < nblocks; ++j) {
                A(i,j) = A(j,i) = 0.0;
            }
            A(i,i) = 1.0;
            B(i) = 1.0;
	    UncalibratedCh.insert(i);
        }
    }
    // add small diagonal term λ·A_ii for singularity
    const double lambda = 0;//1e-4;
    // const double lambda = 1e-4;
    for (int i = 0; i < nblocks; ++i) {
        A(i,i) += lambda * A(i,i);
    }
        
    // -------------------- Invert the Matrix --------------------
    //A.Print();
    std::cout
        <<"[DEBUG] A(0,0)="<<A(0,0)
        <<", B(0,0)="<<B(0)
        <<"\n";
        for (int i = 0; i < std::min(5,nblocks); ++i) {
        std::cout<<"  B(0,"<<i<<")="<<B(i)<<"\n";
        }

    //Double_t coeff[nblocks];
    std::vector<double> coeff(nblocks, 1.0);

    // SVD decomposition of A and solve C = A⁺·B (see https://root.cern.ch/doc/master/classTDecompSVD.html)
    TDecompSVD svd(A);
    if (!svd.Decompose()) {
      cerr<<"ERROR: SVD decomposition failed\n";
      return;
    }
    // Prepare solution vector
    // Copy B into C, then solve C = A⁺·B
    cout << "\n--- Solving the matrix system ---" << endl;
    TVectorD C = B;           // initialize C with RHS
    Bool_t ok = svd.Solve(C); // Solve in place: C ← pseudo-inverse(A)·C
    if (!ok) {
      cerr<<"ERROR: SVD Solve failed\n";
      return;
    }

    for (int i = 0; i < min(nblocks,5); ++i) {
      cout << "C["<<i<<"] = " << C(i) << "\n";
    }
        // the coefficients of the linear system 
        // C = A^-1 * B where A^-1 is the inverse of A, B is the vector of unknowns, and C is the vector of calibration coefficients. originally A * C = B
        //TMatrixD C = B * A;  // 1×nblocksm * nblocksm×nblocksm

    // now overwrite just the “good” ones using occupancy
    for (int compID = 0; compID < nblocks; ++compID) {
        //int rawID = num_mtob[compID];     // reverse map
        //int blockID = idx_to_blockID[compID];
        //double c = C(0,compID);
        double c = C(compID);

        if (occupancy[compID] < N_min || c < 1e-2 || c > 10) {
            coeff[compID] = 1.0;
        } else {
            coeff[compID] = c;
        //coeff[rawID] = C(0, compID);      // use same compID as above
        }
    }
    
   
    // Now, with the coefficients in hand, I can use them to calibrate the Ecal energy in the code.
    // Use the coefficients to calibrate the Ecal energy and plot the pi0 invariant mass
    // with and without the correction

    // counters for passing events
    Double_t pass_clus_corr=0, pass_deltar_corr=0, pass_time_corr=0, pass_mass_corr=0;
    double test_mass_sum = 0;
    int test_mass_count = 0;
 
    for (Long64_t i = 0; i < nEvents; i++) {
        ch->GetEntry(i);
        if (i % 10000 == 0) {
            cout << "Processing event (correction) " << i << " / " << nEvents << "\r";
            cout.flush();
        }
        if (nclus < 2) continue;
        //pass_clus_corr++;
         // Find the two clusters with the best time difference and apply cuts
        double best_dt = 1e9;
        int best_icl = -1, best_jcl = -1;
        best_pair_cuts(ecal_e, clus_nblk, ecal_x, ecal_y, clus_a_time, nclus, best_icl, best_jcl, best_dt);
        //pass_time_corr++;
        // Find indices of the two clusters with the highest energies
        int imax1 = -1, imax2 = -1; // Index of highest and second-highest energy cluster
        double e1 = -1, e2 = -1;  // Energy of highest and second-highest energy cluster
        find_highest_energy(ecal_e, nclus, imax1, imax2, e1, e2);
 
        // **After** correction: rebuild each photon’s energy by summing
        //     per‑block energies * coefficients
        if (best_icl >= 0 && best_jcl >= 0) {
            if (!((best_icl == imax1 && best_jcl == imax2) || (best_icl == imax2 && best_jcl == imax1))) {
                continue; // skip if not the two highest-energy clusters
            }
            double e_corr[2] = {0., 0.};

            int idx[2] = {best_icl, best_jcl};
            for (int ic = 0; ic < 2; ++ic) {
                int cl = idx[ic];

                // Assign energy to each good block
                for (unsigned int b = 0; b < (unsigned int)ngoodADChits; ++b) {
                    int rawID = (int)goodblock_id[b];
                    int cid = (int)goodblock_cid[b];
                    // Use blockID_to_idx to get the matrix index
                    auto it = blockID_to_idx.find(rawID);
                    if (it == blockID_to_idx.end()) continue; // skip if not a known block
                    if (cid != cl) continue; // Only assign energy to blocks in this cluster
                    e_corr[ic] += coeff[it->second] * goodblock_e[b];
                }
            }
            TVector3 pos1(ecal_x[best_icl], ecal_y[best_icl], z_calo);  // in m
            TVector3 pos2(ecal_x[best_jcl], ecal_y[best_jcl], z_calo);  // in m

            TVector3 vertex(0, 0, z_target);
            TVector3 dir1 = (pos1 - vertex).Unit();
            TVector3 dir2 = (pos2 - vertex).Unit();

            // rebuild corrected TLorentzVectors
            TLorentzVector ph1_corr(dir1.X() * e_corr[0], dir1.Y() * e_corr[0], dir1.Z() * e_corr[0], e_corr[0]);
            TLorentzVector ph2_corr(dir2.X() * e_corr[1], dir2.Y() * e_corr[1], dir2.Z() * e_corr[1], e_corr[1]);
            
            double pi0_mass_corr = (ph1_corr + ph2_corr).M();

            Double_t opening_angle = dir1.Angle(dir2) * (180.0 / TMath::Pi());
            // if (opening_angle < 3.5 || opening_angle > 8) continue;  // The lower cut (e.g., < 6°) removes nearly collinear photon pairs → likely merged.
            // if (opening_angle < 3.5 || opening_angle > 15) continue;  // The lower cut (e.g., < 6°) removes nearly collinear photon pairs → likely merged.
            // The upper cut (e.g., > 80°) removes highly unphysical, possibly misreconstructed pairs.

            if (pi0_mass_corr <= 0.02 || pi0_mass_corr >= 0.4) continue;    // Making a cut on the pi0 mass, e.g., between 0.06 and 0.6 GeV

            // fill histograms
            test_mass_sum += pi0_mass_corr;
            test_mass_count++;
            //pass_mass_corr++;
            //h_pi0_mass_corr->Fill(pi0_mass_corr);
        }
    }
    cout << "\n -----  Scaling coefficients to match PDG π⁰ mass.. -----" << endl;
    if (test_mass_sum > 0) {
        double mean_mass = test_mass_sum / test_mass_count;
        double scale_factor = pi0_mass_pdg / mean_mass;

        // Apply the scaling factor to the coefficients
        for (int i = 0; i < nblocks; ++i) {
            coeff[i] *= scale_factor;
        }
    }
    //h_pi0_mass_corr->Reset();

     // -------------------- Saving the coefficients in a text file --------------------
    cout << "\n--- Creating the coefficient file ---" << endl;
    cout<<"Writing coefficients to file..."<<endl;

    ofstream coeff_file("calib_coeff/ecal_block_calibration_factors_noiter.txt");
    coeff_file << "# BlockID\tCalibrationCoeff\n";

    for (int i = 0; i < nblocks; i++) {
        int blockID = idx_to_blockID[i];
        coeff_file << blockID << "\t" << coeff[i] << endl;
        cout << "BlockID: " << i << "\tCalibrationCoeff: " << coeff[i] << endl;
    }

    coeff_file.close();


    for (Long64_t i = 0; i < nEvents; i++) {
        ch->GetEntry(i);
        cout << "Processing event " << i << "\r";
        if (nclus < 2) continue;
        pass_clus_corr++;
         // Find the two clusters with the best time difference and apply cuts
        double best_dt = 1e9;
        int best_icl = -1, best_jcl = -1;
        best_pair_cuts(ecal_e, clus_nblk, ecal_x, ecal_y, clus_a_time, nclus, best_icl, best_jcl, best_dt);
        pass_time_corr++;
        // Find indices of the two clusters with the highest energies
        int imax1 = -1, imax2 = -1; // Index of highest and second-highest energy cluster
        double e1 = -1, e2 = -1;  // Energy of highest and second-highest energy cluster
        find_highest_energy(ecal_e, nclus, imax1, imax2, e1, e2);
 
        // **After** correction: rebuild each photon’s energy by summing
        //     per‑block energies * coefficients
        if (best_icl >= 0 && best_jcl >= 0) {
            if (!((best_icl == imax1 && best_jcl == imax2) || (best_icl == imax2 && best_jcl == imax1))) {
                continue; // skip if not the two highest-energy clusters
            }
            double e_corr[2] = {0., 0.};

            int idx[2] = {best_icl, best_jcl};
            for (int ic = 0; ic < 2; ++ic) {
                int cl = idx[ic];

                // Assign energy to each good block
                for (unsigned int b = 0; b < (unsigned int)ngoodADChits; ++b) {
                    int rawID = (int)goodblock_id[b];
                    int cid = (int)goodblock_cid[b];
                    // Use blockID_to_idx to get the matrix index
                    auto it = blockID_to_idx.find(rawID);
                    if (it == blockID_to_idx.end()) continue; // skip if not a known block
                    if (cid != cl) continue; // Only assign energy to blocks in this cluster
                    e_corr[ic] += coeff[it->second] * goodblock_e[b];
                }
            }
            TVector3 pos1(ecal_x[best_icl], ecal_y[best_icl], z_calo);  // in m
            TVector3 pos2(ecal_x[best_jcl], ecal_y[best_jcl], z_calo);  // in m

            TVector3 vertex(0, 0, z_target);
            TVector3 dir1 = (pos1 - vertex).Unit();
            TVector3 dir2 = (pos2 - vertex).Unit();

            // rebuild corrected TLorentzVectors
            TLorentzVector ph1_corr(dir1.X() * e_corr[0], dir1.Y() * e_corr[0], dir1.Z() * e_corr[0], e_corr[0]);
            TLorentzVector ph2_corr(dir2.X() * e_corr[1], dir2.Y() * e_corr[1], dir2.Z() * e_corr[1], e_corr[1]);
            
            double pi0_mass_corr = (ph1_corr + ph2_corr).M();

            Double_t opening_angle = dir1.Angle(dir2) * (180.0 / TMath::Pi());
            if (opening_angle < 3.5 || opening_angle > 8) continue;  // The lower cut (e.g., < 6°) removes nearly collinear photon pairs → likely merged.
                                                                    // The upper cut (e.g., > 80°) removes highly unphysical, possibly misreconstructed pairs.

            if (pi0_mass_corr <= 0.02 || pi0_mass_corr >= 0.4) continue;    // Making a cut on the pi0 mass, e.g., between 0.06 and 0.6 GeV
            h_pi0_mass_corr->Fill(pi0_mass_corr);
            
            //fill pi0 E after
            double Epi0_corr = (ph1_corr + ph2_corr).E(); // == e_corr[0] + e_corr[1]
            h_Epi0_corr->Fill(Epi0_corr);
        }
    }

    // -------------------- Plotting --------------------
    TCanvas *c = new TCanvas("c","Reconstructed #pi^{0} Invariant Mass Before and After Calibration", 800,600);
    h_pi0_mass->SetStats(0);
    h_pi0_mass_corr->SetStats(0);
    h_pi0_mass->GetXaxis()->SetTitle("M_{#pi^{0}} [GeV]");
    h_pi0_mass->GetYaxis()->SetTitle("Events");
    h_pi0_mass->SetLineWidth(2);
    h_pi0_mass_corr->SetLineWidth(2);
    h_pi0_mass->SetLineColor(kBlue);
    h_pi0_mass_corr->SetLineColor(kRed);
    h_pi0_mass_corr->Draw();
    h_pi0_mass->Draw("SAME");
    auto leg = new TLegend(0.6,0.7,0.9,0.9);
    leg->AddEntry(h_pi0_mass,      "Before calib", "l");
    leg->AddEntry(h_pi0_mass_corr, "After calib",  "l");
    leg->Draw();
    c->SaveAs(Form("plots/ecal_pi0_mass_calib_noiter%i_%i.png",run_start,run_end));


    // Create a 2D histogram to visualize calibration coefficients in the detector geometry
    TH2F *h_coeff_map = new TH2F("h_coeff_map", "ECAL Block Calibration Coefficients;Column;Row",
        ncol, minCol, maxCol + 1l,
        nlin, minRow, maxRow + 1);
    TGraph *gr = new TGraph();
    TH1D *hCoeff = new TH1D("hCoeff","Calibration Coefficient Distribution",100,0,10);
    // Fill the 2D histogram using the coeff[] array 
    for (int i = 0; i < nblocks; ++i) {
        int blockID = idx_to_blockID[i];
        int row = blockID_to_row[blockID];
        int col = blockID_to_col[blockID];
        double val = coeff[i];
        h_coeff_map->SetBinContent(col - minCol + 1, row - minRow + 1, val);  // ROOT bins start from 1
	gr->SetPoint(i,blockID,val);
	hCoeff->Fill(val);
    }
    cout << "\n -----  Summary of the cuts -----" << endl;
    cout<<"Events passing number of clusters cut after calib: "<<pass_clus_corr<<endl;
    //cout<<"Events passing deltaR cut after calib: "<<pass_deltar_corr<<endl;
    cout<<"Events passing time cut after calib: "<<pass_time_corr<<endl;
    cout<<"Events passing all cuts after calib: "<<pass_mass_corr<<endl;

    cout<<"-------------------- Difference between before and after calib -------------------"<<endl;
    cout << "Before calib: mean = " << h_pi0_mass->GetMean() 
        << ", sigma = " << h_pi0_mass->GetRMS() << endl;
    cout << "After calib:  mean = " << h_pi0_mass_corr->GetMean() 
        << ", sigma = " << h_pi0_mass_corr->GetRMS() << endl;

    // Draw the map
    TCanvas *c_map = new TCanvas("c_map", "ECAL Coefficients Heatmap", 1600, 900);
    c_map->Divide(2,2);
    c_map->cd(1);
    h_coeff_map->SetStats(0);
    h_coeff_map->Draw("COLZ");
    c_map->cd(2);
    gr->SetTitle("Calibration Coefficients vs. Block ID");
    gr->SetMarkerStyle(6);
    gr->Draw("ap");
    gr->GetXaxis()->SetTitle("Block ID");
    gr->GetYaxis()->SetTitle("Coefficient");
    gPad->Update();
    c_map->cd(3);
    hCoeff->Draw();
    c_map->cd(4);
    hEvsE->Draw("colz");
    c_map->SaveAs(Form("plots/ecal_coefficients_heatmap_noiter%i_%i.png",run_start,run_end));
    TCanvas *cDt = new TCanvas("cDt","Photon Differences",0,0,1600,900);
    cDt->Divide(2,2);
    cDt->cd(1);
    hBestDt->Draw();
    cDt->cd(2);
    hBestDx->Draw();
    cDt->cd(3);
    hAngle->Draw();
    hAngleCut->Draw("sames");
    cout<<UncalibratedCh.size()<<" uncalibrated channels:"<<endl;
    cDt->cd(4);
    hBestDE->Draw();
    cDt->SaveAs(Form("plots/randomPlots%i_%i.png",run_start,run_end));


    TCanvas *cE = new TCanvas("cE","Pi0 Energy Before/After Calibration",800,600);
    h_Epi0->SetStats(0);        h_Epi0_corr->SetStats(0);
    h_Epi0->SetLineWidth(2);    h_Epi0_corr->SetLineWidth(2);
    h_Epi0->SetLineColor(kBlue);
    h_Epi0_corr->SetLineColor(kRed);

    h_Epi0_corr->Draw();        // draw after-calib first
    h_Epi0->Draw("SAME");

    auto legE = new TLegend(0.6,0.7,0.9,0.9);
    legE->AddEntry(h_Epi0,      "Before calib", "l");
    legE->AddEntry(h_Epi0_corr, "After calib",  "l");
    legE->Draw();

    cE->SaveAs(Form("plots/ecal_pi0_energy_%i_%i.png",run_start,run_end));

    for(int val : UncalibratedCh)cout<<val<<", ";
    cout<<endl;
}
