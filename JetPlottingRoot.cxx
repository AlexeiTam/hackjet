#include <iostream>
#include <cmath>


void JetPlottingRoot(const char* filename = "myfile.root", const char* treename = "TreeName", const char* outimage = "myjetimage.png"){

//open file    
TFile* file = TFile::Open(filename, "READ");

std::cout << filename << std::endl;
if (!file || file->IsZombie()){

    std::cout << "no file opened. breaking out" << std::endl;
    return;
}

//obtain tree
TTree *tree = (TTree*) file->Get(treename);
if (!tree){
    std::cout << "tree not located. breaking out." << std::endl;
}

//read branches:

//variables
double E;
double pT;
double eta;
double phi;

TBranch *bE = 0;
TBranch *bpT = 0;
TBranch *beta = 0;
TBranch *bphi = 0;

//Assigning branches

tree->SetBranchAddress("E", &E, &bE);
tree->SetBranchAddress("pT", &pT, &bpT);
tree->SetBranchAddress("eta", &eta, &beta);
tree->SetBranchAddress("phi", &phi, &bphi);

    //place branches into vectors, so we can obtain ranges for plotting

    std::vector<double> veta;
    std::vector<double> vphi;

for(int i = 0; i < tree->GetEntries();i++){
    tree->LoadTree(i);

    beta->GetEntry(i);
    bphi->GetEntry(i);

    veta.emplace_back(E);
    vphi.emplace_back(phi);
}


//creating histogram

const int NBins = 200;

double eta_max = *max_element(veta.begin(), veta.end());
double eta_min = *min_element(veta.begin(), veta.end());

double phi_max = *max_element(vphi.begin(), vphi.end());
double phi_min = *min_element(vphi.begin(), vphi.end());


TH2D *hist = new TH2D("hist", "pT Distribution", NBins, eta_min, eta_max, NBins, phi_min, phi_max);

hist->GetXaxis()->SetTitle("#eta");
hist->GetYaxis()->SetTitle("#phi");
hist->GetZaxis()->SetTitle("pT");

//preventing histogram from deleting upon file closure
//hist->SetDirectory(nullptr);
    
for(int i = 0; i < tree->GetEntries();i++){
    tree->LoadTree(i);

    beta->GetEntry(i);
    bphi->GetEntry(i);
    bpT->GetEntry(i);

    hist->Fill(eta,phi,pT);

    
}



// creating the canvas and saving to image

TCanvas *c1 = new TCanvas("c1","",900,900);

gStyle->SetOptStat(0);
hist->Draw("LEGO0");

c1->Update();
c1->SaveAs(outimage);
//note to self: may want to save output files to a root file

//file->Close();
return;


}
