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



file->ls();
file->Close();
return;


}
