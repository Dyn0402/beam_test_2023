#include <iostream>
#include <iomanip>
#include <fstream>
#include <bitset>
#include <vector>
#include <map>
#include <array>

#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TTreeReader.h"

#include "definition_banco.h"

// takes the path to the noise root file and returns a vector with the pixId of all noisy pixels
std::vector<int> noiseVect(std::string fname){

	int threshold = 0;
	
	TFile* fin = TFile::Open(fname.c_str(),"read");

	TTreeReader reader("pixTree",fin);
	TTreeReaderValue<UInt_t> trgNum(reader, "fData.trgNum");
	TTreeReaderValue<UInt_t> row(reader, "fData.row");
	TTreeReaderValue<UInt_t> chipCol(reader, "fData.col");
	TTreeReaderValue<UInt_t> ladderId(reader, "fData.deviceType");
	TTreeReaderValue<UInt_t> chipId(reader, "fData.chipId");

	// pixId maps where 
	std::vector<int> noise;
	std::map<int, int> count;

	while(reader.Next()){

		int col = *chipCol + (*chipId - 4)*1024;
		int pixId = *row*1e4 + col;
		
		if(count.find(pixId) == count.end()){
			count[pixId] = 1;
		}else{
			count[pixId]++;
		}
	}

	for(std::map<int,int>::iterator iter = count.begin(); iter != count.end(); ++iter){
		if(iter->second > threshold){
			noise.push_back(iter->first);
			std::cout<<"noise pixId: "<<iter->first<<std::endl;
		}
	}

	return noise;
}

void buildCluster(int pixId, std::vector<hitBanco> *hits, std::vector<hitBanco> *hitsInCl, int clId){
	//loop of rows then cols
	for(int i=pixId-1e4; i<pixId+2e4; i+=1e4){
		for(int j=i-1; j<i+2; j++){
			auto it = std::find_if( hits->begin(), hits->end(), [j](const hitBanco& hit) {return hit.pixId == j;} );

			if(it!=hits->end() and it->clusterId<0){
				it->clusterId = clId;
				hitsInCl->push_back(*it);
				buildCluster(it->pixId, hits, hitsInCl, clId);
			}
		}
	}
}

void findClusters(std::vector<hitBanco> *hits, std::vector<clusterBanco> *cls){
	
	std::vector<hitBanco> *hitsInCl = new std::vector<hitBanco>();
	clusterBanco acl;
	int clId = -1;
	for(auto h = hits->begin(); h != hits->end(); ++h){
		if(h->clusterId > -1) continue;
		hitsInCl->clear();
		clId++;
		h->clusterId = clId;
		hitsInCl->push_back(*h);
		buildCluster(h->pixId, hits, hitsInCl, clId);
		// std::cout<<h.pixId<<" cluser size "<<hitsInCl->size()<<std::endl;

		float colCentroid = 0, rowCentroid = 0;
		
		for(auto hc : *hitsInCl){
			colCentroid += hc.col;
			rowCentroid += hc.row;
		}
		acl.colCentroid = colCentroid/hitsInCl->size();
		acl.rowCentroid = rowCentroid/hitsInCl->size();
		acl.size = hitsInCl->size();
		acl.ladderId = hitsInCl->at(0).ladderId;
		acl.id = clId;
		cls->push_back(acl);
	}
}


// takes the alpide root file and outputs a root file of hitBanco denoised where event number matches trigger number 
void decodeBanco(std::string fnameIn, std::string fnameOut, std::string fnameNoise){

	// noise vector
	std::vector<int> noise = noiseVect(fnameNoise);

	// input file
	TFile* fin = TFile::Open(fnameIn.c_str(),"read");
	TTree* pixTree = (TTree*)fin->Get("pixTree");
	pixTree->SetMakeClass(1);

	std::cout<<"Reading file: "<<fnameIn<<std::endl;
	pixTree->Print();

	TTreeReader reader("pixTree",fin);
	TTreeReaderValue<UInt_t> trgNum(reader, "fData.trgNum");
	TTreeReaderValue<UInt_t> row(reader, "fData.row");
	TTreeReaderValue<UInt_t> chipCol(reader, "fData.col");
	TTreeReaderValue<UInt_t> ladderId(reader, "fData.deviceType");
	TTreeReaderValue<UInt_t> chipId(reader, "fData.chipId");

	// output file
	TFile* fout = TFile::Open(fnameOut.c_str(),"recreate");
	TTree* outTree = new TTree("events","");

	int eventId = 0;
	std::vector<hitBanco> *hits = new std::vector<hitBanco>();
	std::vector<clusterBanco> *cls = new std::vector<clusterBanco>();
	
	outTree->Branch( "eventId", &eventId );
	outTree->Branch( "hits", &hits );
	outTree->Branch( "clusters", &cls );

	hitBanco ahit;
	int currentTrg = 0;
	int clId = 0;

	while(reader.Next()){
		
		int col = *chipCol + (*chipId - 4)*1024;
		int pixId = *row*1e4 + col;
		if(std::find(noise.begin(), noise.end(), pixId) != noise.end()) continue;
		// std::cout<<*row<<" "<<col<<std::	endl;
		ahit.row = *row;
		ahit.chipCol = *chipCol;
		ahit.col = col;
		ahit.pixId = pixId;
		ahit.ladderId = *ladderId;
		ahit.chipId = *chipId;
		ahit.clusterId = -1;
		
		if(currentTrg != *trgNum){
			eventId = currentTrg;
			if(eventId%10000==0) std::cout<<eventId<<std::endl;
			findClusters(hits, cls);
			outTree->Fill();
			hits->clear();
			cls->clear();
			currentTrg ++;
			// in case of empty events
			while(currentTrg != *trgNum){
				currentTrg++;
				eventId = currentTrg;
				outTree->Fill();
			}
		}
		hits->push_back(ahit);
	}

	outTree->Write();
	fout->Close();
}

int main(int argc, char const *argv[])
{
	std::string fnameIn = argv[1];
	std::string fnameOut = fnameIn;
	std::string noise = "multinoiseScan_230610_102306_NOBEAM-B0-ladder163.root";
	fnameOut.insert(fnameIn.find(".root"), "_decoTh0");

	decodeBanco(fnameIn, fnameOut, noise);

	TFile* fin = TFile::Open(fnameOut.c_str(),"read");
	TTreeReader reader("events",fin);
	TTreeReaderValue<std::vector<hitBanco>> hits(reader, "hits");
	TTreeReaderValue<std::vector<clusterBanco>> cls(reader, "clusters");

	TH2F *h2 = new TH2F("h2", "Map hits", 1024,4*1024,5*1024, 512,0,512);
	h2->SetXTitle("rows");
	h2->SetYTitle("columns");

	TH2F *h2c = new TH2F("h2c", "Map clusters", 1024,4*1024,5*1024, 512,0,512);
	h2c->SetXTitle("rows");
	h2c->SetYTitle("columns");

	TH1F *hclsize = new TH1F("hclsize", "cluster size", 16,-0.5,15.5);
	hclsize->SetXTitle("cluster size");

	TH1F *hnbcl = new TH1F("hnbcl", "Number of clusters per event", 6,-0.5,5.5);
	hnbcl->SetXTitle("number of clusters");

	while(reader.Next()){
		for(auto h : *hits){
			h2->Fill(h.col, h.row);
		}
		hnbcl->Fill(cls->size());
		for(auto cl : *cls){
			h2c->Fill(cl.colCentroid, cl.rowCentroid);
			hclsize->Fill(cl.size);
		}
	}

	TCanvas *c2 = new TCanvas("c2", "c2", 1600,1000);
	// h2->SetStats(111);
	// h2c->SetStats(111);
	c2->Divide(1,2);
	c2->cd(1);
	h2->Draw("colz");
	// gPad->SetLogz();
	c2->cd(2);
	h2c->Draw("colz");
	c2->Print("mapTestTh0.png", "png");

	TCanvas *c3 = new TCanvas("c3", "c3", 1600,1000);
	c3->Divide(2,1);
	c3->cd(1);
	hnbcl->Draw();
	c3->cd(2);
	hclsize->Draw();
	gPad->SetLogy();
	c3->Print("clusterTestTh0.png", "png");

	return 0;
}