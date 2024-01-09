//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id$
//
/// \file RunAction.cc
/// \brief Implementation of the RunAction class

#include "RunAction.hh"
#include "TString.h"

#include "G4AnalysisManager.hh"
#include "G4Run.hh"
#include "G4SystemOfUnits.hh"

namespace ED
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(EventAction* eventAction) : fEventAction(eventAction)
{

  G4String label = "eic";
  G4String title = "Hits eic";
  G4String fileName = "beamTestGeant4.root";

  // Create analysis manager
  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetNtupleMerging(true);
  analysisManager->SetDefaultFileType("root");
  analysisManager->SetFileName(fileName);

  // analysisManager->CreateNtuple(label, title);
  // analysisManager->CreateNtupleIColumn("Layer");   // column id = 0
  // analysisManager->CreateNtupleDColumn("Xpos");    // column id = 1
  // analysisManager->CreateNtupleDColumn("Ypos");    // column id = 2
  // analysisManager->CreateNtupleDColumn("Zpos");    // column id = 3
  // analysisManager->FinishNtuple();

  if(fEventAction){
    analysisManager->CreateNtuple(label, title);
    analysisManager->CreateNtupleDColumn("x0");
    analysisManager->CreateNtupleDColumn("y0");
    analysisManager->CreateNtupleDColumn("mx");
    analysisManager->CreateNtupleDColumn("my");
    analysisManager->CreateNtupleDColumn("chi2x");
    analysisManager->CreateNtupleDColumn("chi2y");
    for(int i=0; i<5; i++){
      analysisManager->CreateNtupleDColumn(Form("MMpos%d", i), fEventAction->getMMpos(i));
    }
    analysisManager->FinishNtuple();
    // analysisManager->SetNtupleFileName(0, fileName);
    analysisManager->CreateH1("testGauss","smearing", 100, 0., 1.);
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* /*run*/)
{
  // Get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->Reset();
  analysisManager->OpenFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* /*run*/)
{
  // save histograms
  //
  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->Write();
  analysisManager->CloseFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}