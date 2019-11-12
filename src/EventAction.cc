#include "EventAction.hh"
#include "RunAction.hh"
#include "Analysis.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "DetectorConstruction.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction(RunAction* runAction)
: G4UserEventAction(),
fRunAction(runAction)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event* myEvent)
{
	fDetectedPhotons = 0;
	fDepositedEnergy = 0;
	fProducedPhotons = 0;
	fTopPhoton =0;
	fBottomPhoton=0;
	fSidePhoton=0;
	fEscapedPhoton = 0;
	fAbsorbedPhoton=0;
	if (myEvent->GetEventID() % 1000 == 0)
		G4cout << "event no.: " << myEvent->GetEventID() << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void EventAction::AddWavelength(G4double wavelength){
	auto analysisManager = G4AnalysisManager::Instance();

	analysisManager->FillH1(6, wavelength);
}

void EventAction::EndOfEventAction(const G4Event* myEvent)
{

	auto analysisManager = G4AnalysisManager::Instance();

	 if (fDetectedPhotons > 5){
		 analysisManager->FillNtupleDColumn(0,fDetectedPhotons);
		 analysisManager->FillNtupleDColumn(1,fDepositedEnergy);
		// G4cout << fDepositedEnergy << G4endl;
	 }

	 analysisManager->AddNtupleRow(0);

//}

}
