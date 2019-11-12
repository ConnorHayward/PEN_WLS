#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "SiliconPlateConstruction.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4Cons.hh"
#include "G4Torus.hh"
#include "G4Hype.hh"

#include "G4Transform3D.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4NistManager.hh"

#include "G4MultiFunctionalDetector.hh"
#include "G4SDManager.hh"
#include "G4PSEnergyDeposit.hh"
#include <G4VPrimitiveScorer.hh>

#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4VoxelLimits.hh"

#include "G4RunManager.hh"
#include "G4PhysicalConstants.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4Navigator.hh"
#include "G4TransportationManager.hh"

#include "G4GDMLParser.hh"

#include <G4VisAttributes.hh>
#include <iostream>
#include <fstream>

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/*
Constructs DetectorConstruction, defines default values.
*/
DetectorConstruction::DetectorConstruction()
 : G4VUserDetectorConstruction(),fPBox(nullptr), fLBox(nullptr),
  fBox(nullptr)
{
  fDetectorMessenger = new DetectorMessenger(this);
  fTargetMPT = new G4MaterialPropertiesTable();
  fExpHall_x = fExpHall_y = fExpHall_z = 0.05*m;
  fTargetName = "holder";
  fThickness = 1*mm;
  fTargetThickness = 3*mm;
  fDetectorType = 0;
  fABSL = 1;
  fRES=4.0;
  fLY=10500./MeV;
  fDetectorName = "6pmt_coverage_pe";
  fVolName = "World";
  fSigAlpha = 0.5;
  DefineMaterials();
//  SetTargetMaterial("Scint");
  SetWorldMaterial("Air");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction(){;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
/*
Sets thickness of target.
*/
void DetectorConstruction::SetSize(G4double value){
  fTargetThickness=value;
  if(fBox){
    fBox->SetZHalfLength(fTargetThickness/2);
  }
  UpdateGeometry();

  G4RunManager::GetRunManager()->PhysicsHasBeenModified();
}

void DetectorConstruction::SetLY(G4double value){
  fLY=value;
  UpdateGeometry();
  G4RunManager::GetRunManager()->PhysicsHasBeenModified();
}

void DetectorConstruction::SetRes(G4double value){
  fRES=value;
  UpdateGeometry();
  G4RunManager::GetRunManager()->PhysicsHasBeenModified();
}

/*
Sets which detector geometry is used.
*/
void DetectorConstruction::SetDetectorType(G4int value){
  fDetectorType=value;

  UpdateGeometry();
  G4RunManager::GetRunManager()->PhysicsHasBeenModified();
}

void DetectorConstruction::SetABS(G4double value){
  fABSL=value;

  UpdateGeometry();
  G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  DefineMaterials();
}

void DetectorConstruction::SetSigAlpha(G4double value){
  fSigAlpha=value;

  UpdateGeometry();
  G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  DefineMaterials();
}

void DetectorConstruction::SetDetectorName(G4String name){
  fDetectorName=name;

  UpdateGeometry();
  G4RunManager::GetRunManager()->PhysicsHasBeenModified();
}

/*
Sets material of target.
*/
void DetectorConstruction::SetTargetMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =
     G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

  if (pttoMaterial) {
    fTargetMaterial = pttoMaterial;
    fTargetName = fTargetMaterial->GetName();
    if ( fLBox ) { fLBox->SetMaterial(fTargetMaterial); }
  } else {
    G4cout << "\n--> warning from DetectorConstruction::SetMaterial : "
           << materialChoice << " not found" << G4endl;
  }
  G4RunManager::GetRunManager()->PhysicsHasBeenModified();
}

void DetectorConstruction::SetRI(G4double value){
  fRI = value;
  UpdateGeometry();
  G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  DefineMaterials();
}

/*
Sets material of world volume.
*/
void DetectorConstruction::SetWorldMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =
     G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

  if (pttoMaterial) {
    fWorldMaterial = pttoMaterial;
    if ( fWLBox ) { fWLBox->SetMaterial(fWorldMaterial); }
  } else {
    G4cout << "\n--> warning from DetectorConstruction::SetMaterial : "
           << materialChoice << " not found" << G4endl;
  }
  G4RunManager::GetRunManager()->PhysicsHasBeenModified();
}

/*
Defines materials used in simulation. Sets material properties for PEN and other optical components.
*/
void DetectorConstruction::DefineMaterials(){// ------------- Materials -------------
  G4double a, z, density;
  G4int nelements;

  // fAir
  //
  G4Element* N = new G4Element("Nitrogen", "N", z=7 , a=14.01*g/mole);
  G4Element* O = new G4Element("Oxygen"  , "O", z=8 , a=16.00*g/mole);

  fAir = new G4Material("Air", density=1.29*mg/cm3, nelements=2);
  fAir->AddElement(N, 70.*perCent);fABSL = 1;
  fAir->AddElement(O, 30.*perCent);

  G4NistManager* man = G4NistManager::Instance();
  // Water
  //
  G4Element* H = new G4Element("Hydrogen", "H", z=1 , a=1.01*g/mole);

  G4Material* water = new G4Material("Water", density= 1.0*g/cm3, nelements=2);
  water->AddElement(H, 2);
  water->AddElement(O, 1);

  G4Element* C = new G4Element("Carbon", "C", z=12, a=12*g/mole);
  G4Element* Pb = new G4Element("Lead", "Pb", z=87, a=207*g/mole);
  fGlass = man->FindOrBuildMaterial("G4_Pyrex_Glass");
  fPOM = new G4Material("POM",density=1.41*g/cm3,nelements=3);
  fPOM->AddElement(O,1);
  fPOM->AddElement(C,1);
  fPOM->AddElement(H,2);

  fABS = new G4Material("ABS",density=1.07*g/cm3,nelements=3);
  fABS->AddElement(C,15);
  fABS->AddElement(H,17);
  fABS->AddElement(N,1);

  // Scintillators
  G4int number_of_atoms;
  fPEN = new G4Material("PEN", density= 1.3*g/cm3, nelements=3);
  fPEN->AddElement(O, number_of_atoms=4);
  fPEN->AddElement(H, number_of_atoms=10);
  fPEN->AddElement(C, number_of_atoms=14);

  G4double wavelength;
  char filler;
  G4double varabsorlength;
  G4double ems;
  G4double rindex;fABSL = 1;

  G4double absEnergy[102]  = {0};
  G4double abs[102]={0};
  G4double emission[102]={0};
  G4double rIndex[102]={0};
  G4double rIndex_fAir[102]={0};
  G4double ems_abs[102]={0};

  G4int absEntries = 0;
  ifstream ReadAbs;

  G4String abs_file = "../input_files/Exp4.csv";
  G4double emission_fibre[102]={0};
  ReadAbs.open(abs_file);
  G4double var = GetABS();
  if(ReadAbs.is_open())
  {
    while(!ReadAbs.eof())
    {
      ReadAbs>>wavelength>>filler>>varabsorlength>>filler>>ems>>filler>>rindex;
      if(ReadAbs.eof()){
        break;
      }
      absEnergy[absEntries] = (1240/wavelength)*eV;
      fABSL = 1;
      abs[absEntries] = varabsorlength*mm;
      emission[absEntries] = ems;
      rIndex[absEntries] = 1.65;
      rIndex_fAir[absEntries]=1.0;
      ems_abs[absEntries]=ems;
      emission_fibre[absEntries]=1.0;
      absEntries++;
    }
  }

  else G4cout<<"Error opening file: " <<abs_file<<G4endl;
  ReadAbs.close();
  absEntries--;

  const G4int nEntries1 = sizeof(absEnergy)/sizeof(G4double);
  assert(sizeof(rIndex) == sizeof(absEnergy));
  assert(sizeof(abs) == sizeof(absEnergy));
  assert(sizeof(emission) == sizeof(absEnergy));
  assert(sizeof(rIndex_fAir == sizeof(absEnergy)));

  G4double Fiber_energy[] = {2.00*eV,2.87*eV,3.2*eV, 8*eV, 9*eV};
  G4double Fiber_abslength[]={9.00*m, 9.00*m, 0.01*mm,0.01*mm, 0.01*mm};

  fTargetMPT->AddProperty("RINDEX",       absEnergy, rIndex, nEntries1)->SetSpline(true);
  fTargetMPT->AddProperty("ABSLENGTH",    absEnergy, abs, nEntries1)->SetSpline(true); // *
  fTargetMPT->AddProperty("FASTCOMPONENT",absEnergy, emission, nEntries1)->SetSpline(true);
  fTargetMPT->AddProperty("SLOWCOMPONENT",absEnergy, emission, nEntries1)->SetSpline(true);

  fTargetMPT->AddProperty("WLSCOMPONENT",absEnergy, emission, nEntries1)->SetSpline(true);
  fTargetMPT->AddProperty("WLSABSLENGTH",    absEnergy, abs, nEntries1)->SetSpline(true);

  fTargetMPT->AddConstProperty("SCINTILLATIONYIELD",10500./MeV); // * 2.5 * PEN = PS, 10*PEN=PS
  fTargetMPT->AddConstProperty("RESOLUTIONSCALE",4.0); // * 1, 4, 8
  fTargetMPT->AddConstProperty("FASTTIMECONSTANT", 5.198*ns);
  fTargetMPT->AddConstProperty("SLOWTIMECONSTANT",24.336*ns);
  fTargetMPT->AddConstProperty("YIELDRATIO",0.05);
  fTargetMPT->AddConstProperty("WLSTIMECONSTANT", 0.5*ns);

  fPEN->SetMaterialPropertiesTable(fTargetMPT);

  density = universe_mean_density;    //from PhysicalConstants.h
  fVacuum = new G4Material("Galactic", z=1., a=1.008*g/mole, density,
                           kStateGas,2.73*kelvin,3.e-18*pascal);
  //
  // fAir
  G4MaterialPropertiesTable* worldMPT = new G4MaterialPropertiesTable();
  worldMPT->AddProperty("RINDEX", absEnergy, rIndex_fAir, nEntries1)->SetSpline(true);

  fAir->SetMaterialPropertiesTable(worldMPT);
  fVacuum->SetMaterialPropertiesTable(worldMPT);
  fSi = man->FindOrBuildMaterial("G4_Si");

  fScintilator = new G4Material("Scint", density= 1.03*g/cm3, 2);
  fScintilator->AddElement(C, 0.475);
  fScintilator->AddElement(H, 0.525);

  Pstyrene = new G4Material("Polystyrene", density= 1.03*g/cm3, 2);
  Pstyrene->AddElement(C, 8);
  Pstyrene->AddElement(H, 8);

  G4double rindexEnergy[500] = {0};
  G4double scintIndex[500] = {0};

  G4int rindexEntries = 0;
  ifstream ReadRindex;

  G4String rindex_file="../input_files/rindexScint.txt";
  ReadRindex.open(rindex_file);

  if(ReadRindex.is_open())
  {
      while(!ReadRindex.eof())
      {
          ReadRindex>>wavelength>>filler>>scintIndex[76-rindexEntries];
          rindexEnergy[76 - rindexEntries] = (1240/wavelength)*eV;
          rindexEntries++;
      }
  }else G4cout<<"Error opening file: "<<rindex_file<<G4endl;
  ReadRindex.close();
  rindexEntries--;

  G4double scintEnergy[501] = {0};
  G4double scintEmit[501] = {0};
  G4double scintEmitSlow[501] = {0};

  G4int scintEntries = 0;
  ifstream ReadScint;

  G4String Scint_file="../input_files/pTP_emission.txt";
  ReadScint.open(Scint_file);

  if(ReadScint.is_open())
  {
      while(!ReadScint.eof())
      {
          ReadScint>>wavelength>>filler>>scintEmit[500-scintEntries];
          //convert wavelength to eV:
          scintEnergy[500 - scintEntries] = (1240/wavelength)*eV;
          scintEmitSlow[500 - scintEntries] = scintEmit[500 - scintEntries];
          scintEntries++;
      }
  }else G4cout<<"Error opening file: "<<Scint_file<<G4endl;
  ReadScint.close();
  scintEntries--;

  G4int absorbEntries = 0;
  G4double varabsorblength;
  G4double absorbEnergy[501] = {0};
  G4double Absorb[501] = {0};

  ifstream ReadAbsorb;
  G4String ReadAbsorbLength="../input_files/PlasticBulkAbsorb2.cfg";

  ReadAbsorb.open(ReadAbsorbLength);
  if (ReadAbsorb.is_open())
  {
      while(!ReadAbsorb.eof())
      {
          ReadAbsorb>>wavelength>>filler>>varabsorblength;
          absorbEnergy[500 - absorbEntries]=(1240/wavelength)*eV;
          Absorb[500 - absorbEntries]=varabsorblength*m;
          absorbEntries++;
      }
  }else G4cout<<"Error opening file: "<<ReadAbsorb<<G4endl;
  ReadAbsorb.close();
  absorbEntries--;

  G4double wlsEnergy[501] = {0};
  G4double wlsEmit[501] = {0};

  G4int wlsScintEntries = 0;
  ifstream ReadWLSScint;

  G4String wls_Scint_file="../input_files/full_popop_emission.cfg";
  ReadWLSScint.open(wls_Scint_file);

  if(ReadWLSScint.is_open())
  {
      while(!ReadWLSScint.eof())
      {
          ReadWLSScint>>wavelength>>filler>>wlsEmit[500-wlsScintEntries];
          //convert wavelength to eV:
          wlsEnergy[500 - wlsScintEntries] = (1240/wavelength)*eV;
          wlsScintEntries++;
      }
  }else G4cout<<"Error opening file: "<<wls_Scint_file<<G4endl;
  ReadWLSScint.close();
  wlsScintEntries--;

  G4int wlsAbsorbEntries = 0;
  G4double wlsAbsorbEnergy[501] = {0};
  G4double wlsAbsorb[501] = {0};

  ifstream ReadWLSAbsorb;
  G4String ReadWLSAbsorbLength="../input_files/scintAbsLen.txt";

  ReadWLSAbsorb.open(ReadWLSAbsorbLength);
  if (ReadWLSAbsorb.is_open())
  {
      while(!ReadWLSAbsorb.eof())
      {
          ReadWLSAbsorb>>wavelength>>filler>>varabsorblength;
          wlsAbsorbEnergy[500 - wlsAbsorbEntries]=(1240/wavelength)*eV;
          wlsAbsorb[500 - wlsAbsorbEntries]=varabsorblength*m;
          wlsAbsorbEntries++;
      }
  }else G4cout<<"Error opening file: "<<ReadWLSAbsorb<<G4endl;
  ReadWLSAbsorb.close();
  wlsAbsorbEntries--;

  G4MaterialPropertiesTable* MPT = new G4MaterialPropertiesTable();

  MPT->AddProperty("WLSABSLENGTH",wlsAbsorbEnergy,wlsAbsorb,wlsAbsorbEntries);
  MPT->AddProperty("WLSCOMPONENT",wlsEnergy,wlsEmit,wlsScintEntries);
  MPT->AddConstProperty("WLSTIMECONSTANT", 12*ns);

  MPT->AddProperty("RINDEX",        rindexEnergy,  scintIndex, rindexEntries);
  MPT->AddProperty("ABSLENGTH",     absorbEnergy, Absorb,     absorbEntries);
  MPT->AddProperty("FASTCOMPONENT", scintEnergy,  scintEmit,  scintEntries);
  MPT->AddProperty("SLOWCOMPONENT",scintEnergy, scintEmitSlow,     scintEntries);

  MPT->AddConstProperty("SCINTILLATIONYIELD",11520./MeV);
  //MPT->AddConstProperty("SCINTILLATIONYIELD",100./MeV);
  MPT->AddConstProperty("RESOLUTIONSCALE",4.0);
  MPT->AddConstProperty("FASTTIMECONSTANT", 2.1*ns);
  MPT->AddConstProperty("SLOWTIMECONSTANT",14.2*ns);
  MPT->AddConstProperty("YIELDRATIO",1.0);

  fScintilator->SetMaterialPropertiesTable(MPT);
}

void DetectorConstruction::SetVolName(G4ThreeVector thePoint){
  G4Navigator* theNavigator = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
  G4VPhysicalVolume* myVolume= theNavigator->LocateGlobalPointAndSetup(thePoint);
  fVolName =  myVolume->GetName();
}

void DetectorConstruction::SetPropertyTable(G4Material* material, G4MaterialPropertiesTable* table){
  material->SetMaterialPropertiesTable(table);
}

void DetectorConstruction::UpdateGeometry(){
  G4RunManager::GetRunManager()->DefineWorldVolume(Construct());
}

/*
Clears stored geometry, then constructs all volumes that can be used in the simulation.

Builds and places volumes in world.

Defines detector sensitivities and properties.
*/
G4VPhysicalVolume* DetectorConstruction::Construct()
{
    G4GDMLParser parser;
    G4GeometryManager::GetInstance()->OpenGeometry();
    G4LogicalVolumeStore::GetInstance()->Clean();
    G4PhysicalVolumeStore::GetInstance()->Clean();
    G4SolidStore::GetInstance()->Clean();

// ------------- Volumes --------------

// The experimental Hall
  fWorldBox = new G4Box("World",fExpHall_x,fExpHall_y,fExpHall_z);

  fWLBox = new G4LogicalVolume(fWorldBox,fWorldMaterial,"World",0,0,0);

  fWPBox = new G4PVPlacement(0,G4ThreeVector(),fWLBox,"World",0,false,0);

  double target_width = 15*mm;
  fTargetThickness = 2.5*mm;
  fBox = new G4Box("target", target_width, target_width, fTargetThickness);
  fLBox = new G4LogicalVolume(fBox,fPEN, "target",0,0,0);
  double position = 0;

  double rod_sides = 2.5*mm;
  double rod_length = 50*mm;

  double curve = 2.5*mm;

  G4Box* rod_stick = new G4Box("target",rod_sides,rod_sides,rod_length);

  G4RotationMatrix* rm = new G4RotationMatrix(0,-90*deg,0);
  G4Tubs* curve_edge = new G4Tubs("curve",0,curve,rod_sides,0,180.*deg);

  G4UnionSolid* rod = new G4UnionSolid("final_rod",rod_stick, curve_edge, rm, G4ThreeVector(0,0,-rod_length));
  G4LogicalVolume* rod_log = new G4LogicalVolume(rod,fPEN,"target",0,0,0);

  G4NistManager* man = G4NistManager::Instance();

  G4ThreeVector point = G4ThreeVector(0,0,5*cm);
  G4Navigator* pointNavigator = new G4Navigator();
  pointNavigator->SetWorldVolume(fWPBox);
  pointNavigator->LocateGlobalPointAndSetup(point);

  // --------------Detectors--------------

  char filler;
  G4double wavelength;
  G4double cath_eff;
  G4double photocath_energy[57];
  G4double photocath_EFF[57];
  G4double perfect_EFF[57];
  G4double perfect_REFL[57];
  G4double photocath_REFL[57]={0};
  G4String pmt_file = "../input_files/pmtQE.csv";

  ifstream ReadEff;
  G4int effCounter = 0;
  ReadEff.open(pmt_file);

  if(ReadEff.is_open())
  {
    while(!ReadEff.eof())
    {
      ReadEff>>wavelength>>filler>>cath_eff;
      if(ReadEff.eof()){
        break;
      }
      photocath_energy[57-effCounter] = (1240/wavelength)*eV;
      photocath_EFF[57-effCounter] = cath_eff;
      perfect_EFF[57-effCounter] = 1;
      perfect_REFL[57-effCounter] = 0;
      effCounter++;
    }
  }

  else G4cout<<"Error opening file: " <<pmt_file<<G4endl;
  ReadEff.close();
  effCounter--;

  const G4int nPMT_EFF = sizeof(photocath_energy)/sizeof(G4double);

  G4OpticalSurface* perfect_optsurf = new G4OpticalSurface("perfect",glisur,polished, dielectric_metal);
  G4MaterialPropertiesTable* detector_MT = new G4MaterialPropertiesTable();
  detector_MT->AddProperty("EFFICIENCY", photocath_energy, perfect_EFF,nPMT_EFF);
  detector_MT->AddProperty("REFLECTIVITY", photocath_energy, perfect_REFL,nPMT_EFF);
  perfect_optsurf->SetMaterialPropertiesTable(detector_MT);

  G4LogicalVolume* tileDetectorLog = new G4LogicalVolume(fBox,fSi, "tile_sensor");
  new G4LogicalSkinSurface("pen_det_surf",tileDetectorLog,perfect_optsurf);

  // Spectrometer Sensor

  G4double spectrometer_side = (35.8/2)*mm;
  G4double spectrometer_top = (23.9/2)*mm;
  G4double spectrometer_thickness = 0.01*mm;

  G4Box* spectrometer_sensor = new G4Box("spec_sensor",spectrometer_thickness,spectrometer_top,spectrometer_side);
  G4LogicalVolume* spec_log = new G4LogicalVolume(spectrometer_sensor,fSi,"spec_sensor_log");
  G4LogicalVolume* abs_log = new G4LogicalVolume(spectrometer_sensor,fSi,"spec_abs_log");

  G4double spec_energy[57];
  G4double spec_EFF[57];
  G4double spec_REFL[57]={0};
  G4String spec_file = "../input_files/specQE.csv";

  effCounter = 0;
  ReadEff.open(spec_file);

  if(ReadEff.is_open())
  {
    while(!ReadEff.eof())
    {
      ReadEff>>wavelength>>filler>>cath_eff;
      if(ReadEff.eof()){
        break;
      }
      spec_energy[57-effCounter] = (1240/wavelength)*eV;
      spec_EFF[57-effCounter] = cath_eff;
      effCounter++;
    }
  }

  else G4cout<<"Error opening file: " <<spec_file<<G4endl;
  ReadEff.close();
  effCounter--;

  const G4int nSPEC_EFF = sizeof(spec_energy)/sizeof(G4double);

  G4OpticalSurface* spec_optsurf = new G4OpticalSurface("photocath_opsurf",glisur,polished, dielectric_metal);
  G4MaterialPropertiesTable* spec_MT = new G4MaterialPropertiesTable();
  spec_MT->AddProperty("EFFICIENCY", spec_energy, spec_EFF,nSPEC_EFF);
  spec_MT->AddProperty("REFLECTIVITY", spec_energy, spec_REFL,nSPEC_EFF);
  spec_optsurf->SetMaterialPropertiesTable(spec_MT);
  new G4LogicalSkinSurface("spec_surf",spec_log,perfect_optsurf);

  G4OpticalSurface* AirPEN = new G4OpticalSurface("AirPEN",glisur, ground, dielectric_dielectric);
  AirPEN -> SetPolish(fSigAlpha);
  AirPEN -> SetMaterialPropertiesTable(fTargetMPT);

  G4VPhysicalVolume* siPM_placement;
  G4VPhysicalVolume* siPM_case_placement;
  G4VPhysicalVolume* siPM_board_placement;

  G4VPhysicalVolume* vacPlacement;

  G4VisAttributes* tileAttr = new G4VisAttributes(G4Colour::Blue());
  tileAttr->SetVisibility(true);
  rod_log->SetVisAttributes(tileAttr);

  G4double casing_length = 30*mm/2;
  G4double casing_height = 32.5*mm/2;
  G4double casing_hole = 29*mm / 2;
  G4double casing_hole_height = 31.5*mm/2;

  G4double photo_cath_length = 26*mm/2;
  G4double photo_cath_eff_length = 23*mm/2;
  G4double photo_cath_height = 0.8*mm;

  G4Box* pmt_case = new G4Box("case",casing_length,casing_length,casing_height);
  G4Box* pHole = new G4Box("hole",casing_hole,casing_hole,casing_hole_height);
  G4SubtractionSolid* pmt_shell = new G4SubtractionSolid("pmt_shell",pmt_case,pHole);

  G4LogicalVolume* pmt_void = new G4LogicalVolume(pHole, fVacuum,"vacuum");
  G4LogicalVolume* pmt_case_log = new G4LogicalVolume(pmt_shell,fPOM,"pmt_case_log");

  G4Box* pmt_cath = new G4Box("pmt_cathode",photo_cath_length,photo_cath_length,photo_cath_height);
  G4Box* pmt_active = new G4Box("pmt_active",photo_cath_eff_length,photo_cath_eff_length,photo_cath_height);
  G4SubtractionSolid* pmt_inactive_cath = new G4SubtractionSolid("pmt_inactive_cathode",pmt_cath,pmt_active);
  G4LogicalVolume* pmt_inactive_cath_log = new G4LogicalVolume(pmt_inactive_cath,fGlass,"pmt_inactive_cath_log");
  G4LogicalVolume* pmt_cath_log = new G4LogicalVolume(pmt_active,fGlass,"pmt_cath_log");

  new G4LogicalSkinSurface("spec_surf",pmt_cath_log,perfect_optsurf);

  G4RotationMatrix* rotationMatrix = new G4RotationMatrix(0,0,0);
  G4RotationMatrix* rotationMatrix1 = new G4RotationMatrix(0,0,0);
  G4RotationMatrix* rotationMatrix2 = new G4RotationMatrix(0,0,0);

  G4VPhysicalVolume* pmtPlacement;
  G4VPhysicalVolume* incathPlacement;
  G4VPhysicalVolume* cathPlacement;

  /*
  0 - PMT on base of tile, collimator included.
  */
  fDetectorType = 0;

  fPBox = new G4PVPlacement(0, G4ThreeVector(0,0,0),fLBox,"target",fWLBox,false,0,false);
  cathPlacement = new G4PVPlacement(rotationMatrix1,G4ThreeVector(0,0,-(fTargetThickness+photo_cath_height)),pmt_cath_log,"active_detector1",fWLBox,false,0,true);
  incathPlacement= new G4PVPlacement(rotationMatrix1,G4ThreeVector(0,0,-(fTargetThickness+photo_cath_height)),pmt_inactive_cath_log,"inactive_detector1",fWLBox,false,0,true);
  //pmtPlacement = new G4PVPlacement(rotationMatrix1,G4ThreeVector(0,0,-(casing_height+photo_cath_height)),pmt_case_log,"pmt1", pmt_cath_log,false,0,true);

   G4LogicalBorderSurface* surfaceAirPEN = new G4LogicalBorderSurface("AirPEN",fWPBox,fPBox,AirPEN);
   G4LogicalBorderSurface* surfacePENAir = new G4LogicalBorderSurface("AirPEN",fPBox,fWPBox,AirPEN);

  return fWPBox;
}
