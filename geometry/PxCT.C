

/********************************************************************************
 *    Copyright (C) 2014 GSI Helmholtzzentrum fuer Schwerionenforschung GmbH    *
 *                                                                              *
 *              This software is distributed under the terms of the             *
 *         GNU Lesser General Public Licence version 3 (LGPL) version 3,        *
 *                  copied verbatim in the file "LICENSE"                       *
 ********************************************************************************/
// Note: In root all sizes are given in cm

#include "TFile.h"
#include "TGeoCompositeShape.h"
#include "TGeoManager.h"
#include "TGeoMaterial.h"
#include "TGeoMatrix.h"
#include "TGeoMedium.h"
#include "TGeoPgon.h"
#include "TGeoVolume.h"
#include "TList.h"
#include "TMath.h"
#include "TROOT.h"
#include "TString.h"
#include "TSystem.h"

#include <iostream>

// Name of geometry version and output file
const TString geoVersion = "PxCT";
const TString FileName = geoVersion + ".root";
const TString FileName1 = geoVersion + "_geomanager.root";

// Names of the different materials which are used to build the modules
// The materials are defined in the global media.geo file
const TString MediumGas = "GADGET_P10_800";

const TString MediumVacuum = "vacuum4";

const TString PxCTMatter = "germanium";
const TString WindowMatter = "carbon";
// tpc

const TString CylinderVolumeMedium = "steel";

const TString WindowMedium = "kapton";
const TString MicromegasMedium = "G10";

//PEEK
const TString PEEKMedium = "PEEK";

// Detector Dimensions (cm)
const Float_t tpc_diameter_in = 15.4051;
const Float_t tpc_diameter_out = 16.51;
const Float_t drift_length = 40.;

const Float_t upstream_ring_dia_in = 11.684;
const Float_t upstream_ring_dia_out = tpc_diameter_out;
const Float_t upstream_ring_length = 2.6162 + 0.9398;
const Float_t downstream_ring_dia_in = 15.4051;
const Float_t downstream_ring_dia_out = 22.86;
const Float_t downstream_ring_length = 1.5748;
const Float_t upstream_cap_dia_in = 5.08;
const Float_t upstream_cap_dia_out = tpc_diameter_out;
const Float_t upstream_cap_length = 1.7272;
const Float_t downstream_cap_dia_in = 18.161;
const Float_t downstream_cap_dia_out = 22.86;
const Float_t downstream_cap_length = 2.6162;
const Float_t window_dia_in = 5.08;
const Float_t window_dia_out = 8.4074;
const Float_t window_length = 0.8636;
const Float_t micromegas_dia_in = 0;
const Float_t micromegas_dia_out = 22.86;
const Float_t micromegas_length = 0.2794;
const Float_t micromegas_ring_dia_in = 0;
const Float_t micromegas_ring_dia_out = 22.86;
const Float_t micromegas_ring_length = 0.635;
const Float_t cathode_dia_in = 5.5118;
const Float_t cathode_dia_out = 12.1844;
const Float_t cathode_length = 1.1016;
const Float_t cathode_mount_dia_in = 6.096;
const Float_t cathode_mount_dia_out = 14.8844;
const Float_t cathode_mount_length = 0.9525;
const Float_t cage_dia_in = 12.2074;
const Float_t cage_dia_out = (12.2074 + (12.2174 - 12.2074) * 0.8);
const Float_t cage_length = 38.0956;
const Float_t cage_support_shift = 2.71655;
const Float_t cage_support_dia_in = 12.2174;
const Float_t cage_support_dia_out = 13.8684;
const Float_t cage_support_length = 38.6867;
//PxCT
const Float_t QSD1DetectorToSourceDistance = 8.5;
const Float_t QSD1DetectorThickness = 6.52;
const Float_t QSD1DetectorDiameter = 8.48;
const Float_t QSD1WindowThickness = 0.060;
const Float_t QSD1DetectorToWindowDistance = 0.68;

const Float_t QSD2DetectorToSourceDistance = 8.5;
const Float_t QSD2DetectorThickness = 8.0;
const Float_t QSD2DetectorDiameter = 7.98;
const Float_t QSD2WindowThickness = 0.060;
const Float_t QSD2DetectorToWindowDistance = 0.63;
const Float_t zloc= drift_length/2;

//PEEK
const Float_t fPDetCageSupportInnRad = 12.2174/2.0;
const Float_t fPDetCageSupportOutRad = 13.8684/2.0;
const Float_t fPDetCageSupportHalfLen = 37.6867/2.0;
const Float_t fPDetCenter = drift_length/2; 
const Float_t CutOutWindowAPos =  (31.4427-7.2441)/2.0;//z positions
const Float_t CutOutWindowBPos =  (7.3711-31.3157)/2.0;//z positions
const Float_t CutOutWindowCPos =  (19.083-19.6037)/2.0;//z positions
// some global variables
TGeoManager *gGeoMan = new TGeoManager("ATTPC", "ATTPC");
;                     // Pointer to TGeoManager instance
TGeoVolume *gModules; // Global storage for module types

// Forward declarations
void create_materials_from_media_file();
TGeoVolume *create_detector();
void position_detector();
void add_alignable_volumes();

void PxCT()
{
   // Load the necessary FairRoot libraries
   gSystem->Load("libGeoBase");
   gSystem->Load("libParBase");
   gSystem->Load("libBase");

   // Load needed material definition from media.geo file
   create_materials_from_media_file();

   // Get the GeoManager for later usage
   gGeoMan = (TGeoManager *)gROOT->FindObject("FAIRGeom");
   gGeoMan->SetVisLevel(7);

   // Create the top volume
   TGeoVolume *top = new TGeoVolumeAssembly("TOP");
   gGeoMan->SetTopVolume(top);

   TGeoMedium *gas = gGeoMan->GetMedium(MediumVacuum);
   TGeoVolume *tpcvac = new TGeoVolumeAssembly(geoVersion);
   tpcvac->SetMedium(gas);
   top->AddNode(tpcvac, 1);

   gModules = create_detector();

   cout << "Voxelizing." << endl;
   top->Voxelize("");
   gGeoMan->CloseGeometry();

   gGeoMan->CheckOverlaps(0.001);
   gGeoMan->PrintOverlaps();
   gGeoMan->Test();

   TFile *outfile = new TFile(FileName, "RECREATE");
   top->Write();
   outfile->Close();

   TFile *outfile1 = new TFile(FileName1, "RECREATE");
   gGeoMan->Write();
   outfile1->Close();

   top->Draw("ogl");
}

void create_materials_from_media_file()
{

   // Use the FairRoot geometry interface to load the media which are already defined
   FairGeoLoader *geoLoad = new FairGeoLoader("TGeo", "FairGeoLoader");
   FairGeoInterface *geoFace = geoLoad->getGeoInterface();
   TString geoPath = gSystem->Getenv("VMCWORKDIR");
   TString geoFile = geoPath + "/geometry/media.geo";
   geoFace->setMediaFile(geoFile);
   geoFace->readMedia();

   // Read the required media and create them in the GeoManager
   FairGeoMedia *geoMedia = geoFace->getMedia();
   FairGeoBuilder *geoBuild = geoLoad->getGeoBuilder();

   FairGeoMedium *vacuum4 = geoMedia->getMedium("vacuum4");
   FairGeoMedium *GADGET_5IsoAr_800 = geoMedia->getMedium("GADGET_5IsoAr_800");
   FairGeoMedium *steel = geoMedia->getMedium("steel");

   FairGeoMedium *kapton = geoMedia->getMedium("kapton");
   FairGeoMedium *G10 = geoMedia->getMedium("G10");
   FairGeoMedium *germanium = geoMedia->getMedium("germanium");
   FairGeoMedium *carbon = geoMedia->getMedium("carbon");

   FairGeoMedium *PEEK = geoMedia->getMedium("PEEK");

   // include check if all media are found
   geoBuild->createMedium(GADGET_5IsoAr_800);
   geoBuild->createMedium(steel);

   geoBuild->createMedium(kapton);
   geoBuild->createMedium(G10);
   ;
   geoBuild->createMedium(vacuum4);

   geoBuild->createMedium(germanium);
   geoBuild->createMedium(carbon);

   geoBuild->createMedium(PEEK);
}

TGeoVolume *create_detector()
{

   // needed materials
   TGeoMedium *OuterCylinder = gGeoMan->GetMedium(CylinderVolumeMedium);
   TGeoMedium *gas = gGeoMan->GetMedium(MediumGas);
   TGeoMedium *windowmat = gGeoMan->GetMedium(WindowMedium);
   TGeoMedium *micromegasmat = gGeoMan->GetMedium(MicromegasMedium);
   TGeoMedium *mediumvacuum4 = gGeoMan->GetMedium(MediumVacuum);
   TGeoMedium *pxctmatter = gGeoMan->GetMedium(PxCTMatter);
   TGeoMedium *windowmatter = gGeoMan->GetMedium(WindowMatter);
   TGeoMedium *PEEKmatter = gGeoMan->GetMedium(PEEKMedium);

   TGeoVolume **Cry_vol;
   Cry_vol = new TGeoVolume *[3];

   TString CrystalName = "Crystal_";
   TString name_cry[3] = {"01", "02", "03"};

//PxCT
      Cry_vol[0] = gGeoManager->MakeTube(CrystalName + name_cry[0], pxctmatter, 0,QSD1DetectorDiameter/2,QSD1DetectorThickness/2);
      Cry_vol[0]->SetLineColor(kCyan);
      gGeoMan->GetVolume(geoVersion)
         ->AddNode(Cry_vol[0],  1,
                   new TGeoCombiTrans(tpc_diameter_out / 2+ QSD1WindowThickness+QSD1DetectorThickness/2,0,zloc,
                                      new TGeoRotation(CrystalName + name_cry[0],90,90, 0)));
      Cry_vol[0]->SetTransparency(90);
    
      TGeoVolume *window1 = gGeoManager->MakeTube("window1", windowmatter, 0, QSD1DetectorDiameter/2, QSD1WindowThickness);
      window1->SetLineColor(kRed);
      gGeoMan->GetVolume(geoVersion)
         ->AddNode(window1, 1,
                   new TGeoCombiTrans(tpc_diameter_out / 2+ QSD1WindowThickness,0,zloc,
                                      new TGeoRotation("window1", 90,90, 0.0)));
      window1->SetTransparency(90);

      Cry_vol[1] = gGeoManager->MakeTube(CrystalName + name_cry[1], pxctmatter, 0,QSD2DetectorDiameter/2,QSD2DetectorThickness/2);
      Cry_vol[1]->SetLineColor(kCyan);
      gGeoMan->GetVolume(geoVersion)
         ->AddNode(Cry_vol[1],  1,
                   new TGeoCombiTrans(-(tpc_diameter_out / 2+ QSD2WindowThickness+QSD2DetectorThickness/2),0,zloc,
                                      new TGeoRotation(CrystalName + name_cry[1],-90,-90, 0)));

      Cry_vol[1]->SetTransparency(90);

      TGeoVolume *window2 = gGeoManager->MakeTube("window2", windowmatter, 0, QSD2DetectorDiameter/2, QSD2WindowThickness);
      window2->SetLineColor(kRed);
      gGeoMan->GetVolume(geoVersion)
         ->AddNode(window2, 1,
                   new TGeoCombiTrans(-(tpc_diameter_out / 2+ QSD2WindowThickness),0,zloc,
                                      new TGeoRotation("window2", -90, -90, 0.0)));
      window2->SetTransparency(90);


   //PEEK
   TGeoTube *PEEK_vol = new TGeoTube("PEEK_vol", fPDetCageSupportInnRad, fPDetCageSupportOutRad, fPDetCageSupportHalfLen);

   TGeoTube *HoleA = new TGeoTubeSeg("HoleA", fPDetCageSupportInnRad-0.1, fPDetCageSupportOutRad+0.1, 13.4874/2.0,(-90-65.4408867/2.0),(-90-65.4408867/2.0)+65.4408867);//
   TGeoCombiTrans *HoleAtrans = new TGeoCombiTrans("HoleAtrans", 0, 0, CutOutWindowAPos, new TGeoRotation("HoleAtrans", 0, 0, 0));
   HoleAtrans->RegisterYourself();

   TGeoTube *HoleB = new TGeoTubeSeg("HoleB", fPDetCageSupportInnRad-0.1, fPDetCageSupportOutRad+0.1, 13.4874/2.0,(-90-65.4408867/2.0+65.4408867+54.5587858),(-90-65.4408867/2.0+65.4408867+54.5587858)+65.4408867);
   TGeoCombiTrans *HoleBtrans = new TGeoCombiTrans("HoleBtrans", 0, 0, CutOutWindowBPos, new TGeoRotation("HoleBtrans", 0, 0, 0));
   HoleBtrans->RegisterYourself();

   TGeoTube *HoleC = new TGeoTubeSeg("HoleC", fPDetCageSupportInnRad-0.1, fPDetCageSupportOutRad+0.1, 13.4874/2.0,(-90-65.4408867/2.0+65.4408867+54.5587858+65.4408867+54.5387358),(-90-65.4408867/2.0+65.4408867+54.5587858+65.4408867+54.5387358)+65.4409);
   TGeoCombiTrans *HoleCtrans = new TGeoCombiTrans("HoleCtrans", 0, 0, CutOutWindowCPos, new TGeoRotation("HoleCtrans", 0, 0, 0));
   HoleCtrans->RegisterYourself();

TGeoCompositeShape *PeekCyl = new TGeoCompositeShape("PeekCyl", "(PEEK_vol - HoleA:HoleAtrans -HoleB:HoleBtrans- HoleC:HoleCtrans)");//- HoleB:HoleBtrans- HoleC:HoleCtrans

Cry_vol[2] = new TGeoVolume(CrystalName + name_cry[2], PeekCyl, pxctmatter);
        Cry_vol[2]->SetLineColor(kMagenta);
        gGeoMan->GetVolume(geoVersion)->AddNode(Cry_vol[2], 0, new TGeoCombiTrans(0,0,fPDetCenter, new TGeoRotation(CrystalName + name_cry[2], 0, 0, 0))); 
        
        Cry_vol[2]->SetTransparency(0); 
        
   // GADGET Main drift volume with peek
   TGeoTube *drift_volumes = new TGeoTube("drift_volumes", tpc_diameter_in / 2, tpc_diameter_out / 2, drift_length / 2);
   TGeoCompositeShape *drift_vol = new TGeoCompositeShape("drift_vol", "(drift_volumes-(PEEK_vol - HoleA:HoleAtrans -HoleB:HoleBtrans- HoleC:HoleCtrans))");
   
   double tpc_rot = 0;
   TGeoVolume *drift_volume = new TGeoVolume("drift_volume", drift_vol,gas);
//TGeoVolume *drift_volume = gGeoManager->MakeTube("drift_volume", gas, 0, tpc_diameter_in / 2, drift_length / 2);
   gGeoMan->GetVolume(geoVersion)
      ->AddNode(drift_volume, 1,
                new TGeoCombiTrans(0.0, 0, drift_length / 2.0, new TGeoRotation("drift_volume", 0, tpc_rot, 0)));
   drift_volume->SetTransparency(90);

   // GADGET Outer steel chamber
   TGeoVolume *steel_chamber = gGeoManager->MakeTube("steel_chamber", OuterCylinder, tpc_diameter_in / 2,
                                                     tpc_diameter_out / 2, drift_length / 2);
   gGeoMan->GetVolume(geoVersion)
      ->AddNode(steel_chamber, 1,
                new TGeoCombiTrans(0.0, 0, drift_length / 2.0, new TGeoRotation("steel_chamber", 0, tpc_rot, 0)));
   steel_chamber->SetTransparency(90);

   // GADGET Upstream  Ring
   TGeoVolume *upstream_ring = gGeoManager->MakeTube("upstream_ring", OuterCylinder, upstream_ring_dia_in / 2,
                                                     upstream_ring_dia_out / 2, upstream_ring_length / 2);
   gGeoMan->GetVolume(geoVersion)
      ->AddNode(
         upstream_ring, 1,
         new TGeoCombiTrans(0.0, 0, -upstream_ring_length / 2.0, new TGeoRotation("upstream_ring", 0, tpc_rot, 0)));
   upstream_ring->SetTransparency(90);

   // GADGET Upstream Cap
   TGeoVolume *upstream_cap = gGeoManager->MakeTube("upstream_cap", OuterCylinder, upstream_cap_dia_in / 2,
                                                    upstream_cap_dia_out / 2, upstream_cap_length / 2);
   gGeoMan->GetVolume(geoVersion)
      ->AddNode(upstream_cap, 1,
                new TGeoCombiTrans(0.0, 0, -upstream_cap_length / 2.0 - upstream_ring_length,
                                   new TGeoRotation("upstream_cap", 0, tpc_rot, 0)));
   upstream_cap->SetTransparency(90);

   // GADGET Downstream  Ring
   TGeoVolume *downstream_ring = gGeoManager->MakeTube("downstream_ring", OuterCylinder, downstream_ring_dia_in / 2,
                                                       downstream_ring_dia_out / 2, downstream_ring_length / 2);
   gGeoMan->GetVolume(geoVersion)
      ->AddNode(downstream_ring, 1,
                new TGeoCombiTrans(0.0, 0, drift_length + downstream_ring_length / 2.0,
                                   new TGeoRotation("downstream_ring", 0, tpc_rot, 0)));
   downstream_ring->SetTransparency(90);

   // GADGET Downstream Cap
   TGeoVolume *downstream_cap = gGeoManager->MakeTube("downstream_cap", OuterCylinder, downstream_cap_dia_in / 2,
                                                      downstream_cap_dia_out / 2, downstream_cap_length / 2);
   gGeoMan->GetVolume(geoVersion)
      ->AddNode(downstream_cap, 1,
                new TGeoCombiTrans(0.0, 0, downstream_cap_length / 2.0 + downstream_ring_length + drift_length,
                                   new TGeoRotation("downstream_cap", 0, tpc_rot, 0)));
   downstream_cap->SetTransparency(90);

   // GADGET Entrance Window
   TGeoVolume *tpc_window =
      gGeoManager->MakeTube("tpc_window", OuterCylinder, window_dia_in / 2, window_dia_out / 2, window_length / 2);
   gGeoMan->GetVolume(geoVersion)
      ->AddNode(tpc_window, 1,
                new TGeoCombiTrans(0.0, 0, -upstream_ring_length - upstream_cap_length - window_length / 2,
                                   new TGeoRotation("tpc_window", 0, tpc_rot, 0)));
   tpc_window->SetTransparency(90);

   // GADGET Kapton Window
   TGeoVolume *kapton_window = gGeoManager->MakeTube("kapton_window", windowmat, 0, 5.08 / 2.0, 0.0127 / 2.0);
   kapton_window->SetLineColor(kBlue);
   gGeoMan->GetVolume(geoVersion)
      ->AddNode(kapton_window, 1,
                new TGeoCombiTrans(0.0, 0.0, -upstream_ring_length - upstream_cap_length,
                                   new TGeoRotation("kapton_window", 0, tpc_rot, 0)));
   kapton_window->SetTransparency(90);

   // GADGET Micromeags
   TGeoVolume *micromegas =
      gGeoManager->MakeTube("micromegas", micromegasmat, 0, micromegas_dia_out / 2, micromegas_length / 2);
   micromegas->SetLineColor(kRed);
   gGeoMan->GetVolume(geoVersion)
      ->AddNode(micromegas, 1,
                new TGeoCombiTrans(
                   0.0, 0.0, drift_length + downstream_ring_length + downstream_cap_length + micromegas_length / 2,
                   new TGeoRotation("micromegas", 0, tpc_rot, 0)));
   micromegas->SetTransparency(90);

   // GADGET Micromegas Ring
   TGeoVolume *micromegas_ring = gGeoManager->MakeTube("micromegas_ring", OuterCylinder, micromegas_ring_dia_in / 2,
                                                       micromegas_ring_dia_out / 2, micromegas_ring_length / 2);
   gGeoMan->GetVolume(geoVersion)
      ->AddNode(micromegas_ring, 1,
                new TGeoCombiTrans(0.0, 0,
                                   micromegas_ring_length / 2.0 + drift_length + micromegas_length +
                                      downstream_ring_length + downstream_cap_length,
                                   new TGeoRotation("micromegas_ring", 0, tpc_rot, 0)));
   micromegas_ring->SetTransparency(90);

   // GADGET Cathode
   TGeoVolume *cathode =
      gGeoManager->MakeTube("cathode", OuterCylinder, cathode_dia_in / 2, cathode_dia_out / 2, cathode_length / 2);
   cathode->SetLineColor(kGreen);
   gGeoMan->GetVolume(geoVersion)
      ->AddNode(cathode, 1, new TGeoCombiTrans(0.0, 0, cathode_length / 2, new TGeoRotation("cathode", 0, tpc_rot, 0)));
   downstream_cap->SetTransparency(90);

   // GADGET Cathode Mount
   TGeoVolume *cathode_mount = gGeoManager->MakeTube("cathode_mount", OuterCylinder, cathode_mount_dia_in / 2,
                                                     cathode_mount_dia_out / 2, cathode_mount_length / 2);
   gGeoMan->GetVolume(geoVersion)
      ->AddNode(cathode_mount, 1,
                new TGeoCombiTrans(0.0, 0, cathode_mount_length / 2, new TGeoRotation("cathode_mount", 0, tpc_rot, 0)));
   downstream_cap->SetTransparency(90);
   
   return Cry_vol[0];
}
