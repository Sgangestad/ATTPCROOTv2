void DeGAi_60Co(Int_t nEvents = 50000, TString mcEngine = "TGeant4")
{

   TString dir = getenv("VMCWORKDIR");

   // Output file name
   TString outFile = "./DeGAi_60Co.root";

   // Parameter file name
   TString parFile = "./DeGAipar.root";

   // -----   Timer   --------------------------------------------------------
   TStopwatch timer;
   timer.Start();
   // ------------------------------------------------------------------------

   // -----   Create simulation run   ----------------------------------------
   FairRunSim *run = new FairRunSim();
   run->SetName(mcEngine);      // Transport engine
   run->SetOutputFile(outFile); // Output file
   FairRuntimeDb *rtdb = run->GetRuntimeDb();
   // ------------------------------------------------------------------------

   // -----   Create media   -------------------------------------------------
   run->SetMaterials("media.geo"); // Materials
   // ------------------------------------------------------------------------

   // -----   Create geometry   ----------------------------------------------

   FairModule *cave = new AtCave("CAVE");
   cave->SetGeometryFileName("cave.geo");
   run->AddModule(cave);

   // FairModule* magnet = new AtMagnet("Magnet");
   // run->AddModule(magnet);

   /*FairModule* pipe = new AtPipe("Pipe");
   run->AddModule(pipe);*/

   FairDetector *DeGAi = new AtDeGAi("AtDeGAi", kTRUE);
   DeGAi->SetGeometryFileName("DeGAi.root");
   // ATTPC->SetModifyGeometry(kTRUE);
   run->AddModule(DeGAi);

   // ------------------------------------------------------------------------

   // -----   Create PrimaryGenerator   --------------------------------------
   FairPrimaryGenerator *primGen = new FairPrimaryGenerator();

   Double_t pdgId = 22;       // 22 for gamma emission, 2212 for proton emission
     Double_t theta1 = 0;      // polar angle distribution: lower edge (50)
     Double_t theta2 = 180.;    // polar angle distribution: upper edge (51)
     Double_t momentum = 0.005; // GeV/c
     Int_t multiplicity = 1;
     AtTPCGammaDummyGenerator* gammasGen = new AtTPCGammaDummyGenerator(pdgId, multiplicity);
     gammasGen->SetThetaRange(theta1, theta2);
     gammasGen->SetCosTheta();
     gammasGen->SetPRange(momentum, momentum);
	gammasGen->SetNuclearDecayChain();
     gammasGen->SetDecayChainPoint(0.001173240,0.499337);
     /*gammasGen->SetDecayChainPoint(0.0008261,0.000076);
    gammasGen->SetDecayChainPoint(0.00034714,0.000075);*/
     gammasGen->SetDecayChainPoint(0.001332508,0.500663);
    /* gammasGen->SetDecayChainPoint(0.00215861,0.000012);
     gammasGen->SetDecayChainPoint(0.002505748,0.00000002);*/
     gammasGen->SetPhiRange(0., 360.); //(2.5,4)
    gammasGen->SetXYZ(0.0,0,15);
     gammasGen->SetLorentzBoost(0.0); // for instance beta=0.8197505718204776 for 700 A MeV
     // add the gamma generator
     primGen->AddGenerator(gammasGen);


     run->SetGenerator(primGen);

   // ------------------------------------------------------------------------

   //---Store the visualiztion info of the tracks, this make the output file very large!!
   //--- Use it only to display but not for production!
   run->SetStoreTraj(kTRUE);

   // -----   Initialize simulation run   ------------------------------------
   run->Init();
   // ------------------------------------------------------------------------

   // -----   Runtime database   ---------------------------------------------

   Bool_t kParameterMerged = kTRUE;
   FairParRootFileIo *parOut = new FairParRootFileIo(kParameterMerged);
   parOut->open(parFile.Data());
   rtdb->setOutput(parOut);
   rtdb->saveOutput();
   rtdb->print();
   // ------------------------------------------------------------------------

   // -----   Start run   ----------------------------------------------------
   run->Run(nEvents);

   // You can export your ROOT geometry ot a separate file
   run->CreateGeometryFile("./geofile_gamma_full.root");
   // ------------------------------------------------------------------------

   // -----   Finish   -------------------------------------------------------
   timer.Stop();
   Double_t rtime = timer.RealTime();
   Double_t ctime = timer.CpuTime();
   cout << endl << endl;
   cout << "Macro finished succesfully." << endl;
   cout << "Output file is " << outFile << endl;
   cout << "Parameter file is " << parFile << endl;
   cout << "Real time " << rtime << " s, CPU time " << ctime << "s" << endl << endl;
   // ------------------------------------------------------------------------
}
