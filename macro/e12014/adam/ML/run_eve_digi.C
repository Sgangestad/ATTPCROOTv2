/*#include "TString.h"
#include "AtEventDrawTask.h"
#include "AtEventManager.h"

#include "FairLogger.h"
#include "FairParRootFileIo.h"
#include "FairRunAna.h"
*/

#include "RadiusCheck.h"

void run_eve_digi(TString OutputDataFile = "./data/output.sim_display.root")
{
   TString InputDataFile = "/home/faculty/aanthony/fission/data/e12014/unpacked/Bi200Sim.root";
   InputDataFile = "/home/faculty/aanthony/attpcroot/macro/e12014/adam/ML/data/output_digi01.root";
   std::cout << "Opening: " << InputDataFile << std::endl;

   TString dir = getenv("VMCWORKDIR");
   TString geoFile = "ATTPC_v1.1_geomanager.root";
   TString mapFile = "e12014_pad_mapping.xml";

   TString InputDataPath = InputDataFile;
   TString OutputDataPath = OutputDataFile;
   TString GeoDataPath = dir + "/geometry/" + geoFile;
   TString mapDir = dir + "/scripts/" + mapFile;

   FairRunAna *fRun = new FairRunAna();
   FairRootFileSink *sink = new FairRootFileSink(OutputDataFile);
   FairFileSource *source = new FairFileSource(InputDataFile);
   fRun->SetSource(source);
   fRun->SetSink(sink);
   fRun->SetGeomFile(GeoDataPath);

   FairRuntimeDb *rtdb = fRun->GetRuntimeDb();
   FairParRootFileIo *parIo1 = new FairParRootFileIo();
   // parIo1->open("param.dummy.root");
   rtdb->setFirstInput(parIo1);

   auto fMap = std::make_shared<AtTpcMap>();
   fMap->ParseXMLMap(mapDir.Data());
   AtViewerManager *eveMan = new AtViewerManager(fMap);

     

   auto tabMain = std::make_unique<AtTabMain>();
   tabMain->SetMultiHit(100); // Set the maximum number of multihits in the visualization
   eveMan->AddTab(std::move(tabMain));

   auto tabPad = std::make_unique<AtTabPad>(2, 2);
   tabPad->DrawRawADC(0, 0);
   tabPad->DrawADC(0, 1);
   tabPad->DrawArrayAug("Q", 1, 0);
   tabPad->DrawADC(1, 1);
   eveMan->AddTab(std::move(tabPad));

   
   eveMan->Init();

   std::cout << "Finished init" << std::endl;
   // eveMan->RunEvent(27);
}
