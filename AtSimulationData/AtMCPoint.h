/********************************************************************************
 *    Copyright (C) 2014 GSI Helmholtzzentrum fuer Schwerionenforschung GmbH    *
 *                                                                              *
 *              This software is distributed under the terms of the             *
 *         GNU Lesser General Public Licence version 3 (LGPL) version 3,        *
 *                  copied verbatim in the file "LICENSE"                       *
 ********************************************************************************/
#ifndef ATMCPOINT_H
#define ATMCPOINT_H 1

#include <FairMCPoint.h>

#include <Math/Point3D.h>
#include <Math/Point3Dfwd.h> // for XYZPoint
#include <Math/Vector3D.h>
#include <Math/Vector3Dfwd.h> // for XYZVector
#include <Rtypes.h>
#include <TString.h>
#include <TVector3.h>

class TBuffer;
class TClass;
class TMemberInspector;
class AtSimpleSimulation;

class AtMCPoint : public FairMCPoint {

protected:
   using XYZPoint = ROOT::Math::XYZPoint;
   using XYZVector = ROOT::Math::XYZVector;

   Int_t fDetCopyID = 0;
   TString fVolName;
   Double_t fEnergyIni = 0;
   Double_t fAngleIni = 0;
   Int_t fAiso = 0;
   Int_t fZiso = 0;

public:
   /** Default constructor **/
   AtMCPoint();

   /** Constructor with arguments
    *@param trackID  Index of MCTrack
    *@param detID    Detector ID
    *@param pos      Ccoordinates at entrance to active volume [cm]
    *@param mom      Momentum of track at entrance [GeV]
    *@param tof      Time since event start [ns]
    *@param length   Track length since creation [cm]
    *@param eLoss    Energy deposit [GeV]
    **/
   AtMCPoint(Int_t trackID, Int_t detID, TVector3 pos, TVector3 mom, Double_t tof, Double_t length, Double_t eLoss);
   AtMCPoint(Int_t trackID, Int_t detID, XYZPoint pos, XYZVector mom, Double_t tof, Double_t length, Double_t eLoss);

   AtMCPoint(Int_t trackID, Int_t detID, TVector3 pos, TVector3 mom, Double_t tof, Double_t length, Double_t eLoss,
             TString VolName, Int_t detCopyID, Double_t EIni, Double_t AIni, Int_t A, Int_t Z);

   /** Destructor **/
   virtual ~AtMCPoint() = default;
   /** Copy constructor **/
   AtMCPoint(const AtMCPoint &point) = delete;
   AtMCPoint operator=(const AtMCPoint &point) = delete;

   virtual void Clear(Option_t *) override;

   /** Accessors **/
   Int_t GetDetCopyID() const { return fDetCopyID; } // added by Marc
   TString GetVolName() const { return fVolName; }
   Double_t GetEIni() const { return fEnergyIni; }
   Double_t GetAIni() const { return fAngleIni; }
   Int_t GetMassNum() const { return fAiso; }
   Int_t GetAtomicNum() const { return fZiso; }

   void SetDetCopyID(Int_t id) { fDetCopyID = id; };         // added by Marc
   void SetVolName(TString VolName) { fVolName = VolName; }; // added by Ari
   void SetPosition(const XYZPoint &pos) { SetXYZ(pos.X(), pos.Y(), pos.Z()); }
   void SetMomentum(const XYZVector &mom)
   {
      TVector3 mom2(mom.X(), mom.Y(), mom.Z());
      FairMCPoint::SetMomentum(mom2);
   }

   /** Output to screen **/
   virtual void Print(const Option_t *opt) const override;

   friend AtSimpleSimulation;
   ClassDefOverride(AtMCPoint, 2)
};

#endif
