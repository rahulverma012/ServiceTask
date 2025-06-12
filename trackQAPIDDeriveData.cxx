// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
///
/// \brief  Service Work Code for Alice Data
/// \author Rahul Verma (rahul.verma@iitb.ac.in) :: Marian I Ivanov (mivanov@cern.ch)

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"

#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/FT0Corrected.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/EventSelection.h"

#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/Core/HfHelper.h"

#include "DerivedPIDTables.h"
#include "Common/DataModel/OccupancyTables.h"

#include <typeinfo>
#include <TRandom.h>
#include <chrono>

#include "Framework/AnalysisDataModel.h"
#include "Framework/HistogramRegistry.h"
#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsFT0/Digit.h"
#include "DataFormatsParameters/GRPLHCIFData.h"

#include "CommonConstants/PhysicsConstants.h"
#include "Common/Core/RecoDecay.h"
#include "TPDGCode.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;
using namespace o2::constants::math;



// #include "MetadataHelper.h"

//#00 Make yash's tasks work --> It is not terminating when addded with this task
//#01 Add Phi to the tasks --> Done
//#02 Index matching macros for pointer validaitons ==> Becomes a time consuming part --> choose another strategy
//#03 Creating big data set on GSI for validation and working check
//#04 Add hPt for phi as well -->Done
//#05 Add Igor's Occupancy definitions --> It was already there
//#06 Add df counter ==> Done 
//#07 Add Histograms entries ==>Done
//#08 Add centrality ==> Done 
//#09 Invariant mass histogram for QA debug ==> Done
//#10 What is the interpolation slope and how to add it --> Pol1 fiiting implemented

//#11 How do you know if your phi candidate is a valid candidate or not 
//#12 How do you know if your D0 candiate is a D0 or D0Bar in the tree.
//#13 How do you know that the D0 and D0 bar daughter tracks that you have written corresponds to correct candidate or not ?
//#14 Have you checked and matched the Accepted, Rejected and Full Histogram sum for checking.
//#15 Have you checked the trigger mask mapping, that if it was triggered D0 or D0Bar, how is the flag behavingin that case, 
//    are the counts for D0 and D0Bar same, I mean are the tracks matching, have you taken the duplicacy in the account. 
//#16 Massif, mem-check, heap Profiler and CallGrind.check of the code, 
//#17 Add debugger to process and Learn how to use VS code default debugger

//#18 Add Lossy Compression
//#19 Make Versioned Table out of the DrTrack table
//    V 000 - Baisc info (only necessary occ info)
//    V 001 - Full Occ info()
//    V 010 - Lossy Occ info
//

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

template<typename T>
void PrintTime(T Start, std::string String){
  auto Stop = std::chrono::high_resolution_clock::now();
  auto Duration = duration_cast<std::chrono::microseconds>(Stop-Start);
  LOG(info)<<String<<float(Duration.count())/float(1000000)<<" seconds";//<<endl;
}

int32_t nBCsPerOrbit = o2::constants::lhc::LHCMaxBunches;

int hPt_BinFor1        = -1; //= recoEvent.get<TH1>(HIST("Event/track/h_Pt"  )->FindBin(1))
int hQPt_BinFor_Plus1  = -1; //= recoEvent.get<TH1>(HIST("Event/track/h_QPt" )->FindBin(-1))
int hQPt_BinFor_Minus1 = -1; //= recoEvent.get<TH1>(HIST("Event/track/h_QPt" )->FindBin(-1))

namespace o2::aod
{
namespace resopart
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);
DECLARE_SOA_COLUMN(Px, px, float);
DECLARE_SOA_COLUMN(Py, py, float);
DECLARE_SOA_COLUMN(Pz, pz, float);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(Mass, mass, float);
DECLARE_SOA_COLUMN(DauTrackACharge, dauTrackACharge, short);
DECLARE_SOA_COLUMN(DauTrackBCharge, dauTrackBCharge, short);
DECLARE_SOA_COLUMN(DauTrackAType, dauTrackAType, int8_t);
DECLARE_SOA_COLUMN(DauTrackBType, dauTrackBType, int8_t);
DECLARE_SOA_COLUMN(DauTrackAId, dauTrackAId, int64_t);
DECLARE_SOA_COLUMN(DauTrackBId, dauTrackBId, int64_t);
} // namespace lambdatrack
DECLARE_SOA_TABLE(ResoParts, "AOD", "RESOPARTS", o2::soa::Index<>,
                  resopart::CollisionId,
                  resopart::Px,
                  resopart::Py,
                  resopart::Pz,
                  resopart::Pt,
                  resopart::Eta,
                  resopart::Phi,
                  resopart::Mass,
                  resopart::DauTrackACharge,
                  resopart::DauTrackBCharge,
                  resopart::DauTrackAType,
                  resopart::DauTrackBType,
                  resopart::DauTrackAId,
                  resopart::DauTrackBId);
using ResoPart = ResoParts::iterator;
} // namespace o2::aod


struct trackqapidderivedata{

  Service<o2::ccdb::BasicCCDBManager> ccdb;

  // declare production of tables
  Produces<aod::DrCollisions> DerivedCollisions;
  Produces<aod::DrTracks> DerivedTracks;
  Produces<aod::DrFV0As> DerivedFV0As;
  Produces<aod::DrFT0s> DerivedFT0s;
  Produces<aod::DrV0s> DerivedV0s;
  Produces<aod::DrPhis> DerivedPhis;
  Produces<aod::DrD0s> DerivedD0s;
  Produces<aod::OCCS> BuildOCCS;
  Produces<aod::DrTracksMothers> DerivedDrTracksMothers;

  //Histogram registry;
  HistogramRegistry recoTracks  {"recoTracks"  , {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry recoEvent   {"recoEvent"   , {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  //Configurables
  // Event Selection
    // Configurable<float> cutZvertex{"cutZvertex", 10.0f, "Accepted z-vertex range (cm)"};
  
  Configurable<int> customOrbitOffset{"customOrbitOffset", 0, "customOrbitOffset for MC"};
  Configurable<float> factorMB{"factorMB", 0.01, "factorMB"};
  Configurable<float> factorPt{"factorPt", 0.01, "factorPt"};
  Configurable<bool> cfgDoSkipDF{"cfgDoSkipDF", false, "cfgDoSkipDF"};
  Configurable<int > cfgNSkipDF {"cfgNSkipDF" , 3    , "cfgNSkipDF"};

  Configurable<bool> cfgDoDSOfTrack{"cfgDoDSOfTrack", false, "cfgDoDSOfTrack"};
  Configurable<bool> cfgDoDSOfV0s  {"cfgDoDSOfV0s"  , false, "cfgDoDSOfV0s"};
  Configurable<bool> cfgDoDSOfD0s  {"cfgDoDSOfD0s"  , false, "cfgDoDSOfD0s"};
  Configurable<bool> cfgDoDSOfPhi  {"cfgDoDSOfPhi"  , false, "cfgDoDSOfPhi"};

  Configurable<bool> cfgDoTableOfTrack{"cfgDoDSOfTrack", false, "cfgDoDSOfTrack"};
  Configurable<bool> cfgDoTableOfV0s  {"cfgDoDSOfV0s"  , false, "cfgDoDSOfV0s"};
  Configurable<bool> cfgDoTableOfD0s  {"cfgDoDSOfD0s"  , false, "cfgDoDSOfD0s"};
  Configurable<bool> cfgDoTableOfPhi  {"cfgDoDSOfPhi"  , false, "cfgDoDSOfPhi"};

  //check Fill condition for downsampling

  Configurable<bool> doV0andTrackMatchingCheck{"doV0andTrackMatchingCheck", true, "doV0andTrackMatchingCheck"};
  Configurable<bool> doD0andTrackMatchingCheck{"doD0andTrackMatchingCheck", true, "doD0andTrackMatchingCheck"};

  //Lambda 1.115683 GeV/c2  [1.11,1.12]   
  //K0Short 0.4987 GevV/c2  [0.485,0.51]
  //Gamma   0 GevV/c2  [-0.01,0.02]

  Configurable<float> cfgMLowLambda          {"cfgMLowLambda"          , 1.10, "cfgMLowLambda"          };
  Configurable<float> cfgMLowAntiLambda      {"cfgMLowAntiLambda"      , 1.10, "cfgMLowAntiLambda"      };
  Configurable<float> cfgMLowK0Short         {"cfgMLowK0Short"         , 0.47, "cfgMLowK0Short"         };
  Configurable<float> cfgMLowGamma           {"cfgMLowGamma"           ,-0.01, "cfgMLowGamma"           };
  Configurable<float> cfgMLowHypertriton     {"cfgMLowHypertriton"     ,  2.9, "cfgMLowHypertriton"     };
  Configurable<float> cfgMLowAntiHypertriton {"cfgMLowAntiHypertriton" ,  2.9, "cfgMLowAntiHypertriton" };

  Configurable<float> cfgMHighLambda         {"cfgMHighLambda"         , 1.13, "cfgMHighLambda"         };
  Configurable<float> cfgMHighAntiLambda     {"cfgMHighAntiLambda"     , 1.13, "cfgMHighAntiLambda"     };
  Configurable<float> cfgMHighK0Short        {"cfgMHighK0Short"        , 0.52, "cfgMHighK0Short"        };
  Configurable<float> cfgMHighGamma          {"cfgMHighGamma"          , 0.02, "cfgMHighGamma"          };
  Configurable<float> cfgMHighHypertriton    {"cfgMHighHypertriton"    ,  3.1, "cfgMHighHypertriton"    };
  Configurable<float> cfgMHighAntiHypertriton{"cfgMHighAntiHypertriton",  3.1, "cfgMHighAntiHypertriton"};

  Configurable<float> cfgMLowD0 {"cfgMLowD0", 1.862, "cfgMLowD0"};
  Configurable<float> cfgMLowD0Bar {"cfgMLowD0Bar", 1.862, "cfgMLowD0Bar"};
  Configurable<float> cfgMHighD0   {"cfgMHighD0", 1.868, "cfgMHighD0"};
  Configurable<float> cfgMHighD0Bar {"cfgMHighD0Bar", 1.868, "cfgMHighD0Bar"};

  Configurable<float> cfgMLowPhi {"cfgMLowPhi", 1.013, "cfgMLowPhi"};
  Configurable<float> cfgMHighPhi {"cfgMHighPhi", 1.026, "cfgMHighPhi"};
  
  void init(InitContext const&){
    LOGF(info,"Starting init");

    // auto grplhcif = ccdb->getForTimeStamp<o2::parameters::GRPLHCIFData>("GLO/Config/GRPLHCIF", tsSOR);
    // bcPatternB = grplhcif->getBunchFilling().getBCPattern();
    // auto vector = ccdb->getForTimeStamp<std::vector<Long64_t>>(1731764207045, 1731850607045);
/*
    auto& mgr = o2::ccdb::BasicCCDBManager::instance();
    // // mgr.setURL("http://ourccdbserverver.cern.ch");
    // // http://ccdb-test.cern.ch:8080/browse/Users/r/raverma?report=true
    mgr.setURL("http://ccdb-test.cern.ch:8080");
    mgr.setCaching(true);
    auto ccdb_obj = mgr.getForTimeStamp<TList>("Users/r/raverma/",1);
    if(!ccdb_obj){
      LOG(info)<<"DEBUG :: CCDB OBJECT NOT FOUND";
    } else{
      LOG(info)<<"DEBUG :: CCDB OBJECT FOUND";
    }

    ccdb_obj->Print();

    // auto myVect = reinterpret_cast<TH3F*>(ccdb_obj->FindObject(Form("%s", str.c_str())));
        auto myVect = ccdb_obj->FindObject("TVectorT<float>");
        if(!myVect){
          LOG(info)<<"DEBUG :: Vect OBJECT NOT FOUND";
        } else{
          LOG(info)<<"DEBUG :: Vect OBJECT FOUND";
          myVect->Print();
        }
        

    // Iterate over the TList and print object names
    TIter next(ccdb_obj);
    TObject* obj = nullptr;

    while ((obj = next())) {
        LOG(info) <<"DEBUG :: Obj Name = "<<obj->GetName();
    }

    //RESET the URL otherwise the process functions which contains ccdb lookups will fail
    mgr.setURL("http://alice-ccdb.cern.ch");
*/
    //Axes
      AxisSpec vertexZAxis        = {30 ,-15. ,15.  , "vrtx_{Z} [cm]"};

      AxisSpec Axis_Fine_eta       = {2000 ,-10.0f ,10.0f ,"#eta"};
      AxisSpec Axis_Fine_phi       = {800  ,-1.0f  ,7.0f  ,"#phi"};

    //Event Selection
      recoEvent.add("hCollisionCount","hCollisionCount",{HistType::kTH1D, {{1,0,1}}});
      recoEvent.add("hVertexXRec"    , "hVertexXRec"   ,{HistType::kTH1D,{{10000, -0.2, 0.2}}});
      recoEvent.add("hVertexYRec"    , "hVertexYRec"   ,{HistType::kTH1D,{{10000, -0.2, 0.2}}});
      recoEvent.add("hVertexZRec"    , "hVertexZRec"   ,{HistType::kTH1F, {vertexZAxis}});
      recoEvent.add("h_nBCsPerTF"    , "h_nBCsPerTF"   ,{HistType::kTH1F, {{100,114040,114060}}});//114048
      recoEvent.add("h_nBCinTF"      , "h_nBCinTF"     ,{HistType::kTH1F, {{1000,0,1000000}}});
    //

      recoEvent.add("Event/track/Full/h_Pt" , "h_Pt" , {HistType::kTH1F, {{250,0.1,50}}});
      recoEvent.add("Event/track/Full/h_QPt", "h_QPt", {HistType::kTH1F, {{100,-5,5}}});

      recoEvent.addClone("Event/track/Full/","Event/track/Rejected/");
      recoEvent.addClone("Event/track/Rejected/","Event/track/Accepted/");

      recoEvent.add("Event/track/isDownsampledMB" , "isDownsampledMB" , kTH1F, {{129*2,-1, 128}});  
      recoEvent.add("Event/track/isDownsampledPT" , "isDownsampledPT" , kTH1F, {{129*2,-1, 128}});  
      recoEvent.add("Event/track/isDownsampledQPT", "isDownsampledQPT", kTH1F, {{129*2,-1, 128}}); 
      recoEvent.add("Event/track/FullMask"        , "FullMask"        , kTH1F, {{129*2,-1, 128}});

      recoEvent.add("Event/track/probRatioPt" , "probRatioPt" , kTH1F, {{400,-20, 20}});  
      recoEvent.add("Event/track/probRatioQPT", "probRatioQPT", kTH1F, {{400,-20, 20}}); 

      recoEvent.add("Event/V0/h_Pt" , "h_Pt" , {HistType::kTH1F, {{250,0.1,50}}});

      recoEvent.addClone("Event/V0/","Event/Lambda/"         );
      recoEvent.addClone("Event/V0/","Event/AntiLambda/"     );
      recoEvent.addClone("Event/V0/","Event/K0Short/"        );
      recoEvent.addClone("Event/V0/","Event/Gamma/"          );
      recoEvent.addClone("Event/V0/","Event/Hypertriton/"    );
      recoEvent.addClone("Event/V0/","Event/AntiHypertriton/");

      recoEvent.addClone("Event/V0/","Event/Phi/"   );
      recoEvent.addClone("Event/V0/","Event/D0/"    );
      recoEvent.addClone("Event/V0/","Event/D0Bar/" );
      
      //Invariant Mass QAPlots
      
      const AxisSpec axisV0Mass     = {200, 0.40f, 0.60f, "#it{M}_{inv}^{K_{s}^{0}} [GeV/#it{c}^{2}]"}; 
      const AxisSpec axisGammaMass  = {300, -0.01f, 0.05f, "#it{M}_{inv}^{#gamma} [GeV/#it{c}^{2}]"};
      const AxisSpec axisHypertritonMass = {500, 2.75f, 3.25f, "#it{M}_{inv}^{Hypertriton} [GeV/#it{c}^{2}]"};

      const AxisSpec axisK0sMass    = {200, 0.40f, 0.60f, "#it{M}_{inv}^{K_{s}^{0}} [GeV/#it{c}^{2}]"};
      const AxisSpec axisLambdaMass = {200, 1.f, 1.2f, "#it{M}_{inv}^{#Lambda} [GeV/#it{c}^{2}]"};

      const AxisSpec axisPhiMass = {200, 0.99f, 1.07f, "#it{M}_{inv}^{#phi} [GeV/#it{c}^{2}]"};

      const AxisSpec axisD0Mass  = {200, 1.860f, 1.870f, "#it{M}_{inv}^{D0} [GeV/#it{c}^{2}]"};

      recoEvent.add("Event/K0Short/h_Mass_Full"    ,"Event/K0Short/h_Mass_Full"    , kTH1F, {axisK0sMass});
      recoEvent.add("Event/K0Short/h_Mass_Rejected","Event/K0Short/h_Mass_Rejected", kTH1F, {axisK0sMass});
      recoEvent.add("Event/K0Short/h_Mass_Accepted","Event/K0Short/h_Mass_Accepted", kTH1F, {axisK0sMass});

      recoEvent.add("Event/Lambda/h_Mass_Full"    ,"Event/Lambda/h_Mass_Full"    , kTH1F, {axisLambdaMass});
      recoEvent.add("Event/Lambda/h_Mass_Rejected","Event/Lambda/h_Mass_Rejected", kTH1F, {axisLambdaMass});
      recoEvent.add("Event/Lambda/h_Mass_Accepted","Event/Lambda/h_Mass_Accepted", kTH1F, {axisLambdaMass});
      
      recoEvent.add("Event/AntiLambda/h_Mass_Full"    ,"Event/AntiLambda/h_Mass_Full"    , kTH1F, {axisLambdaMass});
      recoEvent.add("Event/AntiLambda/h_Mass_Rejected","Event/AntiLambda/h_Mass_Rejected", kTH1F, {axisLambdaMass});
      recoEvent.add("Event/AntiLambda/h_Mass_Accepted","Event/AntiLambda/h_Mass_Accepted", kTH1F, {axisLambdaMass});

      recoEvent.add("Event/Gamma/h_Mass_Full"    ,"Event/Gamma/h_Mass_Full"    , kTH1F, {axisGammaMass});
      recoEvent.add("Event/Gamma/h_Mass_Rejected","Event/Gamma/h_Mass_Rejected", kTH1F, {axisGammaMass});
      recoEvent.add("Event/Gamma/h_Mass_Accepted","Event/Gamma/h_Mass_Accepted", kTH1F, {axisGammaMass});

      recoEvent.add("Event/Hypertriton/h_Mass_Full"    ,"Event/Hypertriton/h_Mass_Full"    , kTH1F, {axisHypertritonMass});
      recoEvent.add("Event/Hypertriton/h_Mass_Rejected","Event/Hypertriton/h_Mass_Rejected", kTH1F, {axisHypertritonMass});
      recoEvent.add("Event/Hypertriton/h_Mass_Accepted","Event/Hypertriton/h_Mass_Accepted", kTH1F, {axisHypertritonMass});

      recoEvent.add("Event/AntiHypertriton/h_Mass_Full"    ,"Event/AntiHypertriton/h_Mass_Full"    , kTH1F, {axisHypertritonMass});
      recoEvent.add("Event/AntiHypertriton/h_Mass_Rejected","Event/AntiHypertriton/h_Mass_Rejected", kTH1F, {axisHypertritonMass});
      recoEvent.add("Event/AntiHypertriton/h_Mass_Accepted","Event/AntiHypertriton/h_Mass_Accepted", kTH1F, {axisHypertritonMass});

      recoEvent.add("Event/Phi/h_Mass_Full"    ,"Event/Phi/h_Mass_Full"    , kTH1F, {axisPhiMass});
      recoEvent.add("Event/Phi/h_Mass_Rejected","Event/Phi/h_Mass_Rejected", kTH1F, {axisPhiMass});
      recoEvent.add("Event/Phi/h_Mass_Accepted","Event/Phi/h_Mass_Accepted", kTH1F, {axisPhiMass});

      recoEvent.add("Event/D0/h_Mass_Full"    ,"Event/D0/h_Mass_Full"    , kTH1F, {axisD0Mass});
      recoEvent.add("Event/D0/h_Mass_Rejected","Event/D0/h_Mass_Rejected", kTH1F, {axisD0Mass});
      recoEvent.add("Event/D0/h_Mass_Accepted","Event/D0/h_Mass_Accepted", kTH1F, {axisD0Mass});

      recoEvent.add("Event/D0Bar/h_Mass_Full"    ,"Event/D0Bar/h_Mass_Full"    , kTH1F, {axisD0Mass});
      recoEvent.add("Event/D0Bar/h_Mass_Rejected","Event/D0Bar/h_Mass_Rejected", kTH1F, {axisD0Mass});
      recoEvent.add("Event/D0Bar/h_Mass_Accepted","Event/D0Bar/h_Mass_Accepted", kTH1F, {axisD0Mass});
      
      recoEvent.add("Event/percentV0Filled", "percentV0Filled", kTH1F, {{1030, -1, 102}});
      recoEvent.add("Event/percentD0Filled", "percentD0Filled", kTH1F, {{1030, -1, 102}});
      recoEvent.add("Event/percentPhiFilled", "percentPhiFilled", kTH1F, {{1030, -1, 102}});

      hPt_BinFor1        = recoEvent.get<TH1>(HIST("Event/track/Full/h_Pt" ))->FindBin(1) ;
      hQPt_BinFor_Plus1  = recoEvent.get<TH1>(HIST("Event/track/Full/h_QPt"))->FindBin(1) ;
      hQPt_BinFor_Minus1 = recoEvent.get<TH1>(HIST("Event/track/Full/h_QPt"))->FindBin(-1);

      // Algorithm:
      // * 1D histogram of the pt  for the TPC+ITS tracks  - hisPt(250,0.1,50)
      // * 1D hsitgram of the q/pt for the TPC+ITS tracks  - hisQPt(200,-6,6)
      // * 1D histogram of the pt  for the K0s  - hisPtK0(250,0.1,50)
      // * 1D histogram of the pt  for the Lambda  - hisPtLambda(250,0.1,50)
      // * 1D histogram of the pt  for the Gamma  - hisPtGamma(250,0.1,50)
      // * 1D histogram of the pt  for the Phi  - hisPtPhi(250,0.1,50)
      // * 1D histogram of the pt  for the D0  - hisPtD0(250,0.1,50)

    //
    LOG(info) <<"Printing Event Info ";recoEvent.print();
    LOG(info) << "Finishing Init ";
  }

  enum ParticleTypeEnum{
    kTrack = 0,
    kV0,
    kLambda,
    kAntiLambda,
    kK0Short,
    kGamma,
    kHypertriton,
    kAntiHypertriton,
    kD0,
    kD0Bar,
    kPhi1020
  };

  template<typename T, typename U>
  int BinarySearchTable(T Key, U Table, int low, int high){
    while (low <= high) {
      int mid = low + (high - low) / 2;
      if (Key == Table.iteratorAt(mid).trackId()) {return mid;}

      if (Key >  Table.iteratorAt(mid).trackId()) {low  = mid + 1;} // If Key is greater, ignore left  half, update the low
      else                                        {high = mid - 1;} // If Key is smaller, ignore right half, update the high 
    }
    return -1; // Element is not present
  }

  // o2::aod::ambiguous::TrackId

  template<typename T, typename U>
  int FindInTable(T key, U Table){
    // if ( BinarySearchTable(key, Table, 0, Table.size()-1) != -1) {return true;}
    // return false;
    return BinarySearchTable(key, Table, 0, Table.size()-1);
  }

  int FindTrackIdInList(const int64_t &Key, const std::vector<int64_t> &trackGlobalIndexList, std::vector<int64_t> &DrTrackPositionIndexList, int const& GIListSize){
    int elementPos = -1;
    // if (trackGlobalIndexList.size() > 0) { 
      elementPos = BinarySearchVector(Key,trackGlobalIndexList, 0, GIListSize);
    // } 
    if (elementPos == -1){ return -2000000000;} //element not in track list
    else                 {
      if ( DrTrackPositionIndexList[elementPos] < 0 ) { LOG(info) <<"DEBUG :: ERROR :: ERROR in DrTrackPositionIndexList[elementPos] :: elementPos = "<<elementPos<<" :: DrTrackPositionIndexList[elementPos] = "<<DrTrackPositionIndexList[elementPos];}
      return DrTrackPositionIndexList[elementPos];
    }
  }

  void FillNewListFromOldList( std::vector<int64_t> &NewList, std::vector<int64_t> OldList){
    for(long unsigned int ii = 0; ii < OldList.size(); ii++){
      bool RepeatEntry = false;
      for(long unsigned int jj = 0; jj < NewList.size() ; jj++)
      {if( OldList[ii] == NewList[jj]) { RepeatEntry = true;}}
      if(!RepeatEntry){NewList.push_back(OldList[ii]);}
    }
  }

  void InsertionSortVector(std::vector<int64_t> &UnsortedVector){
    for (long unsigned int i = 1; i < UnsortedVector.size(); i++){
      int currentElement = UnsortedVector[i]; //Element to be Inserted at correct position
      int j; //(j+1) is the correct position of current element
      for (j = i-1; j >= 0 && (UnsortedVector[j] > currentElement) ; j--)
      {UnsortedVector[j+1] = UnsortedVector[j];}
      UnsortedVector[j+1] = currentElement;
    }
  }


  template<typename T>
  bool vectorAscendingSortCheck(const T& vec){
    for(long unsigned int i = 1 ; i < vec.size(); i++){
      if( vec[i] < vec[i-1]) 
      { LOG(info)<<"DEBUG :: Vector unsorted at Position = "<<i;
        return false;
      }
    }
    return true;
  }


  // template <typename int64_t>
  int BinarySearchVector(const int64_t &Key, const std::vector<int64_t> &List, int low, int high)
  {
    while (low <= high) {
      int mid = low + (high - low) / 2;
      if (Key == List[mid]) {
        return mid;
      }

      if (Key > List[mid]) {
        low = mid + 1;
      } // If Key is greater, ignore left  half, update the low
      else {
        high = mid - 1;
      } // If Key is smaller, ignore right half, update the high
    }
    return -1; // Element is not present
  }

  void GetRunInfo(const int& run, int& nBCsPerTF, int64_t& bcSOR)
  {
    auto runDuration = ccdb->getRunDuration(run, true);
    int64_t tsSOR = runDuration.first;
    auto ctpx = ccdb->getForTimeStamp<std::vector<Long64_t>>("CTP/Calib/OrbitReset", tsSOR);
    int64_t tsOrbitReset = (*ctpx)[0];
    uint32_t nOrbitsPerTF = run < 534133 ? 128 : 32;
    int64_t orbitSOR = (tsSOR * 1000 - tsOrbitReset) / o2::constants::lhc::LHCOrbitMUS;
    orbitSOR = orbitSOR / nOrbitsPerTF * nOrbitsPerTF;
    bcSOR = orbitSOR * nBCsPerOrbit + customOrbitOffset * nBCsPerOrbit; // customOrbitOffset is a configurable
    nBCsPerTF = nOrbitsPerTF * nBCsPerOrbit;
  }

  template <typename T>
  void GetTimingInfo(const T& bc, int& lastRun, int32_t& nBCsPerTF, int64_t& bcSOR, uint64_t& time, int64_t& TFidThis, int& bcInTF)
  {
    int run = bc.runNumber();
    if (run != lastRun) { // update run info
      lastRun = run;
      GetRunInfo(run, nBCsPerTF, bcSOR); // update nBCsPerTF && bcSOR
    }
    // update the information
    time = bc.timestamp();
    TFidThis = (bc.globalBC() - bcSOR) / nBCsPerTF;
    bcInTF = (bc.globalBC() - bcSOR) % nBCsPerTF;
  }

  template<typename T>
  bool checkTriggerPt(float pt, float factorDS, float &weight, float &probRatio, T hist){
    float prob  = static_cast<float>(static_cast<double>(hist->GetBinContent(hist->FindBin(pt)))/static_cast<double>(hist->GetEntries()));
    float prob1 = static_cast<float>(static_cast<double>(hist->GetBinContent(hPt_BinFor1))      /static_cast<double>(hist->GetEntries()));
    weight = prob;
    probRatio = prob/prob1;
    bool isTrigger = (gRandom->Rndm()*(prob/prob1))<factorDS;
    return isTrigger;
  }

  template<typename T>
  bool checkTriggerQPt(float Qpt, float factorDS, float &weight, float &probRatio, T hist){
    float prob  = static_cast<float>(static_cast<double>(hist->GetBinContent(hist->FindBin(Qpt)))/static_cast<double>(hist->GetEntries()));
    float prob1 = static_cast<float>(static_cast<double>(hist->GetBinContent(hQPt_BinFor_Plus1)) /static_cast<double>(hist->GetEntries()));
    // float prob1 = static_cast<float>(static_cast<double>(hist->GetBinContent(hQPt_BinFor_Plus1))+static_cast<double>(hist->GetBinContent(hQPt_BinFor_Minus1))/static_cast<double>(hist->GetEntries()));
    weight = prob;
    probRatio = prob/prob1;
    bool isTrigger = (gRandom->Rndm()*(prob/prob1))<factorDS;
    return isTrigger;
  }

  template<typename H>
  float histPDFPol1Evaluation(const H& hist, const float& pointValue, const int& binRangeForFit){
    //Linear fitting and getting the pdf value from the linear fitting  
    int nBins = hist->GetNbinsX();
    int iBinPoint = hist->FindBin(pointValue); //bin of the point
    int iBinLow = iBinPoint - binRangeForFit;
    int iBinUp  = iBinPoint + binRangeForFit;

    //Handling the boundary, overflow and underflow cases 
    if (iBinUp > nBins) { iBinUp = nBins; }  // cout<<"Upper Boundary Point Occurred"<<endl;
    if (iBinLow < 1) { iBinLow = 1; }    // cout<<"Lower Boundary Point Occurred"<<endl;

    // float xFitCenter = hist->GetBinCenter(iBinPoint);
    float xFitLow = hist->GetBinLowEdge(iBinLow);
    float xFitHigh = hist->GetBinLowEdge(iBinUp) + hist->GetBinWidth(iBinUp);

    // cout<<"iBinPoint = "<<iBinPoint<<" ::  bin Range = ["<<iBinLow<<","<<iBinUp<<"]"<<endl;
    // if( xFitCenter < hist->GetBinLowEdge(1)){
    //   cout<<"xFitCenter = "<<xFitCenter<<" :: UnderFlow Case :: fitRange = ["<<xFitLow<<", "<<xFitHigh<<"]"<<endl;
    // } else if( xFitCenter >  hist->GetBinLowEdge(nBins+1)){
    //   cout<<"xFitCenter = "<<xFitCenter<<" :: Overflow Case :: fitRange = ["<<xFitLow<<", "<<xFitHigh<<"]"<<endl;
    // } else {
    //   cout<<"xFitCenter = "<<xFitCenter<<" :: fitRange = ["<<xFitLow<<", "<<xFitHigh<<"]"<<endl;
    // }

    // Check the case if there are no entries in the bin range for fitting
    int firstNonZeroBin = -1;
    int nNonZeroBins = 0; 
    for(int iBin = iBinLow; iBin <= iBinUp ; iBin++ ){
      if  (hist->GetBinContent(iBin) != 0) { nNonZeroBins++;
        if( nNonZeroBins == 1 ) { firstNonZeroBin = iBin;}
        // if( nNonZeroBins >= 2 ) { break;}
      }
    }

    // cout<<"firstNonZeroBin = "<<firstNonZeroBin<<endl;    // cout<<"nNonZeroBins    = "<<nNonZeroBins<<endl;

    float evaluatedFitValue = -999; 
    if ( nNonZeroBins >= 2){
      //Define the fit function
      // TF1 *fitFunction = nullptr;  
      TF1 *fitFunction = new TF1("fitFunction", "pol1", xFitLow, xFitHigh);
      hist->Fit(fitFunction, "RNQ"); //"R" = to fit in the range , "N" = No Drawing, "Q" = Don't print the fitting information
      evaluatedFitValue = fitFunction->Eval(pointValue);
      delete fitFunction;
    } else if ( nNonZeroBins == 1) { // cout<<"Trigger :: nNonZeroBins == 1"<<endl;
      if( firstNonZeroBin == iBinPoint) { //the bin under study is itself non empty // choose that bin content as its value
        evaluatedFitValue = hist->GetBinContent(iBinPoint); // cout<<"Trigger :: firstNonZeroBin == iBinPoint"<<endl;
      }
    } else {
        evaluatedFitValue = 0;
    }
    return evaluatedFitValue;
  }


  template<int motherMode, typename T>
  void getMotherListFromDaughterList( const int64_t &Key, const std::vector<int64_t> &dauList, std::vector<int64_t> &trigV0List, std::vector<int64_t> &newMotherList, const T& V0s){
    for(uint32_t i = 0; i < dauList.size(); i++){
      if(dauList[i] == Key) {
        // LOG(info)<<"DEBUG :: Key :: "<<Key<<" :: pos = "<<i<<" :: trigV0 = "<<trigV0List[i]<<" :: GI of V0 at it = "<<V0s.iteratorAt(trigV0List[i]).globalIndex()<<" :: posTrackId() = "<<V0s.iteratorAt(trigV0List[i]).posTrackId() ;
        if constexpr (motherMode == kV0)
        {
          if(doV0andTrackMatchingCheck){ if(V0s.iteratorAt(trigV0List[i]).posTrackId() != Key && V0s.iteratorAt(trigV0List[i]).negTrackId() != Key){
            LOG(info)<<"DEBUG :: v0 Information checking :: something is wrong :: check the errors :: Indices not mathcing";
            LOG(info)<<"Key :: "<<Key<<" :: ";
          }}
          newMotherList.push_back(trigV0List[i]);
        }
        else if constexpr (motherMode == kD0)
        {
          if(doD0andTrackMatchingCheck){ if(V0s.iteratorAt(trigV0List[i]).prong0Id() != Key && V0s.iteratorAt(trigV0List[i]).prong1Id() != Key){
            LOG(info)<<"DEBUG :: D0 Information checking :: something is wrong :: check the errors :: Indices not mathcing";
            LOG(info)<<"Key :: "<<Key<<" :: ";
          }}
          newMotherList.push_back(trigV0List[i]);
        }
      }
    }
  }

  template<typename T>
  bool selKaon(const T& track){
    if(track.hasTOF()){ //if track has good tof
      if( (std::pow(track.tpcNSigmaKa(),2) + std::pow(track.tofNSigmaKa(),2)) < 9.0){
        return true;
      } else {
        return false;
      }
    } else { //it has bad tof
      if(std::abs(track.tpcNSigmaKa()) < 3.0) {
        return true;
      }else {
        return false;
      }
    }
  }


  using myCollisions = soa::Join<aod::Collisions, aod::Mults, aod::MultsExtra, aod::FT0sCorrected, aod::EvSels
                                ,aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::CentFDDMs, aod::CentNTPVs>;

  using myTracks = soa::Join<aod::Tracks, o2::aod::TracksCov, aod::TracksExtra, aod::TracksDCA, aod::TracksDCACov //, aod::TracksQA>;//, , aod::TrackSelection
                    ,aod::TOFSignal,  aod::pidTOFbeta, aod::pidTOFmass, aod::EvTimeTOFOnly, aod::pidTOFFlags, aod::pidEvTimeFlags
                    ,aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTPCFullEl, aod::pidTPCFullDe
                    ,aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFFullEl, aod::pidTOFFullDe>;

  // using myOccsTable = soa::Join<aod::OccsBCsList, aod::OccsPrim    , aod::OccsT0V0    , aod::ORaod::OccsMultExtra, aod::OccsRobustT0V0Prim, aod::OccsRobustFDDT0V0Prim, aod::OccsRobustNtrackDet, aod::OccsRobustMultExtraTable,
  //                                                 aod::OccsMeanPrim, aod::OccsMeanT0V0, aod::ORaod::OccsMnMultExtra, aod::OccsMeanRobustT0V0Prim, aod::OccsMeanRobustFDDT0V0Prim, aod::OccsMeanRobustNtrackDet, aod::OccsMeanRobustMultExtraTable>;

  using myOccsTable = soa::Join<aod::OccsBCsList, aod::OccsPrim, aod::OccsMeanPrim, 
                    aod::OccsT0V0, aod::OccsMeanT0V0, aod::ORT0V0Prim, aod::OMRT0V0Prim,
                    aod::OccsFDD, aod::OccsMeanFDD, aod::ORFDDT0V0Prim, aod::OMRFDDT0V0Prim,
                    aod::OccsNTrackDet, aod::OccsMeanNTrkDet, aod::ORNtrackDet, aod::OMRNtrackDet,
                    aod::OccsMultExtra, aod::OccsMnMultExtra, aod::ORMultExtra, aod::OMRMultExtra>;

  // using myTrackMeanOccs = soa::Join<aod::TrackMeanOccs0, aod::TrackMeanOccs1, aod::TrackMeanOccs2, aod::TrackMeanOccs3, aod::TrackMeanOccs4,
  //                                                        aod::TrackMeanOccs5, aod::TrackMeanOccs6, aod::TrackMeanOccs7, aod::TrackMeanOccs8>;

  using myTrackMeanOccs = soa::Join<aod::TmoTrackId, aod::TmoPrim , aod::TmoT0V0 , aod::TmoFDD , aod::TmoNTrackDet , aod::TmoMultExtra , aod::TmoRT0V0Prim , aod::TmoRFDDT0V0Prim, aod::TmoRNtrackDet , aod::TmoRMultExtra,
                                                    aod::TwmoPrim, aod::TwmoT0V0, aod::TwmoFDD, aod::TwmoNTrackDet, aod::TwmoMultExtra, aod::TwmoRT0V0Prim, aod::TwmoRFDDT0V0Pri, aod::TwmoRNtrackDet, aod::TwmoRMultExtra>;
                                                         
  using myBCTable = soa::Join<aod::BCsWithTimestamps,aod::OccIndexTable>;
  // using MyTracksQA = aod::TracksQA_000;// aod::TracksQAVersion; //aod::TracksQA
  // using MyTracksQA = aod::TracksQA_002;  //aod::TracksQA_000;// aod::TracksQAVersion; //aod::TracksQA
  using MyTracksQA = aod::TracksQA_003;

  // using tracksQA = aod::TracksQA;
  // For manual sliceBy
  Preslice<myTracks>  TracksPerCollisionPreslice = o2::aod::track::collisionId;
  Preslice<aod::TracksQA>  trackQA_Preslice = o2::aod::trackqa::trackId;
  // Preslice<myTracks>  Tracks_PreSlice = o2::aod::track::globalIndex;
  // Preslice<myTracks>  Tracks_PreSlice = o2::soa::globalIndex;

  //Use Partition after definition of filtered Tracks
  SliceCache cache;
  Partition<myTracks> posTracks = aod::track::signed1Pt > 0.f;// track.sign() is dynamic column so use signed1Pt
  Partition<myTracks> negTracks = aod::track::signed1Pt < 0.f;


  // int lastRun = -1;
  // int64_t bcSOR = -1;     // global bc of the start of the first orbit
  // int32_t nBCsPerTF = -1; // duration of TF in bcs, should be 128*3564 or 32*3564
  // int32_t offsetITSROF = 64;
  // int32_t nBCsPerITSROF = 198;


  //Process the Data
  int iEvent = -1;
  int64_t dfCount = 0;
  int debugCounter = 0;
  int CountCheck = 0;
  int FV0A_Debug = 0;

  uint64_t gCollisionChecker = 0;
  uint64_t gFV0AChecker = 0;
  uint64_t gFT0Checker  = 0;
  uint64_t gTrackChecker = 0;
  uint64_t gV0Checker    = 0;
  uint64_t gPhiChecker = 0;
  uint64_t gD0Checker = 0;
  
  // int      dfCount = 0;
  int lastRun = -999;
  int64_t bcSOR = -1;
  int32_t nBCsPerTF = -1;
  uint64_t time = -1;
  int64_t TFidThis = -1;
  int bcInTF = -1;
  
  bool startDownsampling = false;

  std::chrono::high_resolution_clock::time_point Start0 = std::chrono::high_resolution_clock::now();

  void process( myBCTable const& BCs
               ,myCollisions const& collisions
               ,myTracks const& tracks
               ,MyTracksQA const& tracksQA
               ,o2::aod::Origins const& Origins
               ,o2::aod::AmbiguousTracks const& ambgTracks
               ,o2::aod::FV0As const& FV0As
              //  //,o2::aod::FV0Cs const& FV0Cs
               ,o2::aod::FT0s  const& FT0s
               ,o2::aod::V0Datas  const& V0s
               ,soa::Join<aod::HfCand2Prong, aod::HfSelD0> const& D0s
               ,o2::aod::OccIndexTable const& occIdxTable
               ,myOccsTable const& occTables
               ,myTrackMeanOccs const& trackMeanOccs
              //  ,aod::ResoParts const& resoParts
              )
  {
    lastRun = -999;
    nBCsPerTF = -999;
    bcSOR = -999;
    dfCount++;

    auto Start1 = std::chrono::high_resolution_clock::now();    
    PrintTime(Start1, Form("DEBUG :: df_%ld :: DF Reading :: DF Read Time :: ",dfCount));
    PrintTime(Start0, Form("DEBUG :: df_%ld :: DF Reading :: Elapsed Time :: ",dfCount));

    // if(dfCount > 10) {return;}
    // LOG(info)<<"DEBUG :: df_"<<dfCount;
    LOG(info)<<"DEBUG :: df_"<<dfCount<<" :: DF_"<<Origins.iteratorAt(0).dataframeID()<<" :: collisions.size() = "<<collisions.size()<<" :: tracks.size() = "<<tracks.size()<<" :: tracksQA.size() = "<<tracksQA.size();

    LOG(info)<<"DEBUG :: df_"<<dfCount<<" :: DF_"<<Origins.iteratorAt(0).dataframeID()
                                                                                      <<" :: BCs.size() = "<<BCs.size()
                                                                                      <<" :: occIdxTable.size() = "<<occIdxTable.size()
                                                                                      <<" :: occTables.size() = "<<occTables.size()
                                                                                      <<" :: tracKMeanOccs.size() = "<<trackMeanOccs.size()
                                                                                      <<" :: V0s.size() = "<<V0s.size()
                                                                                      <<" :: D0s.size() = "<<D0s.size()
                                                                                      // <<" :: resoParts.size() = "<<resoParts.size()
                                                                                      // <<" :: FV0As.size() = "<<FV0As.size()
                                                                                      // // <<" :: FV0Cs.size() = "<<FV0Cs.size()
                                                                                      // <<" :: FT0s.size() = "<<FT0s.size()
                                                                                      ;

    // if( dfCount < 10 ) { 
    //   LOG(info)<<"DEBUG :: df_"<<dfCount<<" :: DF_"<<Origins.iteratorAt(0).dataframeID()
    //                                     <<" :: FT0s :: posZ = "        <<FT0s.iteratorAt(0).posZ()
    //                                             <<" :: CollTime = "    <<FT0s.iteratorAt(0).collTime()
    //                                             <<" :: IsValidTimeA = "<<FT0s.iteratorAt(0).isValidTimeA()
    //                                             <<" :: IsValidTimeC = "<<FT0s.iteratorAt(0).isValidTimeC()
    //                                             <<" :: IsValidTime = " <<FT0s.iteratorAt(0).isValidTime()
    //                                             <<" :: SumAmpA = "     <<FT0s.iteratorAt(0).sumAmpA()
    //                                             <<" :: SumAmpC = "     <<FT0s.iteratorAt(0).sumAmpC()
    //                                             ;
    // }
    auto StartOcc = std::chrono::high_resolution_clock::now();
    for(auto const& myocc : occTables){
      std::vector<int64_t> myVector = std::vector<int64_t>(myocc.bcsInTFList().begin(),myocc.bcsInTFList().end());
      BuildOCCS(
                  myocc.tfId(),
                  dfCount,
                  std::vector<int64_t>(myocc.bcsInTFList().begin(),myocc.bcsInTFList().end()),
                  std::vector<float>(myocc.occPrimUnfm80().begin(), myocc.occPrimUnfm80().end()),
                  std::vector<float>(myocc.occFV0AUnfm80().begin(), myocc.occFV0AUnfm80().end()),
                  std::vector<float>(myocc.occFV0CUnfm80().begin(), myocc.occFV0CUnfm80().end()),
                  std::vector<float>(myocc.occFT0AUnfm80().begin(), myocc.occFT0AUnfm80().end()),
                  std::vector<float>(myocc.occFT0CUnfm80().begin(), myocc.occFT0CUnfm80().end()),
                  std::vector<float>(myocc.occFDDAUnfm80().begin(), myocc.occFDDAUnfm80().end()),
                  std::vector<float>(myocc.occFDDCUnfm80().begin(), myocc.occFDDCUnfm80().end()),

                  std::vector<float>(myocc.occNTrackITSUnfm80().begin(), myocc.occNTrackITSUnfm80().end()),
                  std::vector<float>(myocc.occNTrackTPCUnfm80().begin(), myocc.occNTrackTPCUnfm80().end()),
                  std::vector<float>(myocc.occNTrackTRDUnfm80().begin(), myocc.occNTrackTRDUnfm80().end()),
                  std::vector<float>(myocc.occNTrackTOFUnfm80().begin(), myocc.occNTrackTOFUnfm80().end()),
                  std::vector<float>(myocc.occNTrackSizeUnfm80().begin(), myocc.occNTrackSizeUnfm80().end()),
                  std::vector<float>(myocc.occNTrackTPCAUnfm80().begin(), myocc.occNTrackTPCAUnfm80().end()),
                  std::vector<float>(myocc.occNTrackTPCCUnfm80().begin(), myocc.occNTrackTPCCUnfm80().end()),
                  std::vector<float>(myocc.occNTrackITSTPCUnfm80().begin(), myocc.occNTrackITSTPCUnfm80().end()),
                  std::vector<float>(myocc.occNTrackITSTPCAUnfm80().begin(), myocc.occNTrackITSTPCAUnfm80().end()),
                  std::vector<float>(myocc.occNTrackITSTPCCUnfm80().begin(), myocc.occNTrackITSTPCCUnfm80().end()),

                  std::vector<float>(myocc.occMultNTracksHasITSUnfm80().begin(), myocc.occMultNTracksHasITSUnfm80().end()), 
                  std::vector<float>(myocc.occMultNTracksHasTPCUnfm80().begin(), myocc.occMultNTracksHasTPCUnfm80().end()), 
                  std::vector<float>(myocc.occMultNTracksHasTOFUnfm80().begin(), myocc.occMultNTracksHasTOFUnfm80().end()), 
                  std::vector<float>(myocc.occMultNTracksHasTRDUnfm80().begin(), myocc.occMultNTracksHasTRDUnfm80().end()), 
                  std::vector<float>(myocc.occMultNTracksITSOnlyUnfm80().begin(), myocc.occMultNTracksITSOnlyUnfm80().end()), 
                  std::vector<float>(myocc.occMultNTracksTPCOnlyUnfm80().begin(), myocc.occMultNTracksTPCOnlyUnfm80().end()), 
                  std::vector<float>(myocc.occMultNTracksITSTPCUnfm80().begin(), myocc.occMultNTracksITSTPCUnfm80().end()), 
                  std::vector<float>(myocc.occMultAllTracksTPCOnlyUnfm80().begin(), myocc.occMultAllTracksTPCOnlyUnfm80().end()), 

                  std::vector<float>(myocc.occRobustT0V0PrimUnfm80().begin(), myocc.occRobustT0V0PrimUnfm80().end()),
                  std::vector<float>(myocc.occRobustFDDT0V0PrimUnfm80().begin(), myocc.occRobustFDDT0V0PrimUnfm80().end()),
                  std::vector<float>(myocc.occRobustNtrackDetUnfm80().begin(), myocc.occRobustNtrackDetUnfm80().end()),
                  std::vector<float>(myocc.occRobustMultExtraTableUnfm80().begin(), myocc.occRobustMultExtraTableUnfm80().end()),

                  myocc.meanOccPrimUnfm80(),
                  myocc.meanOccFV0AUnfm80(),
                  myocc.meanOccFV0CUnfm80(),
                  myocc.meanOccFT0AUnfm80(),
                  myocc.meanOccFT0CUnfm80(),
                  myocc.meanOccFDDAUnfm80(),
                  myocc.meanOccFDDCUnfm80(),

                  myocc.meanOccNTrackITSUnfm80(),
                  myocc.meanOccNTrackTPCUnfm80(),
                  myocc.meanOccNTrackTRDUnfm80(),
                  myocc.meanOccNTrackTOFUnfm80(),
                  myocc.meanOccNTrackSizeUnfm80(),
                  myocc.meanOccNTrackTPCAUnfm80(),
                  myocc.meanOccNTrackTPCCUnfm80(),
                  myocc.meanOccNTrackITSTPCUnfm80(),
                  myocc.meanOccNTrackITSTPCAUnfm80(),
                  myocc.meanOccNTrackITSTPCCUnfm80(),

                  myocc.meanOccMultNTracksHasITSUnfm80(),
                  myocc.meanOccMultNTracksHasTPCUnfm80(),
                  myocc.meanOccMultNTracksHasTOFUnfm80(),
                  myocc.meanOccMultNTracksHasTRDUnfm80(),
                  myocc.meanOccMultNTracksITSOnlyUnfm80(),
                  myocc.meanOccMultNTracksTPCOnlyUnfm80(),
                  myocc.meanOccMultNTracksITSTPCUnfm80(),
                  myocc.meanOccMultAllTracksTPCOnlyUnfm80(),

                  myocc.meanOccRobustT0V0PrimUnfm80(),
                  myocc.meanOccRobustFDDT0V0PrimUnfm80(),
                  myocc.meanOccRobustNtrackDetUnfm80(),
                  myocc.meanOccRobustMultExtraTableUnfm80()
      );
    }
    PrintTime(StartOcc, Form("DEBUG :: df_%ld :: Occ Table Time :: ",dfCount));

    auto StartColl = std::chrono::high_resolution_clock::now();

    int iColl = -1;
    if(collisions.size() == 0){
    return;}

    int lastRun = -1;
    auto bc = BCs.begin();
    int  run = bc.runNumber();
    // auto time = bc.timestamp();

    auto runDuration = ccdb->getRunDuration(run, true);    
    int64_t tsSOR = runDuration.first;
    // auto grplhcif = ccdb->getForTimeStamp<o2::parameters::GRPLHCIFData>("GLO/Config/GRPLHCIF", tsSOR);
    // bcPatternB = grplhcif->getBunchFilling().getBCPattern();
    auto ctpx = ccdb->getForTimeStamp<std::vector<Long64_t>>("CTP/Calib/OrbitReset", tsSOR);
    int64_t tsOrbitReset = (*ctpx)[0];
    uint32_t nOrbitsPerTF = run < 534133 ? 128 : 32;
    int64_t orbitSOR = (tsSOR * 1000 - tsOrbitReset) / o2::constants::lhc::LHCOrbitMUS;
    orbitSOR = orbitSOR / nOrbitsPerTF * nOrbitsPerTF;
    bcSOR = orbitSOR * nBCsPerOrbit + customOrbitOffset * nBCsPerOrbit;
    nBCsPerTF = nOrbitsPerTF * nBCsPerOrbit;
    // nBCsPerITSROF = (run >= 543437 && run <= 545367) ? 594 : 198;

    TFidThis = (bc.globalBC() - bcSOR) / nBCsPerTF;
    bcInTF = (bc.globalBC() - bcSOR) % nBCsPerTF;      

    const int nBCinTF = 114048;  
    // const int nBCinDrift=114048/32;  /// to get from ccdb

    std::vector<int64_t> TFIDList;
    for(const auto& collision : collisions){
      iColl++;
      const auto& bc = collision.bc_as<myBCTable>();
      GetTimingInfo(bc, lastRun, nBCsPerTF, bcSOR, time, TFidThis, bcInTF);

      if (nBCsPerTF > nBCinTF) {
        LOG(error) << "DEBUG :: FATAL ERROR :: nBCsPerTF > nBCinTF i.e " << nBCsPerTF << " > " << nBCinTF << " will cause crash in further process";
        return;
      }
    
      // std::vector<int64_t> &TF_BC_Map = occ_Prim_Unfm_80[TFidThis] ;
      // BC_TF_Map
      
      const uint64_t collIdx = collision.globalIndex();
      const auto TracksTable_perColl =   tracks.sliceBy( TracksPerCollisionPreslice, collIdx);
      
      int nTrack_PVC = 0;
      int nTrack_ITS = 0;
      int nTrack_TPC = 0;
      int nTrack_TRD = 0;
      int nTrack_TOF = 0;
      int nTrackTPC_A = 0;
      int nTrackTPC_C = 0;
      int nTrackITS_TPC = 0;
      int nTrackITS_TPC_A = 0;
      int nTrackITS_TPC_C = 0;
      for(const auto& track : TracksTable_perColl){
        if(track.isPVContributor()) {nTrack_PVC++;} //D 	isPVContributor 	bool 	Run 3: Has this track contributed to the collision vertex fit
        if(track.hasITS()         ) {nTrack_ITS++;} //D 	hasITS 	bool 	Flag to check if track has a ITS match
        if(track.hasTPC()         ) {nTrack_TPC++;
          if( track.eta() <= 0.0  ) {nTrackTPC_A++;} //includes tracks at eta zero as well.
          else                      {nTrackTPC_C++;}
        } //D 	hasTPC 	bool 	Flag to check if track has a TPC match
        if(track.hasTRD()         ) {nTrack_TRD++;} //D 	hasTRD 	bool 	Flag to check if track has a TRD match
        if(track.hasTOF()         ) {nTrack_TOF++;} //D 	hasTOF 	bool 	Flag to check if track has a TOF measurement
        if(track.hasITS() && track.hasTPC()) {nTrackITS_TPC++;
          if( track.eta() <= 0.0  ) {nTrackITS_TPC_A++;} //includes tracks at eta zero as well.
          else                      {nTrackITS_TPC_C++;}
        }
      }// track loop

      DerivedCollisions(
       collision.globalIndex()+gCollisionChecker   //Global Index
      ,collision.globalIndex()                     //Original Index
      ,bc.globalBC()  // globalBC
      ,dfCount
      ,collision.bcId()
      ,collision.posX()
      ,collision.posY()
      ,collision.posZ()
      ,collision.covXX()
      ,collision.covXY()
      ,collision.covXZ()
      ,collision.covYY()
      ,collision.covYZ()
      ,collision.covZZ()
      ,collision.flags()
      ,collision.chi2()
      ,collision.numContrib()
      ,collision.collisionTime()
      ,collision.collisionTimeRes()
      ,time
      ,TFidThis
      ,bcInTF      
      ,collision.multFV0A()
      ,collision.multFV0C()
      // ,collision.multFV0M()
      ,collision.multFT0A()
      ,collision.multFT0C()
      // ,collision.multFT0M()
      ,collision.multFDDA() // 	float 	
      ,collision.multFDDC() // 	float 	
      // collision.multFDDM() // 	float
      ,collision.multZNA () //	float 	
      ,collision.multZNC () //	float 	
      ,collision.multZEM1() // 	float 	
      ,collision.multZEM2() // 	float 	
      ,collision.multZPA () //	float 	
      ,collision.multZPC () //	float 	
      ,nTrack_PVC // int 
      ,nTrack_ITS // int 
      ,nTrack_TPC // int 
      ,nTrack_TRD // int 
      ,nTrack_TOF // int
      ,nTrackTPC_A
      ,nTrackTPC_C
      ,nTrackITS_TPC
      ,nTrackITS_TPC_A
      ,nTrackITS_TPC_C
      ,TracksTable_perColl.size()
      ,collision.multNTracksHasITS()
      ,collision.multNTracksHasTPC()
      ,collision.multNTracksHasTOF()
      ,collision.multNTracksHasTRD()
      ,collision.multNTracksITSOnly()
      ,collision.multNTracksTPCOnly()
      ,collision.multNTracksITSTPC()
      ,collision.multAllTracksTPCOnly()

      ,collision.sel7()  //	bool	Event selection decision based on V0A & V0C
      ,collision.sel8()  //	bool	Event selection decision based on TVX
      ,collision.foundBCId()  //	int	BC entry index in BCs table (-1 if doesn't exist)
      ,collision.foundFT0Id()  //	int	FT0 entry index in FT0s table (-1 if doesn't exist)
      ,collision.foundFV0Id()  //	int	FV0 entry index in FV0As table (-1 if doesn't exist)
      ,collision.foundFDDId()  //	int	FDD entry index in FDDs table (-1 if doesn't exist)
      ,collision.foundZDCId()  //	int	ZDC entry index in ZDCs table (-1 if doesn't exist)
      ,bcInTF
      // ,collision.bcInTF()  //	int	Position of a (found) bunch crossing inside a given timeframe
      ,collision.trackOccupancyInTimeRange()  //	int	Occupancy in specified time interval by a number of tracks from nearby collisions
      ,collision.ft0cOccupancyInTimeRange()  //	float	Occupancy in specified time interval by a sum of FT0C amplitudes from nearby collisions

      ,collision.t0ACorrected()	    //float	Collision time A-side, corrected with primary vertex
      ,collision.t0CCorrected()	    //float	Collision time C-side, corrected with primary vertex
      ,collision.t0AC()	            //float	Collision time (A+C)/2
      ,collision.t0ACorrectedValid()	//bool	Was T0ACorrected computable?
      ,collision.t0CCorrectedValid()	//bool	Was T0CCorrected computable?
      ,collision.t0ACValid()	        //bool	Was T0AC computable?
      ,collision.t0resolution()	    //float	Was T0CCorrected computable?
      ,collision.centFV0A()
      ,collision.centFT0M()
      ,collision.centFT0A()
      ,collision.centFT0C()
      ,collision.centFDDM()
      ,collision.centNTPV()
      //,      int occupancy = col.trackOccupancyInTimeRange();
      );
    }

    //collision Loop is over
    PrintTime(StartColl, Form("DEBUG :: df_%ld :: Col Table Time :: ",dfCount));

    auto StartV0HistFill = std::chrono::high_resolution_clock::now();

    bool isDownsampledMB[20]; //For all V0s, D0s, K0Short, Lambda, AntiLambda, Gamma, Hypertriton, AntiHypertriton
    bool isDownsampledPT[20];

    float weightDSPt[20]; 
    float probRatioPt[20];    

    //First loop to Fill histograms
    for (const auto& v0 : V0s){
      bool isLambda          = false;
      bool isAntiLambda      = false;
      bool isK0Short         = false;
      bool isGamma           = false;
      bool isHypertriton     = false;
      bool isAntiHypertriton = false;

      if (cfgMLowLambda          < v0.mLambda()          && v0.mLambda()          < cfgMHighLambda          ){ isLambda          = true;}
      if (cfgMLowAntiLambda      < v0.mAntiLambda()      && v0.mAntiLambda()      < cfgMHighAntiLambda      ){ isAntiLambda      = true;}
      if (cfgMLowK0Short         < v0.mK0Short()         && v0.mK0Short()         < cfgMHighK0Short         ){ isK0Short         = true;}
      if (cfgMLowGamma           < v0.mGamma()           && v0.mGamma()           < cfgMHighGamma           ){ isGamma           = true;}
      if (cfgMLowHypertriton     < v0.mHypertriton()     && v0.mHypertriton()     < cfgMHighHypertriton     ){ isHypertriton     = true;}
      if (cfgMLowAntiHypertriton < v0.mAntiHypertriton() && v0.mAntiHypertriton() < cfgMHighAntiHypertriton ){ isAntiHypertriton = true;}

      recoEvent.fill(HIST("Event/V0/h_Pt"),v0.pt());
      if( isLambda          ){ recoEvent.fill(HIST("Event/Lambda/h_Pt"          ),v0.pt());}
      if( isAntiLambda      ){ recoEvent.fill(HIST("Event/AntiLambda/h_Pt"      ),v0.pt());}
      if( isK0Short         ){ recoEvent.fill(HIST("Event/K0Short/h_Pt"         ),v0.pt());}
      if( isGamma           ){ recoEvent.fill(HIST("Event/Gamma/h_Pt"           ),v0.pt());}
      if( isHypertriton     ){ recoEvent.fill(HIST("Event/Hypertriton/h_Pt"     ),v0.pt());}
      if( isAntiHypertriton ){ recoEvent.fill(HIST("Event/AntiHypertriton/h_Pt" ),v0.pt());}
    }

    PrintTime(StartV0HistFill, Form("DEBUG :: df_%ld :: V0HistFill Time :: ",dfCount));

    auto StartV0DS = std::chrono::high_resolution_clock::now();
    std::vector<int64_t> posV0DauList;
    std::vector<int64_t> negV0DauList;
    std::vector<int64_t> fullV0DauList;
    std::vector<int64_t> trigV0IndexList;

    std::vector<float> weightDSPtV0             ;
    std::vector<float> weightDSPtLambda         ;
    std::vector<float> weightDSPtAntiLambda     ;
    std::vector<float> weightDSPtK0Short        ;
    std::vector<float> weightDSPtGamma          ;
    std::vector<float> weightDSPtHypertriton    ;
    std::vector<float> weightDSPtAntiHypertriton;

    std::vector<int> vecTriggerMaskV0;

    int nV0Triggered = 0;

    for (const auto& v0 : V0s){
      bool isLambda          = false;
      bool isAntiLambda      = false;
      bool isK0Short         = false;
      bool isGamma           = false;
      bool isHypertriton     = false;
      bool isAntiHypertriton = false;

      if (cfgMLowLambda          < v0.mLambda()          && v0.mLambda()          < cfgMHighLambda          ){ isLambda          = true;}
      if (cfgMLowAntiLambda      < v0.mAntiLambda()      && v0.mAntiLambda()      < cfgMHighAntiLambda      ){ isAntiLambda      = true;}
      if (cfgMLowK0Short         < v0.mK0Short()         && v0.mK0Short()         < cfgMHighK0Short         ){ isK0Short         = true;}
      if (cfgMLowGamma           < v0.mGamma()           && v0.mGamma()           < cfgMHighGamma           ){ isGamma           = true;}
      if (cfgMLowHypertriton     < v0.mHypertriton()     && v0.mHypertriton()     < cfgMHighHypertriton     ){ isHypertriton     = true;}
      if (cfgMLowAntiHypertriton < v0.mAntiHypertriton() && v0.mAntiHypertriton() < cfgMHighAntiHypertriton ){ isAntiHypertriton = true;}

      recoEvent.fill(HIST("Event/K0Short/h_Mass_Full"), v0.mK0Short());
      recoEvent.fill(HIST("Event/Lambda/h_Mass_Full"), v0.mLambda());
      recoEvent.fill(HIST("Event/AntiLambda/h_Mass_Full"), v0.mAntiLambda());
      recoEvent.fill(HIST("Event/Gamma/h_Mass_Full"), v0.mGamma());
      recoEvent.fill(HIST("Event/Hypertriton/h_Mass_Full"), v0.mHypertriton());
      recoEvent.fill(HIST("Event/AntiHypertriton/h_Mass_Full"), v0.mAntiHypertriton());

      if(cfgDoSkipDF && (dfCount < cfgNSkipDF)) {continue;}

      weightDSPt[kV0]              = -999.0 , probRatioPt[kV0]              = -999.0;
      weightDSPt[kLambda]          = -999.0 , probRatioPt[kLambda]          = -999.0;
      weightDSPt[kAntiLambda]      = -999.0 , probRatioPt[kAntiLambda]      = -999.0;
      weightDSPt[kK0Short]         = -999.0 , probRatioPt[kK0Short]         = -999.0;
      weightDSPt[kGamma]           = -999.0 , probRatioPt[kGamma]           = -999.0;
      weightDSPt[kHypertriton]     = -999.0 , probRatioPt[kHypertriton]     = -999.0;
      weightDSPt[kAntiHypertriton] = -999.0 , probRatioPt[kAntiHypertriton] = -999.0;

      isDownsampledMB[kV0]              = false;
      isDownsampledMB[kLambda]          = false;
      isDownsampledMB[kAntiLambda]      = false;
      isDownsampledMB[kK0Short]         = false;
      isDownsampledMB[kGamma]           = false;
      isDownsampledMB[kHypertriton]     = false;
      isDownsampledMB[kAntiHypertriton] = false;

      isDownsampledMB[kV0]              = gRandom->Rndm()<factorMB;
      if( isLambda          ){isDownsampledMB[kLambda]          = gRandom->Rndm()<factorMB;}
      if( isAntiLambda      ){isDownsampledMB[kAntiLambda]      = gRandom->Rndm()<factorMB;}
      if( isK0Short         ){isDownsampledMB[kK0Short]         = gRandom->Rndm()<factorMB;}
      if( isGamma           ){isDownsampledMB[kGamma]           = gRandom->Rndm()<factorMB;}
      if( isHypertriton     ){isDownsampledMB[kHypertriton]     = gRandom->Rndm()<factorMB;}
      if( isAntiHypertriton ){isDownsampledMB[kAntiHypertriton] = gRandom->Rndm()<factorMB;}

      isDownsampledPT[kV0]              = false;
      isDownsampledPT[kLambda]          = false;
      isDownsampledPT[kAntiLambda]      = false;
      isDownsampledPT[kK0Short]         = false;
      isDownsampledPT[kGamma]           = false;
      isDownsampledPT[kHypertriton]     = false;
      isDownsampledPT[kAntiHypertriton] = false;

      isDownsampledPT[kV0]              = checkTriggerPt(v0.pt(), factorPt, weightDSPt[kV0]              , probRatioPt[kV0]              , recoEvent.get<TH1>(HIST("Event/V0/h_Pt"              )));
      if( isLambda          ){isDownsampledPT[kLambda]          = checkTriggerPt(v0.pt(), factorPt, weightDSPt[kLambda]          , probRatioPt[kLambda]          , recoEvent.get<TH1>(HIST("Event/Lambda/h_Pt"          )));}
      if( isAntiLambda      ){isDownsampledPT[kAntiLambda]      = checkTriggerPt(v0.pt(), factorPt, weightDSPt[kAntiLambda]      , probRatioPt[kAntiLambda]      , recoEvent.get<TH1>(HIST("Event/AntiLambda/h_Pt"      )));}
      if( isK0Short         ){isDownsampledPT[kK0Short]         = checkTriggerPt(v0.pt(), factorPt, weightDSPt[kK0Short]         , probRatioPt[kK0Short]         , recoEvent.get<TH1>(HIST("Event/K0Short/h_Pt"         )));}
      // if( isGamma           ){isDownsampledPT[kGamma]           = checkTriggerPt(v0.pt(), factorPt, weightDSPt[kGamma]           , probRatioPt[kGamma]           , recoEvent.get<TH1>(HIST("Event/Gamma/h_Pt"           )));}
      // if( isHypertriton     ){isDownsampledPT[kHypertriton]     = checkTriggerPt(v0.pt(), factorPt, weightDSPt[kHypertriton]     , probRatioPt[kHypertriton]     , recoEvent.get<TH1>(HIST("Event/Hypertriton/h_Pt"     )));}
      // if( isAntiHypertriton ){isDownsampledPT[kAntiHypertriton] = checkTriggerPt(v0.pt(), factorPt, weightDSPt[kAntiHypertriton] , probRatioPt[kAntiHypertriton] , recoEvent.get<TH1>(HIST("Event/AntiHypertriton/h_Pt" )));}
      
      int triggerMaskV0 = (1<<0 )*isDownsampledMB[kV0]+(1<<1 )*isDownsampledMB[kLambda]+(1<<2 )*isDownsampledMB[kAntiLambda]+(1<<3 )*isDownsampledMB[kK0Short]
                        + (1<<10)*isDownsampledPT[kV0]+(1<<11)*isDownsampledPT[kLambda]+(1<<12)*isDownsampledPT[kAntiLambda]+(1<<13)*isDownsampledPT[kK0Short];

      if( cfgDoDSOfV0s && (triggerMaskV0 == 0)) {   
        if( isLambda          ){recoEvent.fill(HIST("Event/K0Short/h_Mass_Rejected")        , v0.mK0Short()        );}
        if( isAntiLambda      ){recoEvent.fill(HIST("Event/Lambda/h_Mass_Rejected")         , v0.mLambda()         );}
        if( isK0Short         ){recoEvent.fill(HIST("Event/AntiLambda/h_Mass_Rejected")     , v0.mAntiLambda()     );}
        if( isGamma           ){recoEvent.fill(HIST("Event/Gamma/h_Mass_Rejected")          , v0.mGamma()          );}
        if( isHypertriton     ){recoEvent.fill(HIST("Event/Hypertriton/h_Mass_Rejected")    , v0.mHypertriton()    );}
        if( isAntiHypertriton ){recoEvent.fill(HIST("Event/AntiHypertriton/h_Mass_Rejected"), v0.mAntiHypertriton());}
        continue;
      }
      nV0Triggered++;
      trigV0IndexList.push_back(v0.globalIndex());
      posV0DauList.push_back(v0.posTrackId());
      negV0DauList.push_back(v0.negTrackId());
      weightDSPtV0.push_back(weightDSPt[kV0]);
      weightDSPtLambda.push_back(weightDSPt[kLambda]);
      weightDSPtAntiLambda.push_back(weightDSPt[kAntiLambda]);
      weightDSPtK0Short.push_back(weightDSPt[kK0Short]);
      weightDSPtGamma.push_back(weightDSPt[kGamma]);
      weightDSPtHypertriton.push_back(weightDSPt[kHypertriton]);
      weightDSPtAntiHypertriton.push_back(weightDSPt[kAntiHypertriton]);

      vecTriggerMaskV0.push_back(triggerMaskV0);

      //Creating mother daughter pointers 
      if( isLambda          ){recoEvent.fill(HIST("Event/K0Short/h_Mass_Accepted")        , v0.mK0Short()        );}
      if( isAntiLambda      ){recoEvent.fill(HIST("Event/Lambda/h_Mass_Accepted")         , v0.mLambda()         );}
      if( isK0Short         ){recoEvent.fill(HIST("Event/AntiLambda/h_Mass_Accepted")     , v0.mAntiLambda()     );}
      if( isGamma           ){recoEvent.fill(HIST("Event/Gamma/h_Mass_Accepted")          , v0.mGamma()          );}
      if( isHypertriton     ){recoEvent.fill(HIST("Event/Hypertriton/h_Mass_Accepted")    , v0.mHypertriton()    );}
      if( isAntiHypertriton ){recoEvent.fill(HIST("Event/AntiHypertriton/h_Mass_Accepted"), v0.mAntiHypertriton());}
    }

    recoEvent.fill(HIST("Event/percentV0Filled"), 100. * static_cast<double>(nV0Triggered)/static_cast<double>(V0s.size()));

    PrintTime(StartV0DS, Form("DEBUG :: df_%ld :: V0      DS Time :: ",dfCount));

    auto StartPhiHist = std::chrono::high_resolution_clock::now();

    float mass1 = 0., mass2 = 0.;
    std::array<float, 3> pvec;
    float p = 0., e = 0., minv = 0., pt = 0., eta = 0., phi = 0.;

    int64_t phiRecoTableIdx = 0;
    std::vector<int64_t> phiOIlist;
    std::vector<int64_t> posKaOIlist;
    std::vector<int64_t> negKaOIlist;

    std::vector<float> phiMasslist;    
    std::vector<float> phiPtlist;
    std::vector<float> phiCollisionIdlist;    

    for(const auto& collision : collisions){
      // Get the Partitions
      auto posTracks_perColl = posTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
      auto negTracks_perColl = negTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

      for(const auto& posTrack : posTracks_perColl){
        for(const auto& negTrack : negTracks_perColl){
          // if(posTrack.collisionId() != negTrack.collisionId()) {continue;} //both must belong to same collision.
          if(posTrack.globalIndex() == negTrack.globalIndex()) {
            LOG(info)<<"DEBUG :: ERROR ERROR ERRO :: A tracks is both positive and negative :: track.globalIndex() = "<<posTrack.globalIndex();
            continue; //
          }
        //   //do track selection

          //do Kaon selection
          if(!selKaon(posTrack)){continue;}
          if(!selKaon(negTrack)){continue;}

          //do Phi Reconstruction
          // get resonance track attributes
          mass1 = MassKaonCharged;
          mass2 = MassKaonCharged;
          pvec[0] = posTrack.px() + negTrack.px();
          pvec[1] = posTrack.py() + negTrack.py();
          pvec[2] = posTrack.pz() + negTrack.pz();
          pt = RecoDecay::pt(pvec[0], pvec[1]);
          eta = RecoDecay::eta(pvec);
          phi = RecoDecay::phi(pvec[0], pvec[1]);
          
          // get invariant mass of pair
          p = RecoDecay::p((posTrack.px() + negTrack.px()), (posTrack.py() + negTrack.py()), (posTrack.pz() + negTrack.pz()));
          e = RecoDecay::e(posTrack.px(), posTrack.py(), posTrack.pz(), mass1) + RecoDecay::e(negTrack.px(), negTrack.py(), negTrack.pz(), mass2);
          minv = std::sqrt(RecoDecay::m2(p, e));

          if (cfgMLowPhi < minv && minv < cfgMHighPhi) { 
            // Add QA plots too
            recoEvent.fill(HIST("Event/Phi/h_Pt"), pt);

            phiOIlist.push_back(phiRecoTableIdx);
            posKaOIlist.push_back(posTrack.globalIndex());
            negKaOIlist.push_back(negTrack.globalIndex());

            phiMasslist.push_back(minv);
            phiPtlist.push_back(pt);
            phiCollisionIdlist.push_back(posTrack.collisionId());
            phiRecoTableIdx++;
          }
        }
      }
    }
    PrintTime(StartPhiHist, Form("DEBUG :: df_%ld :: PhiHistFill Time :: ",dfCount));

    auto StartPhiDS = std::chrono::high_resolution_clock::now();

    int nPhiTriggered = 0;
    std::vector<int64_t> trigPhiIndexList;
    std::vector<int64_t> posKaonPhiList;
    std::vector<int64_t> negKaonPhiList;
    std::vector<float> weightDSPtPhi;
    std::vector<int> vecTriggerMaskPhi;
    
    // auto trackIt = tracks.begin();

    for (const auto& phiCand : phiOIlist) {
      bool isPhi = false;
      const auto& posTrack = tracks.rawIteratorAt(posKaOIlist[phiCand]);
      const auto& negTrack = tracks.rawIteratorAt(negKaOIlist[phiCand]);

      // const auto& posTrack = tracks.iteratorAt(posKaOIlist[phiCand]);
      // const auto& negTrack = tracks.iteratorAt(negKaOIlist[phiCand]);

      if( posTrack.globalIndex() != posKaOIlist[phiCand] ) { LOG(error)<<"DEBUG :: ERROR ERROR ERROR :: Raw pointer to track is misleading :: posTrack.globalIndex() = "<<posTrack.globalIndex()<<" != "<<posKaOIlist[phiCand]<<" = posKaOIlist["<<phiCand<<"]"; }
      if( negTrack.globalIndex() != negKaOIlist[phiCand] ) { LOG(error)<<"DEBUG :: ERROR ERROR ERROR :: Raw pointer to track is misleading :: negTrack.globalIndex() = "<<negTrack.globalIndex()<<" != "<<negKaOIlist[phiCand]<<" = negKaOIlist["<<phiCand<<"]"; }

      if (cfgMLowPhi < phiMasslist[phiCand] && phiMasslist[phiCand] < cfgMHighPhi) { isPhi = true;}

      recoEvent.fill(HIST("Event/Phi/h_Mass_Full"), phiMasslist[phiCand]);

      if(cfgDoSkipDF && (dfCount < cfgNSkipDF)) {continue;}

      weightDSPt[kPhi1020] = -999.0 , probRatioPt[kPhi1020] = -999.0;

      isDownsampledMB[kPhi1020] = false;
      if (isPhi   ) { isDownsampledMB[kPhi1020] = gRandom->Rndm()<factorMB;}

      isDownsampledPT[kPhi1020] = false;
      if (isPhi   ) {isDownsampledPT[kPhi1020] = checkTriggerPt(phiPtlist[phiCand], factorPt, weightDSPt[kPhi1020], probRatioPt[kPhi1020], recoEvent.get<TH1>(HIST("Event/Phi/h_Pt")));}

      int triggerMaskPhi = (1<<0)*isDownsampledMB[kPhi1020]+(1<<1 )*isDownsampledPT[kPhi1020];

      if(cfgDoDSOfPhi && (triggerMaskPhi == 0)) { 
        recoEvent.fill(HIST("Event/Phi/h_Mass_Rejected"), phiMasslist[phiCand]);
        continue;}
      nPhiTriggered++;
      trigPhiIndexList.push_back(phiCand);
      posKaonPhiList.push_back(posTrack.globalIndex());
      negKaonPhiList.push_back(negTrack.globalIndex());
      weightDSPtPhi.push_back(weightDSPt[kPhi1020]);

      vecTriggerMaskPhi.push_back(triggerMaskPhi);
      // Add QA plots too
      recoEvent.fill(HIST("Event/Phi/h_Mass_Accepted"), phiMasslist[phiCand]);
    }
    recoEvent.fill(HIST("Event/percentPhiFilled"), 100. * static_cast<double>(nPhiTriggered)/static_cast<double>(phiOIlist.size()));
    PrintTime(StartPhiDS, Form("DEBUG :: df_%ld :: Phi DS Time :: ",dfCount));

    auto StartD0Hist = std::chrono::high_resolution_clock::now();    
    // ----------------------------------- D0 ----------------------------------------- //
    HfHelper hfHelper;

    for (const auto& d0 : D0s) {
      bool isD0 = false;
      bool isD0Bar = false;

      float massD0 = hfHelper.invMassD0ToPiK(d0);
      float massD0Bar = hfHelper.invMassD0barToKPi(d0);

      if (cfgMLowD0 < massD0 && massD0 < cfgMHighD0) { isD0 = true;}
      if (cfgMLowD0Bar < massD0Bar && massD0Bar < cfgMHighD0Bar) { isD0Bar = true;}
      if (isD0   ) { recoEvent.fill(HIST("Event/D0/h_Pt"), d0.pt());}
      if (isD0Bar) { recoEvent.fill(HIST("Event/D0Bar/h_Pt"), d0.pt());}
    }

    PrintTime(StartD0Hist, Form("DEBUG :: df_%ld :: D0HistFill Time :: ",dfCount));

    auto StartD0DS = std::chrono::high_resolution_clock::now();

    int nD0Triggered = 0;
    std::vector<int64_t> trigD0IndexList;
    std::vector<int64_t> prong0D0List;
    std::vector<int64_t> prong1D0List;
    std::vector<float> weightDSPtD0;
    std::vector<float> weightDSPtD0Bar;
    std::vector<int> vecTriggerMaskD0;

    for (const auto& d0 : D0s) {
      bool isD0 = false;
      bool isD0Bar = false;

      float massD0 = hfHelper.invMassD0ToPiK(d0);
      float massD0Bar = hfHelper.invMassD0barToKPi(d0);

      if (cfgMLowD0 < massD0 && massD0 < cfgMHighD0) { isD0 = true;}
      if (cfgMLowD0Bar < massD0Bar && massD0Bar < cfgMHighD0Bar) { isD0Bar = true;}

      recoEvent.fill(HIST("Event/D0/h_Mass_Full"), massD0);
      recoEvent.fill(HIST("Event/D0Bar/h_Mass_Full"), massD0Bar);

      if(cfgDoSkipDF && (dfCount < cfgNSkipDF)) {continue;}

      weightDSPt[kD0] = -999.0 , probRatioPt[kD0] = -999.0;
      weightDSPt[kD0Bar] = -999.0, probRatioPt[kD0Bar] = -999.0;

      isDownsampledMB[kD0] = false;
      isDownsampledMB[kD0Bar] = false;
      if (isD0   ) { isDownsampledMB[kD0] = gRandom->Rndm()<factorMB;}
      if (isD0Bar) { isDownsampledMB[kD0Bar] = gRandom->Rndm()<factorMB;}

      isDownsampledPT[kD0] = false;
      isDownsampledPT[kD0Bar] = false;
      if (isD0   ) {isDownsampledPT[kD0] = checkTriggerPt(d0.pt(), factorPt, weightDSPt[kD0], probRatioPt[kD0], recoEvent.get<TH1>(HIST("Event/D0/h_Pt")));}
      if (isD0Bar) {isDownsampledPT[kD0Bar] = checkTriggerPt(d0.pt(), factorPt, weightDSPt[kD0Bar], probRatioPt[kD0Bar], recoEvent.get<TH1>(HIST("Event/D0Bar/h_Pt")));}

      int triggerMaskD0 = (1<<0)*isDownsampledMB[kD0]+(1<<1 )*isDownsampledMB[kD0Bar]
                        +(1<<10)*isDownsampledPT[kD0]+(1<<11)*isDownsampledPT[kD0Bar];

      if(cfgDoDSOfD0s && (triggerMaskD0 == 0)) { 
        recoEvent.fill(HIST("Event/D0/h_Mass_Rejected"), massD0);
        recoEvent.fill(HIST("Event/D0Bar/h_Mass_Rejected"), massD0Bar);
        continue;
      }
      nD0Triggered++;
      trigD0IndexList.push_back(d0.globalIndex());
      prong0D0List.push_back(d0.prong0Id());
      prong1D0List.push_back(d0.prong1Id());
      weightDSPtD0.push_back(weightDSPt[kD0]);
      weightDSPtD0Bar.push_back(weightDSPt[kD0Bar]);

      vecTriggerMaskD0.push_back(triggerMaskD0);
      recoEvent.fill(HIST("Event/D0/h_Mass_Accepted"), massD0);
      recoEvent.fill(HIST("Event/D0Bar/h_Mass_Accepted"), massD0Bar);
    }
    recoEvent.fill(HIST("Event/percentD0Filled"), 100. * static_cast<double>(nD0Triggered)/static_cast<double>(D0s.size()));
    // ----------------------------------- D0 ----------------------------------------- //

    PrintTime(StartD0DS, Form("DEBUG :: df_%ld :: D0      DS Time :: ",dfCount));

    auto StartTrackVars = std::chrono::high_resolution_clock::now();

  //Now reading tracks and Calculating the required stuff
    int64_t oldCollisionIndex = -100;
    int CollisionId_negError = 0;
    int CollisionId_posError = 0;
    int CollisionId_NoError  = 0;
    int nAmbgTracks = 0;
    bool hasCollision = false;
    bool isAmbgTrack = false;
    bool LastTrackHadCollision = false;
    bool doCollisionUpdate = false;
    bool doAmbgUpdate = false;
    int trackCounter = -1;

    std::vector<float> Occ_Prim_Unfm_80            ;// = Occs.iteratorAt(bc.occIndex()).occ_Prim_Unfm_80 ();
    std::vector<float> Occ_FV0A_Unfm_80            ;// = Occs.iteratorAt(bc.occIndex()).occ_FV0A_Unfm_80 ();
    std::vector<float> Occ_FV0C_Unfm_80            ;// = Occs.iteratorAt(bc.occIndex()).occ_FV0C_Unfm_80 ();
    std::vector<float> Occ_FT0A_Unfm_80            ;// = Occs.iteratorAt(bc.occIndex()).occ_FT0A_Unfm_80 ();
    std::vector<float> Occ_FT0C_Unfm_80            ;// = Occs.iteratorAt(bc.occIndex()).occ_FT0C_Unfm_80 ();
    std::vector<float> Occ_FDDA_Unfm_80            ;// = Occs.iteratorAt(bc.occIndex()).occ_FDDA_Unfm_80 ();
    std::vector<float> Occ_FDDC_Unfm_80            ;// = Occs.iteratorAt(bc.occIndex()).occ_FDDC_Unfm_80 ();

    std::vector<float> Occ_NTrack_PVC_Unfm_80      ;//= Occs.iteratorAt(bc.occIndex()).occ_NTrack_PVC_Unfm_80   ();
    std::vector<float> Occ_NTrack_ITS_Unfm_80      ;//= Occs.iteratorAt(bc.occIndex()).occ_NTrack_ITS_Unfm_80   ();
    std::vector<float> Occ_NTrack_TPC_Unfm_80      ;//= Occs.iteratorAt(bc.occIndex()).occ_NTrack_TPC_Unfm_80   ();
    std::vector<float> Occ_NTrack_TRD_Unfm_80      ;//= Occs.iteratorAt(bc.occIndex()).occ_NTrack_TRD_Unfm_80   ();
    std::vector<float> Occ_NTrack_TOF_Unfm_80      ;//= Occs.iteratorAt(bc.occIndex()).occ_NTrack_TOF_Unfm_80   ();
    std::vector<float> Occ_NTrackSize_Unfm_80      ;//= Occs.iteratorAt(bc.occIndex()).occ_NTrackSize_Unfm_80   ();
    std::vector<float> Occ_NTrackTPC_A_Unfm_80     ;//= Occs.iteratorAt(bc.occIndex()).occ_NTrackTPC_A_Unfm_80  ();
    std::vector<float> Occ_NTrackTPC_C_Unfm_80     ;//= Occs.iteratorAt(bc.occIndex()).occ_NTrackTPC_C_Unfm_80  ();
    std::vector<float> Occ_NTrackITS_TPC_Unfm_80   ;//= Occs.iteratorAt(bc.occIndex()).occ_NTrackITS_TPC_Unfm_80();
    std::vector<float> Occ_NTrackITS_TPC_A_Unfm_80 ;//= Occs.iteratorAt(bc.occIndex()).occ_NTrackITS_TPC_A_Unfm_80();
    std::vector<float> Occ_NTrackITS_TPC_C_Unfm_80 ;//= Occs.iteratorAt(bc.occIndex()).occ_NTrackITS_TPC_C_Unfm_80();      

    // int64_t oldTFid = -1;

    float   tpcTime0           = -999 ;// 	 	tpc only time0 (mTime0 in TPC track)
    int16_t tpcdcaR            = -999 ;// tpc only DCAr
    int16_t tpcdcaZ            = -999 ;// tpc only DCAz
    uint8_t tpcClusterByteMask = 0    ;// tracklet bitmask - track defining 8 tracklets (152=8*19 rows) bit set if nCluster>thr (default 5)
    uint8_t tpcdEdxMax0R       = 0    ;// TPC dEdxQMax -ROC0/dEdx
    uint8_t tpcdEdxMax1R       = 0    ;// TPC dEdxQMax -ROC1/dEdx
    uint8_t tpcdEdxMax2R       = 0    ;// TPC dEdxQMax -ROC2/dEdx
    uint8_t tpcdEdxMax3R       = 0    ;// TPC dEdxQMax -ROC3/dEdx
    uint8_t tpcdEdxTot0R       = 0    ;// TPC dEdxQtot -ROC0/dEdx
    uint8_t tpcdEdxTot1R       = 0    ;// TPC dEdxQtot -ROC1/dEdx
    uint8_t tpcdEdxTot2R       = 0    ;// TPC dEdxQtot -ROC2/dEdx
    uint8_t tpcdEdxTot3R       = 0    ;// TPC dEdxQtot -ROC3/dEdx
  
    uint8_t deltaRefContParamY = 0;
    uint8_t deltaRefContParamZ = 0;
    uint8_t deltaRefContParamSnp = 0;
    uint8_t deltaRefContParamTgl = 0;
    uint8_t deltaRefContParamQ2Pt = 0;
    uint8_t deltaRefGloParamY = 0;
    uint8_t deltaRefGloParamZ = 0;
    uint8_t deltaRefGloParamSnp = 0;
    uint8_t deltaRefGloParamTgl = 0;
    uint8_t deltaRefGloParamQ2Pt = 0;
    uint8_t deltaTOFdX = 0;
    uint8_t deltaTOFdZ = 0;

    std::vector<int64_t> trackGlobalIndexList;
    std::vector<int64_t> DrTrackPositionIndexList; //Derived track tree position list;
    int trackQAcounter = -1;
    int trackSkippedBecauseOfMask = 0;

    auto tracksWithPid = soa::Attach<myTracks,
                              aod::pidits::ITSNSigmaEl, aod::pidits::ITSNSigmaMu, aod::pidits::ITSNSigmaPi,
                              aod::pidits::ITSNSigmaKa, aod::pidits::ITSNSigmaPr, aod::pidits::ITSNSigmaDe,
                              aod::pidits::ITSNSigmaTr, aod::pidits::ITSNSigmaHe, aod::pidits::ITSNSigmaAl>(tracks);

    PrintTime(StartTrackVars, Form("DEBUG :: df_%ld :: Track Var Time :: ",dfCount));

    auto StartTrackHist = std::chrono::high_resolution_clock::now();

    for(const auto& trackQA : tracksQA){      

      auto const& track = tracksWithPid.iteratorAt(trackQA.trackId()); //Find the track in track table
      if( track.globalIndex() != trackQA.trackId()      ){ LOG(info)<<"DEBUG :: ERROR :: Track and TrackQA Mismatch";} 

      auto const& collision = collisions.iteratorAt(track.collisionId()); //Find the collision in collision table
      if( track.collisionId() != collision.globalIndex()){ LOG(info)<<"DEBUG :: ERROR :: track collId and collID Mismatch";}

      //Checking out of the range errors
      if( trackQA.trackId()   < 0 || tracks.size()     <= trackQA.trackId()  ) { LOG(info)<<"DEBUG :: ERROR :: trackQA has index out of scope :: trackQA.trackId() = "<<trackQA.trackId()<<" :: track.collisionId() = "<<track.collisionId()<<" :: track.signed1Pt() = "<<track.signed1Pt();}
      hasCollision = false;
      isAmbgTrack = false;
      if( track.collisionId() < 0 || collisions.size() <= track.collisionId()) {
        // LOG(info)<<"DEBUG :: ERROR :: track   has index out of scope :: trackQA.trackId() = "<<trackQA.trackId()<<" :: track.collisionId() = "<<track.collisionId()<<" :: track.signed1Pt() = "<<track.signed1Pt();
        if(track.collisionId()< 0){CollisionId_negError++;
          //check if it is an ambiguous track
          int ambgPos = FindInTable(track.globalIndex(), ambgTracks);
          if( ambgPos >= 0 && (ambgTracks.iteratorAt(ambgPos).trackId() == track.globalIndex()) ){
            nAmbgTracks++;
            isAmbgTrack = true;
            // LOG(info)<<"DEBUG :: Track is an ambiguous Track :: ambgId = "<<ambgTracks.iteratorAt(ambgPos).trackId()<<" :: trackId = "<<track.globalIndex();
          }
          else{
            LOG(info)<<"DEBUG :: ERROR :: Not an ambiguous track either ::";
          }
        } // track doesn't have collision
        else                      {CollisionId_posError++;}
      }
      else { CollisionId_NoError++; hasCollision = true;}

      if( !hasCollision && !isAmbgTrack ) {LOG(info)<<"DEBUG :: ERROR :: A track with no collsiion and is not Ambiguous";}
      if(  hasCollision &&  isAmbgTrack ){ LOG(info)<<"DEBUG :: ERROR :: A track has collision and is also ambiguous";}

      if( hasCollision ) { LastTrackHadCollision = true;}
      doCollisionUpdate = false; //default is false;      
      doAmbgUpdate = false;
      if( hasCollision ){
        if(LastTrackHadCollision)
        { 
          if(collision.globalIndex() == oldCollisionIndex) { doCollisionUpdate = false;} //if collisions are same
          else                                             { doCollisionUpdate = true; } //if collisions are different
        }
        else { doCollisionUpdate = true; }//LastTrackWasAmbiguous
      }
      else if(isAmbgTrack){
        doAmbgUpdate = true;
        // if(LastTrackIsAmbg){
        //   if( haveSameInfo ) { doAmbgUpdate = false;}
        //   else              { doAmbgUpdate = true; }
        // }
        // else { doAmbgUpdate = true;} //Last track had Collisions
      }
      if( doAmbgUpdate ){ continue;}
      if( doCollisionUpdate || doAmbgUpdate) { //collision.globalIndex() != oldCollisionIndex){ //don't update if info is same as old collision
        if(doCollisionUpdate) {
          oldCollisionIndex = collision.globalIndex();
          bc = collision.bc_as<myBCTable>();
        }
        if(doAmbgUpdate){
          //do nothing
          // bc = collisions.iteratorAt(2).bc_as<myBCTable>();
          // bc = ambgTracks.iteratorAt(0).bc_as<myBCTable>();
        }
        // LOG(info)<<" What happens in the case when the collision id is = -1 and it tries to obtain bc"
        GetTimingInfo(bc, lastRun, nBCsPerTF, bcSOR, time, TFidThis, bcInTF);
      }

      recoEvent.fill(HIST("Event/track/Full/h_Pt" ),track.pt());
      recoEvent.fill(HIST("Event/track/Full/h_QPt"),track.signed1Pt());
    }
    PrintTime(StartTrackHist, Form("DEBUG :: df_%ld :: TrackHist Time :: ",dfCount));

    auto StartTrackDS = std::chrono::high_resolution_clock::now();

    uint64_t hEntriesTrkPt  = recoEvent.get<TH1>(HIST("Event/track/Full/h_Pt" ))->GetEntries();
    uint64_t hEntriesTrkQPt = recoEvent.get<TH1>(HIST("Event/track/Full/h_QPt"))->GetEntries();

    for(const auto& trackQA : tracksQA){
      // auto track = trackQA.track_as<myTracks>;  // dereferncing not working // only option is either use iterator way or Table slicing    
      
      auto const& track = tracksWithPid.iteratorAt(trackQA.trackId()); //Find the track in track table
      if( track.globalIndex() != trackQA.trackId()      ){ LOG(info)<<"DEBUG :: ERROR :: Track and TrackQA Mismatch";} 

      auto const& collision = collisions.iteratorAt(track.collisionId()); //Find the collision in collision table
      if( track.collisionId() != collision.globalIndex()){ LOG(info)<<"DEBUG :: ERROR :: track collId and collID Mismatch";}

      //Checking out of the range errors
      if( trackQA.trackId()   < 0 || tracks.size()     <= trackQA.trackId()  ) { LOG(info)<<"DEBUG :: ERROR :: trackQA has index out of scope :: trackQA.trackId() = "<<trackQA.trackId()<<" :: track.collisionId() = "<<track.collisionId()<<" :: track.signed1Pt() = "<<track.signed1Pt();}
      hasCollision = false;
      isAmbgTrack = false;
      if( track.collisionId() < 0 || collisions.size() <= track.collisionId()) {
        // LOG(info)<<"DEBUG :: ERROR :: track   has index out of scope :: trackQA.trackId() = "<<trackQA.trackId()<<" :: track.collisionId() = "<<track.collisionId()<<" :: track.signed1Pt() = "<<track.signed1Pt();
        if(track.collisionId()< 0){CollisionId_negError++;
          //check if it is an ambiguous track
          int ambgPos = FindInTable(track.globalIndex(), ambgTracks);
          if( ambgPos >= 0 && (ambgTracks.iteratorAt(ambgPos).trackId() == track.globalIndex()) ){
            nAmbgTracks++;
            isAmbgTrack = true;
            // LOG(info)<<"DEBUG :: Track is an ambiguous Track :: ambgId = "<<ambgTracks.iteratorAt(ambgPos).trackId()<<" :: trackId = "<<track.globalIndex();
          }
          else{
            LOG(info)<<"DEBUG :: ERROR :: Not an ambiguous track either ::";
          }
        } // track doesn't have collision
        else                      {CollisionId_posError++;}
      }
      else { CollisionId_NoError++; hasCollision = true;}

      if( !hasCollision && !isAmbgTrack ) {LOG(info)<<"DEBUG :: ERROR :: A track with no collsiion and is not Ambiguous";}
      if(  hasCollision &&  isAmbgTrack ){ LOG(info)<<"DEBUG :: ERROR :: A track has collision and is also ambiguous";}


      if( hasCollision ) { LastTrackHadCollision = true;}
      doCollisionUpdate = false; //default is false;      
      doAmbgUpdate = false;
      if( hasCollision ){
        if(LastTrackHadCollision)
        { 
          if(collision.globalIndex() == oldCollisionIndex) { doCollisionUpdate = false;} //if collisions are same
          else                                             { doCollisionUpdate = true; } //if collisions are different
        }
        else { doCollisionUpdate = true; }//LastTrackWasAmbiguous
      }
      else if(isAmbgTrack){
        doAmbgUpdate = true;
        // if(LastTrackIsAmbg){
        //   if( haveSameInfo ) { doAmbgUpdate = false;}
        //   else              { doAmbgUpdate = true; }
        // }
        // else { doAmbgUpdate = true;} //Last track had Collisions
      }
      if( doAmbgUpdate ){ continue;}
      if( doCollisionUpdate || doAmbgUpdate) { //collision.globalIndex() != oldCollisionIndex){ //don't update if info is same as old collision
        if(doCollisionUpdate) {
          oldCollisionIndex = collision.globalIndex();
          bc = collision.bc_as<myBCTable>();
        }
        if(doAmbgUpdate){
          //do nothing
          // bc = collisions.iteratorAt(2).bc_as<myBCTable>();
          // bc = ambgTracks.iteratorAt(0).bc_as<myBCTable>();
        }
        // LOG(info)<<" What happens in the case when the collision id is = -1 and it tries to obtain bc"
      GetTimingInfo(bc, lastRun, nBCsPerTF, bcSOR, time, TFidThis, bcInTF);
      }

      float weightDSPt = -1.0;
      float weightDSQPt = -1.0;
      int  triggerMaskDS = 0;

      bool isDownsampledMB  = false;
      bool isDownsampledPT  = false;
      bool isDownsampledQPT = false;
      bool isDownsampledV0  = false;
      bool isDownsampledPhi = false;
      bool isDownsampledD0  = false;

      float probRatioPt = -999;
      float probRatioQPT = -999;
      
      if(cfgDoSkipDF && (dfCount < cfgNSkipDF)) {continue;} // what if data frames are empty // In data you did encountered empty df's

      // //check non zero count in the bins for ratio
      if(startDownsampling == false){
        int count1 = recoEvent.get<TH1>(HIST("Event/track/Full/h_Pt" ))->GetBinContent(hPt_BinFor1);
        int count2 = recoEvent.get<TH1>(HIST("Event/track/Full/h_QPt"))->GetBinContent(hQPt_BinFor_Plus1) + recoEvent.get<TH1>(HIST("Event/track/Full/h_QPt"))->GetBinContent(hQPt_BinFor_Minus1);
        if ( count1 > 0  && count2 > 0 ) {
          startDownsampling = true;
        }
      }
      if(startDownsampling == false){
        recoEvent.fill(HIST("Event/track/Rejected/h_Pt" ), track.pt());
        recoEvent.fill(HIST("Event/track/Rejected/h_QPt"), track.signed1Pt());
        continue;}

      isDownsampledMB  = gRandom->Rndm()<factorMB;
      isDownsampledPT  =  checkTriggerPt(track.pt()       , factorPt, weightDSPt , probRatioPt , recoEvent.get<TH1>(HIST("Event/track/Full/h_Pt" )));
      isDownsampledQPT = checkTriggerQPt(track.signed1Pt(), factorPt, weightDSQPt, probRatioQPT, recoEvent.get<TH1>(HIST("Event/track/Full/h_QPt")));

      // check If It is a V0 daughter and getItsV0s. 
      if(track.sign() > 0) { isDownsampledV0 = std::binary_search(posV0DauList.begin(), posV0DauList.end(), track.globalIndex());}
      else                 { isDownsampledV0 = std::binary_search(negV0DauList.begin(), negV0DauList.end(), track.globalIndex());}

      // check if It is a prong of Phi
      if(track.sign() > 0) { isDownsampledPhi = std::binary_search(posKaonPhiList.begin(), posKaonPhiList.end(), track.globalIndex());}
      else                 { isDownsampledPhi = std::binary_search(negKaonPhiList.begin(), negKaonPhiList.end(), track.globalIndex());}

      // check If it is a prong of D0
      if      (std::binary_search(prong0D0List.begin(), prong0D0List.end(), track.globalIndex())){ isDownsampledD0 = true;} 
      else if (std::binary_search(prong1D0List.begin(), prong1D0List.end(), track.globalIndex())){ isDownsampledD0 = true;}

      triggerMaskDS = (1<<0)*isDownsampledMB + (1<<1)*isDownsampledPT + (1<<2)*isDownsampledQPT + (1<<3)*isDownsampledV0 + (1<<4)*isDownsampledPhi + (1<<5)*isDownsampledD0;
      
      recoEvent.fill(HIST("Event/track/isDownsampledMB" ), isDownsampledMB );  
      recoEvent.fill(HIST("Event/track/isDownsampledPT" ), isDownsampledPT );  
      recoEvent.fill(HIST("Event/track/isDownsampledQPT"), isDownsampledQPT); 
      recoEvent.fill(HIST("Event/track/FullMask"        ), triggerMaskDS);

      recoEvent.fill(HIST("Event/track/probRatioPt" ), probRatioPt );
      recoEvent.fill(HIST("Event/track/probRatioQPT"), probRatioQPT);

      if( cfgDoDSOfTrack && (triggerMaskDS == 0)) {
        recoEvent.fill(HIST("Event/track/Rejected/h_Pt" ), track.pt());
        recoEvent.fill(HIST("Event/track/Rejected/h_QPt"), track.signed1Pt());
        trackSkippedBecauseOfMask++;
        continue;}

      recoEvent.fill(HIST("Event/track/Accepted/h_Pt" ), track.pt());
      recoEvent.fill(HIST("Event/track/Accepted/h_QPt"), track.signed1Pt());

      tpcTime0           = trackQA.tpcTime0          ();
      tpcdcaR            = trackQA.tpcdcaR           ();
      tpcdcaZ            = trackQA.tpcdcaZ           ();
      tpcClusterByteMask = trackQA.tpcClusterByteMask();
      tpcdEdxMax0R       = trackQA.tpcdEdxMax0R      ();
      tpcdEdxMax1R       = trackQA.tpcdEdxMax1R      ();
      tpcdEdxMax2R       = trackQA.tpcdEdxMax2R      ();
      tpcdEdxMax3R       = trackQA.tpcdEdxMax3R      ();
      tpcdEdxTot0R       = trackQA.tpcdEdxTot0R      ();
      tpcdEdxTot1R       = trackQA.tpcdEdxTot1R      ();
      tpcdEdxTot2R       = trackQA.tpcdEdxTot2R      ();
      tpcdEdxTot3R       = trackQA.tpcdEdxTot3R      ();


      deltaRefContParamY    = trackQA.deltaRefContParamY(); 
      deltaRefContParamZ    = trackQA.deltaRefITSParamZ();
      deltaRefContParamSnp  = trackQA.deltaRefContParamSnp();
      deltaRefContParamTgl  = trackQA.deltaRefContParamTgl();
      deltaRefContParamQ2Pt = trackQA.deltaRefContParamQ2Pt();
      deltaRefGloParamY     = trackQA.deltaRefGloParamY();
      deltaRefGloParamZ     = trackQA.deltaRefGloParamZ();
      deltaRefGloParamSnp   = trackQA.deltaRefGloParamSnp();
      deltaRefGloParamTgl   = trackQA.deltaRefGloParamTgl();
      deltaRefGloParamQ2Pt  = trackQA.deltaRefGloParamQ2Pt();
      deltaTOFdX            = trackQA.deltaTOFdX();
      deltaTOFdZ            = trackQA.deltaTOFdZ();

      trackQAcounter++;
      trackGlobalIndexList.push_back(track.globalIndex());
      DrTrackPositionIndexList.push_back(trackQAcounter);
      
      auto const& trackMeanOcc = trackMeanOccs.iteratorAt(trackQAcounter+trackSkippedBecauseOfMask); //(trackQA.globalIndex()); //because of ambgTracks, you cant use globalIndex()
      if( trackMeanOcc.trackId() != trackQA.trackId()   ){ LOG(info)<<"DEBUG :: ERROR :: TrackQA and trackMeanOcc mismatch";}

      uint8_t detectorMask = (1<<0)*track.hasITS() + (1<<1)*track.hasTPC() + (1<<2)*track.hasTRD() + (1<<3)*track.hasTOF();

      //Store Track Data in trees
      DerivedTracks(
         trackQAcounter+gTrackChecker     //GlobalIndex
        ,track.globalIndex()                   //OrignalIndex
        ,bc.globalBC()                         //GlobalBC
        ,track.collisionId()+gCollisionChecker //GCollId
        ,dfCount
        ,track.collisionId()                   //DrCollisionId
        ,track.collisionId() 	      // Original Collision to which this track belongs
        ,track.trackType() 	        // Type of track. See enum TrackTypeEnum. This cannot be used to decide which detector has contributed to this track. Use hasITS, hasTPC, etc.
        ,track.x() 	                // 
        ,track.alpha() 	            // 
        ,track.y() 	                // 
        ,track.z() 	                // 
        ,track.snp() 	              // 
        ,track.tgl() 	              // 
        ,track.signed1Pt() 	        // (sign of charge)/Pt in c/GeV. Use pt() and sign() instead
        ,track.pt() 	              // Transverse momentum of the track in GeV/c
        ,track.p() 	                // Momentum in Gev/c
        ,track.eta() 	              // Pseudorapidity
        ,track.phi() 	              // Phi of the track, in radians within [0, 2pi)

        ,track.isWithinBeamPipe() 	  // Is the track within the beam pipe (= successfully propagated to a collision vertex)
        ,track.px() 	                // Momentum in x-direction in GeV/c
        ,track.py() 	                // Momentum in y-direction in GeV/c
        ,track.pz() 	                // Momentum in z-direction in GeV/c
        // ,drTrack::PVector         // D 	pVector             std::array<float,3> 	Momentum vector in x,y,z-directions in GeV/c
        ,track.energy(0.13957) 	            // Track energy, computed under the mass assumption given as input
        ,track.rapidity(0.13957) 	          // Track rapidity, computed under the mass assumption given as input
        ,track.sign() 	              // Charge: positive: 1, negative: -1

        // // TracksCov
        ,track.sigmaY() 	    //Covariance matrix
        ,track.sigmaZ() 	    //Covariance matrix
        ,track.sigmaSnp() 	  //Covariance matrix
        ,track.sigmaTgl() 	  //Covariance matrix
        ,track.sigma1Pt() 	  //Covariance matrix
        ,track.rhoZY() 	    //Covariance matrix in compressed form
        ,track.rhoSnpY() 	 	//Covariance matrix in compressed form
        ,track.rhoSnpZ() 	 	//Covariance matrix in compressed form
        ,track.rhoTglY() 	 	//Covariance matrix in compressed form
        ,track.rhoTglZ() 	 	//Covariance matrix in compressed form
        ,track.rhoTglSnp() 	//Covariance matrix in compressed form
        ,track.rho1PtY() 	 	//Covariance matrix in compressed form
        ,track.rho1PtZ() 	 	//Covariance matrix in compressed form
        ,track.rho1PtSnp() 	//Covariance matrix in compressed form
        ,track.rho1PtTgl() 	//Covariance matrix in compressed form
        ,track.cYY() 	  //
        ,track.cZY() 	  
        ,track.cZZ() 	  
        ,track.cSnpY() 	  
        ,track.cSnpZ() 	  
        ,track.cSnpSnp()
        ,track.cTglY() 	  
        ,track.cTglZ() 	  
        ,track.cTglSnp()
        ,track.cTglTgl()
        ,track.c1PtY() 	  
        ,track.c1PtZ() 	  
        ,track.c1PtSnp()
        ,track.c1PtTgl()
        ,track.c1Pt21Pt2()

        // // TracksExtra
        ,track.tpcInnerParam() 	                    // Momentum at inner wall of the TPC
        ,track.flags() 	                            // Track flags. Run 2: see TrackFlagsRun2Enum | Run 3: see TrackFlags
        ,track.itsClusterSizes() 	                  // Clusters sizes, four bits per a layer, starting from the innermost
        ,track.tpcNClsFindable() 	                  // Findable TPC clusters for this track geometry
        ,track.tpcNClsFindableMinusFound() 	        // TPC Clusters: Findable - Found
        ,track.tpcNClsFindableMinusCrossedRows() 	  // TPC Clusters: Findable - crossed rows
        ,track.tpcNClsShared() 	                    // Number of shared TPC clusters
        ,track.trdPattern() 	                      // Contributor to the track on TRD layer in bits 0-5, starting from the innermost, bit 6 indicates a potentially split tracklet, bit 7 if the track crossed a padrow
        ,track.itsChi2NCl() 	                      // Chi2 / cluster for the ITS track segment
        ,track.tpcChi2NCl() 	                      // Chi2 / cluster for the TPC track segment
        ,track.trdChi2() 	                          // Chi2 for the TRD track segment
        ,track.tofChi2() 	                          // Chi2 for the TOF track segment
        ,track.tpcSignal() 	                        // dE/dx signal in the TPC
        ,track.trdSignal() 	                        // PID signal in the TRD
        ,track.length() 	                          // Track length
        ,track.tofExpMom() 	                        // TOF expected momentum obtained in tracking, used to compute the expected times
        ,track.trackEtaEmcal() 	                    // 
        ,track.trackPhiEmcal() 	                    // 
        ,track.trackTime() 	                        // Estimated time of the track in ns wrt collision().bc() or ambiguoustrack.bcSlice()[0]
        ,track.trackTimeRes() 	                    // Resolution of the track time in ns (see TrackFlags::TrackTimeResIsRange)
        
        // ,track.tpcDeltaTFwd() 	                    // Delta Forward of track time in TPC time bis
        // ,track.tpcDeltaTBwd() 	                    // Delta Backward of track time in TPC time bis
        ,track.pidForTracking() 	                  // PID hypothesis used during tracking. See the constants in the class PID in PID.h
        ,track.isPVContributor() 	                  // Run 3: Has this track contributed to the collision vertex fit

        ,detectorMask
        // ,track.hasITS() 	                          // Flag to check if track has a ITS match
        // ,track.hasTPC() 	                          // Flag to check if track has a TPC match
        // ,track.hasTRD() 	                          // Flag to check if track has a TRD match
        // ,track.hasTOF() 	                          // Flag to check if track has a TOF measurement

        ,track.tpcNClsFound() 	                    // Number of found TPC clusters
        ,track.tpcNClsCrossedRows() 	              // Number of crossed TPC Rows
        ,track.itsClusterMap() 	                    // ITS cluster map, one bit per a layer, starting from the innermost
        ,track.itsNCls() 	                          // Number of ITS clusters
        ,track.itsNClsInnerBarrel() 	              // Number of ITS clusters in the Inner Barrel
        // ,track.itsClsSizeInLayer() 	                // Size of the ITS cluster in a given layer
        // ,track.isITSAfterburner() 	                // If the track used the afterburner in the ITS
        // ,track.tofExpTime() 	                      // Expected time for the track to reach the TOF
        // ,track.tofExpTimePi() 	                    // Expected time for the track to reach the TOF under the pion hypothesis
        // ,track.tofExpTimeKa() 	                    // Expected time for the track to reach the TOF under the kaon hypothesis
        // ,track.tofExpTimePr() 	                    // Expected time for the track to reach the TOF under the proton hypothesis
        ,track.tpcCrossedRowsOverFindableCls() 	    // Ratio crossed rows over findable clusters
        ,track.tpcFoundOverFindableCls()            // Ratio of found over findable clusters
        ,track.tpcFractionSharedCls() 	            // Fraction of shared TPC clusters
        ,track.detectorMap() 	                      // Detector map version 1, see enum DetectorMapEnum
      
        // TracksQA
        // // ,o2::soa::Index 	GI 	globalIndex 	    int64_t 	
        // // ,o2::aod::trackqa::TrackId 	            //I 	    trackId 	        //int32 	track to which this QA information belongs
        ,tpcTime0             //tpc only time0 (mTime0 in TPC track)
        ,tpcdcaR              //tpc only DCAr
        ,tpcdcaZ              //tpc only DCAz
        ,tpcClusterByteMask   //tracklet bitmask - track defining 8 tracklets (152=8*19 rows) bit set if nCluster>thr (default 5)
        ,tpcdEdxMax0R         //TPC dEdxQMax -ROC0/dEdx
        ,tpcdEdxMax1R         //TPC dEdxQMax -ROC1/dEdx
        ,tpcdEdxMax2R         //TPC dEdxQMax -ROC2/dEdx
        ,tpcdEdxMax3R         //TPC dEdxQMax -ROC3/dEdx
        ,tpcdEdxTot0R         //TPC dEdxQtot -ROC0/dEdx
        ,tpcdEdxTot1R         //TPC dEdxQtot -ROC1/dEdx
        ,tpcdEdxTot2R         //TPC dEdxQtot -ROC2/dEdx
        ,tpcdEdxTot3R         //TPC dEdxQtot -ROC3/dEdx
        //
        ,track.beta()
        ,track.betaerror()
        ,track.mass()

        ,track.tpcNSigmaPi()
        ,track.tofNSigmaPi()
        ,track.tpcNSigmaKa()
        ,track.tofNSigmaKa()
        ,track.tpcNSigmaPr()
        ,track.tofNSigmaPr()
        ,track.tpcNSigmaEl()
        ,track.tofNSigmaEl()
        ,track.tpcNSigmaDe()
        ,track.tofNSigmaDe()
        
        ,track.tofSignal()
        // ,track.eventCollisionTime()
        ,track.goodTOFMatch()

        ,track.tofFlags()        // o2::aod::pidflags::TOFFlags		     //  tofFlags	uint8_t	Flag for the complementary TOF PID information for the event time
        ,track.isEvTimeDefined() //  D	isEvTimeDefined	bool	True if the Event Time was computed with any method i.e. there is a usable event time
        ,track.isEvTimeTOF()	   //  D	isEvTimeTOF	bool	True if the Event Time was computed with the TOF
        ,track.isEvTimeT0AC()	   //  D	isEvTimeT0AC	bool	True if the Event Time was computed with the T0AC
        ,track.isEvTimeTOFT0AC() //  

        ,time
        ,TFidThis
        ,bcInTF

        ,trackMeanOcc.trackId()

        ,trackMeanOcc.tmoPrimUnfm80() 
        ,trackMeanOcc.tmoFV0AUnfm80() 
        ,trackMeanOcc.tmoFV0CUnfm80() 
        ,trackMeanOcc.tmoFT0AUnfm80() 
        ,trackMeanOcc.tmoFT0CUnfm80() 
        ,trackMeanOcc.tmoFDDAUnfm80() 
        ,trackMeanOcc.tmoFDDCUnfm80() 
        ,trackMeanOcc.tmoNTrackITSUnfm80()   
        ,trackMeanOcc.tmoNTrackTPCUnfm80()   
        ,trackMeanOcc.tmoNTrackTRDUnfm80()   
        ,trackMeanOcc.tmoNTrackTOFUnfm80()   
        ,trackMeanOcc.tmoNTrackSizeUnfm80()   
        ,trackMeanOcc.tmoNTrackTPCAUnfm80()  
        ,trackMeanOcc.tmoNTrackTPCCUnfm80()  
        ,trackMeanOcc.tmoNTrackITSTPCUnfm80()
        ,trackMeanOcc.tmoNTrackITSTPCAUnfm80()
        ,trackMeanOcc.tmoNTrackITSTPCCUnfm80()
        ,trackMeanOcc.tmoMultNTracksHasITSUnfm80()
        ,trackMeanOcc.tmoMultNTracksHasTPCUnfm80()
        ,trackMeanOcc.tmoMultNTracksHasTOFUnfm80()
        ,trackMeanOcc.tmoMultNTracksHasTRDUnfm80()
        ,trackMeanOcc.tmoMultNTracksITSOnlyUnfm80()
        ,trackMeanOcc.tmoMultNTracksTPCOnlyUnfm80()
        ,trackMeanOcc.tmoMultNTracksITSTPCUnfm80()
        ,trackMeanOcc.tmoMultAllTracksTPCOnlyUnfm80()
        ,trackMeanOcc.tmoRobustT0V0PrimUnfm80()   
        ,trackMeanOcc.tmoRobustFDDT0V0PrimUnfm80()
        ,trackMeanOcc.tmoRobustNtrackDetUnfm80()  
        ,trackMeanOcc.tmoRobustMultExtraTableUnfm80()  
        ,trackMeanOcc.twmoPrimUnfm80() 
        ,trackMeanOcc.twmoFV0AUnfm80() 
        ,trackMeanOcc.twmoFV0CUnfm80() 
        ,trackMeanOcc.twmoFT0AUnfm80() 
        ,trackMeanOcc.twmoFT0CUnfm80() 
        ,trackMeanOcc.twmoFDDAUnfm80() 
        ,trackMeanOcc.twmoFDDCUnfm80() 
        ,trackMeanOcc.twmoNTrackITSUnfm80()   
        ,trackMeanOcc.twmoNTrackTPCUnfm80()   
        ,trackMeanOcc.twmoNTrackTRDUnfm80()   
        ,trackMeanOcc.twmoNTrackTOFUnfm80()   
        ,trackMeanOcc.twmoNTrackSizeUnfm80()   
        ,trackMeanOcc.twmoNTrackTPCAUnfm80()  
        ,trackMeanOcc.twmoNTrackTPCCUnfm80()  
        ,trackMeanOcc.twmoNTrackITSTPCUnfm80()
        ,trackMeanOcc.twmoNTrackITSTPCAUnfm80()
        ,trackMeanOcc.twmoNTrackITSTPCCUnfm80()
        ,trackMeanOcc.twmoMultNTracksHasITSUnfm80()
        ,trackMeanOcc.twmoMultNTracksHasTPCUnfm80()
        ,trackMeanOcc.twmoMultNTracksHasTOFUnfm80()
        ,trackMeanOcc.twmoMultNTracksHasTRDUnfm80()
        ,trackMeanOcc.twmoMultNTracksITSOnlyUnfm80()
        ,trackMeanOcc.twmoMultNTracksTPCOnlyUnfm80()
        ,trackMeanOcc.twmoMultNTracksITSTPCUnfm80()
        ,trackMeanOcc.twmoMultAllTracksTPCOnlyUnfm80()
        ,trackMeanOcc.twmoRobustT0V0PrimUnfm80()   
        ,trackMeanOcc.twmoRobustFDDT0V0PrimUnfm80()
        ,trackMeanOcc.twmoRobustNtrackDetUnfm80()  
        ,trackMeanOcc.twmoRobustMultExtraTableUnfm80()
        
        ,track.dcaXY()   //	dcaXY	float	Impact parameter in XY of the track to the primary vertex
        ,track.dcaZ()	   //	dcaZ	float	Impact parameter in Z of the track to the primary vertex
        ,track.sigmaDcaXY2()  //		sigmaDcaXY2	float	Impact parameter sigma^2 in XY of the track to the primary vertex
        ,track.sigmaDcaZ2()   //		sigmaDcaZ2	float	Impact parameter sigma^2 in Z of the track to the primary verte

        ,track.usedForTOFEvTime()
        ,track.evTimeTOF()
        ,track.evTimeTOFErr()
        ,track.evTimeTOFMult()

        ,track.tofExpTimeEl() //D	tofExpTimeEl	float	Expected time for the track to reach the TOF under the electron hypothesis
        ,track.tofExpTimeMu() //D	tofExpTimeMu	float	Expected time for the track to reach the TOF under the muon hypothesis
        ,track.tofExpTimePi() //D	tofExpTimePi	float	Expected time for the track to reach the TOF under the pion hypothesis
        ,track.tofExpTimeKa() //D	tofExpTimeKa	float	Expected time for the track to reach the TOF under the kaon hypothesis
        ,track.tofExpTimePr() //D	tofExpTimePr	float	Expected time for the track to reach the TOF under the proton hypothesis
        ,track.tofExpTimeDe() //D	tofExpTimeDe	float	Expected time for the track to reach the TOF under the deuteron hypothesis
        ,track.tofExpTimeTr() //D	tofExpTimeTr	float	Expected time for the track to reach the TOF under the triton hypothesis
        ,track.tofExpTimeHe() //D	tofExpTimeHe	float	Expected time for the track to reach the TOF under the helium3 hypothesis
        ,track.tofExpTimeAl() //D	tofExpTimeAl

        //o2::aod::track::v001::ITSClsSizeInLayer	D	itsClsSizeInLayer	uint8_t	Size of the ITS cluster in a given layer
        // ,track.itsClsSizeInLayer() //	uint8_t	Size of the ITS cluster in a given layer

        ,track.itsNSigmaEl()
        ,track.itsNSigmaMu()
        ,track.itsNSigmaPi()
        ,track.itsNSigmaKa()
        ,track.itsNSigmaPr()
        ,track.itsNSigmaDe()
        ,track.itsNSigmaTr()
        ,track.itsNSigmaHe()
        ,track.itsNSigmaAl()

        ,deltaRefContParamY   
        ,deltaRefContParamZ   
        ,deltaRefContParamSnp 
        ,deltaRefContParamTgl 
        ,deltaRefContParamQ2Pt
        ,deltaRefGloParamY    
        ,deltaRefGloParamZ    
        ,deltaRefGloParamSnp  
        ,deltaRefGloParamTgl  
        ,deltaRefGloParamQ2Pt 
        ,deltaTOFdX
        ,deltaTOFdZ           

        ,collision.centFV0A()
        ,collision.centFT0M()
        ,collision.centFT0A()
        ,collision.centFT0C()
        ,collision.centFDDM()
        ,collision.centNTPV()
        ,collision.trackOccupancyInTimeRange()  //	int	Occupancy in specified time interval by a number of tracks from nearby collisions
        ,collision.ft0cOccupancyInTimeRange()  //	float	Occupancy in specified time interval by a sum of FT0C amplitudes from nearby collisions

        ,weightDSPt
        ,weightDSQPt
        ,triggerMaskDS
        ,hEntriesTrkPt 
        ,hEntriesTrkQPt

        ,histPDFPol1Evaluation(recoEvent.get<TH1>(HIST("Event/track/Full/h_Pt" )), track.pt()       , 2)
        ,histPDFPol1Evaluation(recoEvent.get<TH1>(HIST("Event/track/Full/h_QPt")), track.signed1Pt(), 2)

      );
      trackCounter++;
    }// end of trackQA loop

    PrintTime(StartTrackDS, Form("DEBUG :: df_%ld :: TrackTable Time :: ",dfCount));

    auto StartV0Fill = std::chrono::high_resolution_clock::now();

    int lastCollId = -10;
    int v0Counter = 0;

    int PosIndexPosition = -1;
    int NegIndexPosition = -1;
    int indexPosition = -1;
    std::vector<std::vector<int64_t>> v0PositionList;
    uint64_t hEntriesV0Pt              = recoEvent.get<TH1>(HIST("Event/V0/h_Pt"             ))->GetEntries();
    uint64_t hEntriesLambdaPt          = recoEvent.get<TH1>(HIST("Event/Lambda/h_Pt"         ))->GetEntries();
    uint64_t hEntriesAntiLambdaPt      = recoEvent.get<TH1>(HIST("Event/AntiLambda/h_Pt"     ))->GetEntries();
    uint64_t hEntriesK0ShortPt         = recoEvent.get<TH1>(HIST("Event/K0Short/h_Pt"        ))->GetEntries();
    // uint64_t hEntriesGammaPt           = recoEvent.get<TH1>(HIST("Event/Gamma/h_Pt"          ))->GetEntries();
    // uint64_t hEntriesHypertritonPt     = recoEvent.get<TH1>(HIST("Event/Hypertriton/h_Pt"    ))->GetEntries();
    // uint64_t hEntriesAntiHypertritonPt = recoEvent.get<TH1>(HIST("Event/AntiHypertriton/h_Pt"))->GetEntries();

    for (const auto& v0 : V0s){
      auto myV0Coll = collisions.iteratorAt(v0.collisionId());
      if(v0.collisionId() != collisions.iteratorAt(v0.collisionId()).globalIndex()) {
        LOG(info)<<"DEBUG :: ERROR :: iteratorAt for collision dereferencing not working";
        continue;
      }
      if( v0.collisionId() != lastCollId){
        lastCollId = v0.collisionId();
        bc = myV0Coll.bc_as<myBCTable>();
        GetTimingInfo(bc, lastRun, nBCsPerTF, bcSOR, time, TFidThis, bcInTF);
      };
      
      if(!std::binary_search(trigV0IndexList.begin(), trigV0IndexList.end(), v0.globalIndex())){ continue;}
      if(trackGlobalIndexList.size() > 0 ){
        PosIndexPosition = FindTrackIdInList(v0.posTrackId(),trackGlobalIndexList, DrTrackPositionIndexList, trackGlobalIndexList.size()); //Position in created DrTrackTable per DataFrame 
        NegIndexPosition = FindTrackIdInList(v0.negTrackId(),trackGlobalIndexList, DrTrackPositionIndexList, trackGlobalIndexList.size()); //Position in created DrTrackTable per DataFrame 
        if(PosIndexPosition < 0 && NegIndexPosition < 0) { continue;}
      }
      v0PositionList.push_back({v0.globalIndex(), static_cast<int64_t>(v0Counter+gV0Checker)});
      indexPosition = BinarySearchVector(v0.globalIndex(), trigV0IndexList, 0, trigV0IndexList.size()-1);

      DerivedV0s(
        v0Counter+gV0Checker                //fGlobalIndex
        ,v0.globalIndex()                   //fOriginalIndex
        ,bc.globalBC()                      //fGBCid
        ,dfCount
        ,v0.collisionId()+gCollisionChecker //fGCollId
        ,v0Counter                          //fGV0Id
        ,PosIndexPosition+gTrackChecker     //fGPosTrackId       //Position in created DrTrackTable in merged DataFrame
        ,NegIndexPosition+gTrackChecker     //fGNegTrackId       //Position in created DrTrackTable in merged DataFrame
        ,PosIndexPosition                   //fIndexDrTracks_Pos //Position in created DrTrackTable per DataFrame
        ,NegIndexPosition                   //fIndexDrTracks_Neg //Position in created DrTrackTable per DataFrame
        ,v0.posTrackId() 	                  //fIndexTracks_Pos   //Position in original table//int 	Pointer into Tracks
        ,v0.negTrackId() 	                  //fIndexTracks_Neg   //Position in original table//int 	Pointer into Tracks
        ,v0.collisionId() 	//int32 	Pointer into Collisions
        ,v0.v0Id() 	        //int32 	Pointer into V0s
        ,v0.posX() 	        //float 	positive track X at min
        ,v0.negX() 	        //float 	negative track X at min
        ,v0.x() 	          //float 	decay position X
        ,v0.y() 	          //float 	decay position Y
        ,v0.z() 	          //float 	decay position Z
        ,v0.dcaV0daughters() 	//float 	DCA between V0 daughters
        ,v0.dcapostopv() 	    //float 	DCA positive prong to PV
        ,v0.dcanegtopv() 	    //float 	DCA negative prong to PV
        ,v0.v0cosPA() 	      //float 	V0 CosPA
        ,v0.dcav0topv() 	    //float 	DCA V0 to PV
        ,v0.v0radius() 	       //float 	V0 decay radius (2D, centered at zero)
        ,v0.mLambda() 	       //float 	mass under lambda hypothesis
        ,v0.mAntiLambda() 	   //float 	mass under antilambda hypothesis
        ,v0.mK0Short() 	       //float 	mass under K0short hypothesis
        ,v0.mGamma() 	         //float 	mass under gamma hypothesis
        ,v0.mHypertriton() 	   //float 	mass under hypertriton hypothesis
        ,v0.mAntiHypertriton() //float 	mass under antihypertriton hypothesis
        ,v0.negativept() 	     //float 	negative daughter pT
        ,v0.positivept() 	     //float 	positive daughter pT
        ,v0.negativeeta() 	   //float 	negative daughter eta
        ,v0.negativephi() 	   //float 	negative daughter phi
        ,v0.positiveeta() 	   //float 	positive daughter eta
        ,v0.positivephi() 	   //float 	positive daughter phi
        ,v0.isStandardV0() 	   //bool 	is standard V0
        ,v0.isPhotonTPConly()  //bool 	is tpc-only photon V0
        ,time
        ,TFidThis
        ,bcInTF
        ,v0.px()
        ,v0.py()
        ,v0.pz()        
        ,v0.pt()
        ,v0.eta()
        ,v0.phi()
        ,v0.alpha()
        ,v0.qtarm()
        // ,v0.distovertotmom()
        ,v0.psipair()
        ,v0.pfracpos()
        ,v0.pfracneg()
        
        ,weightDSPtV0[indexPosition]
        ,weightDSPtLambda[indexPosition]
        ,weightDSPtAntiLambda[indexPosition]
        ,weightDSPtK0Short[indexPosition]
        ,weightDSPtGamma[indexPosition]
        ,weightDSPtHypertriton[indexPosition]
        ,weightDSPtAntiHypertriton[indexPosition]
        ,vecTriggerMaskV0[indexPosition]
        ,hEntriesV0Pt             
        ,hEntriesLambdaPt         
        ,hEntriesAntiLambdaPt     
        ,hEntriesK0ShortPt        
        // ,hEntriesGammaPt          
        // ,hEntriesHypertritonPt    
        // ,hEntriesAntiHypertritonPt

        ,histPDFPol1Evaluation(recoEvent.get<TH1>(HIST("Event/V0/h_Pt"             )), v0.pt(), 2)
        ,histPDFPol1Evaluation(recoEvent.get<TH1>(HIST("Event/Lambda/h_Pt"         )), v0.pt(), 2)
        ,histPDFPol1Evaluation(recoEvent.get<TH1>(HIST("Event/AntiLambda/h_Pt"     )), v0.pt(), 2)
        ,histPDFPol1Evaluation(recoEvent.get<TH1>(HIST("Event/K0Short/h_Pt"        )), v0.pt(), 2)

      );
      v0Counter++;
    }
    PrintTime(StartV0Fill, Form("DEBUG :: df_%ld :: V0 Table Time :: ",dfCount));

    // ---------------------- Phi ------------------------ //
    auto StartPhiFill = std::chrono::high_resolution_clock::now();
    lastCollId = -10;
    int phiCounter = 0;

    int posKaonPhiIndexPosition = -1;
    int negKaonPhiIndexPosition = -1;

    std::vector<std::vector<int64_t>> phiPositionList;

    uint64_t hEntriesPhiPt     = recoEvent.get<TH1>(HIST("Event/Phi/h_Pt" ))->GetEntries();
    float rapidity = 0.0;

    for (const auto& phiCand : phiOIlist) {
      const auto& posTrack = tracks.rawIteratorAt(posKaOIlist[phiCand]);
      const auto& negTrack = tracks.rawIteratorAt(negKaOIlist[phiCand]);

      auto myPhiColl = collisions.iteratorAt(phiCollisionIdlist[phiCand]);
      if(phiCollisionIdlist[phiCand] != myPhiColl.globalIndex()) {
        LOG(info)<<"DEBUG :: ERROR :: iteratorAt for collision dereferencing not working for D0";
        continue;
      }

      if (phiCollisionIdlist[phiCand] != lastCollId) {
        lastCollId = phiCollisionIdlist[phiCand];
        bc = myPhiColl.bc_as<myBCTable>();
        GetTimingInfo(bc, lastRun, nBCsPerTF, bcSOR, time, TFidThis, bcInTF);
      }

      if (!std::binary_search(trigPhiIndexList.begin(), trigPhiIndexList.end(), phiCand)){continue;} //check if it is a triggered D0;
      if ( trackGlobalIndexList.size() > 0){
        posKaonPhiIndexPosition = FindTrackIdInList(posTrack.globalIndex(), trackGlobalIndexList, DrTrackPositionIndexList, trackGlobalIndexList.size());
        negKaonPhiIndexPosition = FindTrackIdInList(negTrack.globalIndex(), trackGlobalIndexList, DrTrackPositionIndexList, trackGlobalIndexList.size());
        if (posKaonPhiIndexPosition < 0 && negKaonPhiIndexPosition < 0) {continue;}
      }
      
      phiPositionList.push_back({phiCand, static_cast<int64_t>(phiCounter + gPhiChecker)});
      indexPosition = BinarySearchVector(phiCand, trigPhiIndexList, 0, trigPhiIndexList.size()-1);

      // if(indexPosition != phiCand) { LOG(error)<<"DEBUG :: ERROR ERROR ERROR :: Phi Index Mismatch Error";}

      mass1 = MassKaonCharged;
      mass2 = MassKaonCharged;
      pvec[0] = posTrack.px() + negTrack.px();
      pvec[1] = posTrack.py() + negTrack.py();
      pvec[2] = posTrack.pz() + negTrack.pz();
      pt = RecoDecay::pt(pvec[0], pvec[1]);
      eta = RecoDecay::eta(pvec);
      phi = RecoDecay::phi(pvec[0], pvec[1]);
      rapidity = RecoDecay::y(pvec, 1.019461);
      
      // get invariant mass of pair
      p = RecoDecay::p((posTrack.px() + negTrack.px()), (posTrack.py() + negTrack.py()), (posTrack.pz() + negTrack.pz()));
      e = RecoDecay::e(posTrack.px(), posTrack.py(), posTrack.pz(), mass1) + RecoDecay::e(negTrack.px(), negTrack.py(), negTrack.pz(), mass2);
      minv = std::sqrt(RecoDecay::m2(p, e));

      // // Fill Phi Derived Tables
      DerivedPhis(
        phiCounter+gPhiChecker
        ,phiCand
        ,bc.globalBC()
        ,dfCount
        ,phiCollisionIdlist[phiCand]+gCollisionChecker
        ,phiCounter

        ,posKaonPhiIndexPosition+gTrackChecker //global Id in merged DFs in derived Tracks
        ,negKaonPhiIndexPosition+gTrackChecker //global Id in merged DFs in derived Tracks

        ,posKaonPhiIndexPosition //global Id in single DF in derived Tracks
        ,negKaonPhiIndexPosition //global Id in single DF in derived Tracks

        ,posTrack.globalIndex() //original Id in single DF in original Tracks
        ,negTrack.globalIndex() //original Id in single DF in original Tracks

        ,phiCollisionIdlist[phiCand]

        ,pvec[0]    // phiPx
        ,pvec[1]    // phiPy
        ,pvec[2]    // phiPz

        ,pt
        ,eta
        ,phi
        ,rapidity
        ,minv
        ,e

        ,posTrack.px()
        ,posTrack.py()
        ,posTrack.pz()
        ,posTrack.pt()
        ,posTrack.eta()
        // ,posTrack.y() //(MassKaonCharged)
        // ,posTrack.isPVContributor()
        // ,posTrack.isGlobalTrack()

        ,negTrack.px()
        ,negTrack.py()
        ,negTrack.pz()
        ,negTrack.pt()
        ,negTrack.eta()
        // ,negTrack.y() //(MassKaonCharged)
        // ,negTrack.isPVContributor()
        // ,negTrack.isGlobalTrack()        

        ,time
        ,TFidThis
        ,bcInTF

        ,weightDSPtPhi[indexPosition]
        ,vecTriggerMaskPhi[indexPosition]
        ,hEntriesPhiPt
        ,histPDFPol1Evaluation(recoEvent.get<TH1>(HIST("Event/Phi/h_Pt" )), pt , 2)
      );

      phiCounter++;
    }

    PrintTime(StartPhiFill, Form("DEBUG :: df_%ld :: Phi Table Time :: ",dfCount));
    // ---------------------- D0 ------------------------ //
    auto StartD0Table = std::chrono::high_resolution_clock::now();
    lastCollId = -10;
    int d0Counter = 0;

    int Prong0IndexPosition = -1;
    int Prong1IndexPosition = -1;

    std::vector<std::vector<int64_t>> d0PositionList;

    uint64_t hEntriesD0Pt    = recoEvent.get<TH1>(HIST("Event/D0/h_Pt" ))  ->GetEntries();
    uint64_t hEntriesD0BarPt = recoEvent.get<TH1>(HIST("Event/D0Bar/h_Pt"))->GetEntries();

    for (const auto& d0 : D0s) {
      auto myD0Coll = collisions.iteratorAt(d0.collisionId());
      if(d0.collisionId() != collisions.iteratorAt(d0.collisionId()).globalIndex()) {
        LOG(info)<<"DEBUG :: ERROR :: iteratorAt for collision dereferencing not working for D0";
        continue;
      }

      if (d0.collisionId() != lastCollId) {
        lastCollId = d0.collisionId();
        bc = myD0Coll.bc_as<myBCTable>();
        GetTimingInfo(bc, lastRun, nBCsPerTF, bcSOR, time, TFidThis, bcInTF);
      }

      if (!std::binary_search(trigD0IndexList.begin(), trigD0IndexList.end(), d0.globalIndex())){continue;} //check if it is a triggered D0;
      if ( trackGlobalIndexList.size() > 0){
        Prong0IndexPosition = FindTrackIdInList(d0.prong0Id(), trackGlobalIndexList, DrTrackPositionIndexList, trackGlobalIndexList.size());
        Prong1IndexPosition = FindTrackIdInList(d0.prong1Id(), trackGlobalIndexList, DrTrackPositionIndexList, trackGlobalIndexList.size());
        if (Prong0IndexPosition < 0 && Prong1IndexPosition < 0) {continue;}
      }
      
      d0PositionList.push_back({d0.globalIndex(), static_cast<int64_t>(d0Counter + gD0Checker)});
      indexPosition = BinarySearchVector(d0.globalIndex(), trigD0IndexList, 0, trigD0IndexList.size()-1);

      // // Fill D0 Derived Tables
      DerivedD0s(
        d0Counter+gD0Checker
        ,d0.globalIndex()
        ,bc.globalBC()
        ,dfCount
        ,d0.collisionId()+gCollisionChecker
        ,d0Counter
        ,Prong0IndexPosition+gTrackChecker
        ,Prong1IndexPosition+gTrackChecker
        ,Prong0IndexPosition
        ,Prong1IndexPosition
        ,d0.prong0Id() 	      
        ,d0.prong1Id() 	      
        ,d0.collisionId()    
        ,d0.posX() 	          
        ,d0.posY() 	          
        ,d0.posZ() 	          
        ,d0.xSecondaryVertex()
        ,d0.ySecondaryVertex()
        ,d0.zSecondaryVertex()
        ,d0.errorDecayLength()
        ,d0.errorDecayLengthXY()
        ,d0.chi2PCA()
        ,d0.rSecondaryVertex() 	          //GI 		? 	
        ,d0.decayLength() 	
        ,d0.decayLengthXY() 	
        ,d0.decayLengthNormalised() 	
        ,d0.decayLengthXYNormalised() 	
        ,d0.impactParameterNormalised0() 	//GI 		? 	
        ,d0.ptProng0() 	
        ,d0.pt2Prong0() 	
        // ,d0.pVectorProng0()
        ,d0.impactParameterNormalised1() 	//GI 		? 	
        ,d0.ptProng1() 	
        ,d0.pt2Prong1() 	
        // ,d0.pVectorProng1()
        ,d0.pxProng0() 	
        ,d0.pyProng0() 	
        ,d0.pzProng0() 	
        ,d0.pxProng1() 	
        ,d0.pyProng1() 	
        ,d0.pzProng1() 	
        ,d0.impactParameter0() 	
        ,d0.impactParameter1() 	
        ,d0.errorImpactParameter0() 	
        ,d0.errorImpactParameter1() 	
        ,d0.impactParameterZ0() 	
        ,d0.impactParameterZ1() 	
        ,d0.errorImpactParameterZ0() 	
        ,d0.errorImpactParameterZ1() 	
        ,d0.nProngsContributorsPV()
        ,d0.hfflag()
        // ,d0.m() 	                          //    GI 		? 	
        // ,d0.m2() 	
        ,d0.impactParameterProduct() 	
        // ,d0.cosThetaStar() 	
        ,d0.impactParameterProngSqSum() 	
        ,d0.pt() 	                                //GI 		? 	
        ,d0.pt2() 	
        ,d0.p() 	
        ,d0.p2() 	
        // // // ,d0.pVector()
        ,d0.cpa() 	
        ,d0.cpaXY() 	
        // // // ,d0.ct() 	
        ,d0.impactParameterXY() 	
        ,d0.maxNormalisedDeltaIP()
        ,d0.px()
        ,d0.py()
        ,d0.pz()
        ,d0.eta() 	
        ,d0.phi() 	
        ,d0.y(1.86484) 	//Mass in GeV/c2
        ,d0.e(1.86484)  //Mass in GeV/c2

        ,time
        ,TFidThis
        ,bcInTF

        ,weightDSPtD0[indexPosition]
        ,weightDSPtD0Bar[indexPosition]
        ,vecTriggerMaskD0[indexPosition]

        ,hEntriesD0Pt             
        ,hEntriesD0BarPt         

        ,histPDFPol1Evaluation(recoEvent.get<TH1>(HIST("Event/D0/h_Pt" ))  , d0.pt() , 2)
        ,histPDFPol1Evaluation(recoEvent.get<TH1>(HIST("Event/D0Bar/h_Pt")), d0.pt() , 2)
      );

      d0Counter++;
    }
    PrintTime(StartD0Table, Form("DEBUG :: df_%ld :: D0 Table Time :: ",dfCount));

    auto StartPointer = std::chrono::high_resolution_clock::now();

    //Mother Listing
    std::vector<int64_t> v0GIlist; //v0GI of all possible v0 mothers of a track
    std::vector<int64_t> v0MotherList;
    std::vector<int64_t> d0GIlist;
    std::vector<int64_t> d0MotherList;
    for (const auto& track : tracks){
      if(!std::binary_search(trackGlobalIndexList.begin(), trackGlobalIndexList.end(), track.globalIndex())){continue;} //check if the track was triggered/selected or not

      //For  V0s
      v0MotherList.clear();
      v0GIlist.clear();
      // if(track.sign() > 0) { getMotherListFromDaughterList<kV0>(track.globalIndex(), posV0DauList, trigV0IndexList, v0MotherList, V0s);}
      // else                 { getMotherListFromDaughterList<kV0>(track.globalIndex(), negV0DauList, trigV0IndexList, v0MotherList, V0s);}
      // for(const auto& v0Mother: v0MotherList){
      //   for(uint i = 0; i < v0PositionList.size(); i++){if(v0Mother == v0PositionList[i][0]){ v0GIlist.push_back(v0PositionList[i][1]);}}
      // }

      //For D0s
      d0GIlist.clear();
      d0MotherList.clear();
      // // check If it is a prong of D0
      // if      (std::binary_search(prong0D0List.begin(), prong0D0List.end(), track.globalIndex())){ getMotherListFromDaughterList<kD0>(track.globalIndex(), prong0D0List, trigD0IndexList, d0MotherList, D0s);} 
      // else if (std::binary_search(prong1D0List.begin(), prong1D0List.end(), track.globalIndex())){ getMotherListFromDaughterList<kD0>(track.globalIndex(), prong1D0List, trigD0IndexList, d0MotherList, D0s);}
      // for(const auto& d0Mother: d0MotherList){
      //   for(uint i = 0; i < d0PositionList.size(); i++){if(d0Mother == d0PositionList[i][0]){ d0GIlist.push_back(d0PositionList[i][1]);}}
      // }

      // // Filling the table
      // DerivedDrTracksMothers(v0GIlist, v0MotherList, d0GIlist, d0MotherList);
    }
    PrintTime(StartPointer, Form("DEBUG :: df_%ld :: Pointer Table Time :: ",dfCount));

    // auto StartColl = std::chrono::high_resolution_clock::now();


    //Global Counters to get correct global index for merging the dataframes.
    gCollisionChecker += collisions.size();
    gFV0AChecker      += FV0As.size();
    gFT0Checker       += FT0s.size();
    gTrackChecker     += (trackQAcounter+1);
    gV0Checker        += (v0Counter);
    gPhiChecker       += (phiCounter);
    gD0Checker        += (d0Counter);
/**/
    PrintTime(Start1, Form("DEBUG :: df_%ld :: DF End    :: DF Read Time :: ",dfCount));
    PrintTime(Start0, Form("DEBUG :: df_%ld :: DF End    :: Elapsed Time :: ",dfCount));
    LOG(info)<<"DEBUG ::";

  }//Process function ends

};


WorkflowSpec defineDataProcessing(ConfigContext const& cfgc){
  return WorkflowSpec{ 
                      adaptAnalysisTask<trackqapidderivedata>(cfgc)
                      // ,adaptAnalysisTask<occupancyTableConsumer>(cfgc)
  };
}


/*

    // // // -----------------------------------Phi ------------------------------------------//

    // // for (const auto& phiCand : resoParts) {
    // //   bool isPhi = false;

    // //   float massPhi = phiCand.mass();

    // //   if (cfgMLowPhi < massPhi && massPhi < cfgMHighPhi) { isPhi = true;}
    // //   if (isPhi   ) { recoEvent.fill(HIST("Event/Phi/h_Pt"), phiCand.pt());}
    // // }
    // // PrintTime(StartPhiHist, Form("DEBUG :: df_%ld :: PhiHistFill Time :: ",dfCount));

    // // auto StartPhiDS = std::chrono::high_resolution_clock::now();

    // int nPhiTriggered = 0;
    // std::vector<int64_t> trigPhiIndexList;
    // std::vector<int64_t> prong0PhiList;
    // std::vector<int64_t> prong1PhiList;
    // std::vector<float> weightDSPtPhi;
    // std::vector<float> weightDSPtPhiBar;
    // std::vector<int> vecTriggerMaskPhi;

    // for (const auto& phiCand : Phis) {
    //   bool isPhi = false;
    //   bool isPhiBar = false;

    //   float massPhi = hfHelper.invMassPhiToPiK(phiCand);
    //   float massPhiBar = hfHelper.invMassPhibarToKPi(phiCand);

    //   if (cfgMLowPhi < massPhi && massPhi < cfgMHighPhi) { isPhi = true;}
    //   if (cfgMLowPhiBar < massPhiBar && massPhiBar < cfgMHighPhiBar) { isPhiBar = true;}

    //   if(cfgDoSkipDF && (dfCount < cfgNSkipDF)) {continue;}

    //   weightDSPt[kPhi1020] = -999.0 , probRatioPt[kPhi1020] = -999.0;
    //   weightDSPt[kPhiBar] = -999.0, probRatioPt[kPhiBar] = -999.0;

    //   isDownsampledMB[kPhi1020] = false;
    //   isDownsampledMB[kPhiBar] = false;
    //   if (isPhi   ) { isDownsampledMB[kPhi1020] = gRandom->Rndm()<factorMB;}
    //   if (isPhiBar) { isDownsampledMB[kPhiBar] = gRandom->Rndm()<factorMB;}

    //   isDownsampledPT[kPhi1020] = false;
    //   isDownsampledPT[kPhiBar] = false;
    //   if (isPhi   ) {isDownsampledPT[kPhi1020] = checkTriggerPt(phiCand.pt(), factorPt, weightDSPt[kPhi1020], probRatioPt[kPhi1020], recoEvent.get<TH1>(HIST("Event/Phi/h_Pt")));}
    //   if (isPhiBar) {isDownsampledPT[kPhiBar] = checkTriggerPt(phiCand.pt(), factorPt, weightDSPt[kPhiBar], probRatioPt[kPhiBar], recoEvent.get<TH1>(HIST("Event/PhiBar/h_Pt")));}

    //   int triggerMaskPhi = (1<<0)*isDownsampledMB[kPhi1020]+(1<<1 )*isDownsampledMB[kPhiBar]
    //                     +(1<<10)*isDownsampledPT[kPhi1020]+(1<<11)*isDownsampledPT[kPhiBar];

    //   if(cfgDoDSOfPhi && (triggerMaskPhi == 0)) { continue;}
    //   nPhiTriggered++;
    //   trigPhiIndexList.push_back(phiCand.globalIndex());
    //   prong0PhiList.push_back(phiCand.prong0Id());
    //   prong1PhiList.push_back(phiCand.prong1Id());
    //   weightDSPtPhi.push_back(weightDSPt[kPhi1020]);
    //   weightDSPtPhiBar.push_back(weightDSPt[kPhiBar]);

    //   vecTriggerMaskPhi.push_back(triggerMaskPhi);
    // }

    // recoEvent.fill(HIST("Event/percentPhiFilled"), static_cast<double>(nPhiTriggered)/static_cast<double>(Phis.size()));



    lastRun = -1; 
    for(const auto& FV0A : FV0As){
      bc = FV0A.bc_as<myBCTable>();
      run = bc.runNumber();
      GetTimingInfo(bc, lastRun, nBCsPerTF, bcSOR, time, TFidThis, bcInTF);

      // std::vector<float> vecAmplitude;// = FV0A.amplitude();
      std::vector<float> vecAmplitude(FV0A.amplitude().begin(), FV0A.amplitude().end());// = FV0A.amplitude();
      // for(auto ampVal : FV0A.amplitude()){
      //   vecAmplitude.push_back(ampVal);
      // }
      if( vecAmplitude.size() != FV0A.amplitude().size()){ LOG(info)<<"DEBUG :: Error in vecAmplitude";}

      std::vector<uint8_t> vecChannel(FV0A.channel().begin(), FV0A.channel().end());
      if( vecChannel.size() != FV0A.channel().size()) { LOG(info)<<"DEBUG :: Error in vecChannel";}

      FV0A_Debug++;
      // if(FV0A_Debug < 100) {
      //   LOG(info)<<"DEBUG :: df_"<<dfCount<<" :: AmpSize = "<<FV0A.amplitude().size()<<" :: channelSize = "<<FV0A.channel().size();
      //   for (uint i = 0 ; i < FV0A.amplitude().size(); i++){
      //     LOG(info)<<"DEBUG :: df_"<<dfCount
      //             <<" :: FV0A_Debug "<<FV0A_Debug<<" :: Amp = "<<i<<" :: "<<FV0A.amplitude()[i]<<" :: "<<vecAmplitude[i]
      //             ;
      //   }
      //   for (uint i = 0 ; i < FV0A.channel().size(); i++){
      //     LOG(info)<<"DEBUG :: df_"<<dfCount
      //             // <<" :: FV0A_Debug "<<FV0A_Debug<<" :: Chn = "<<i<<" :: "<<FV0A.channel()[i]<<" :: "<<vecChannel[i]
      //             ;
      //   }
      // } 

      DerivedFV0As(
        FV0A.globalIndex()+gFV0AChecker ,//GlobalIndex
        FV0A.globalIndex(),
        bc.globalBC()

        ,FV0A.bcId() 	      //bcId 	int32 	BC index
        ,vecAmplitude //FV0A.amplitude()    //amplitude 	std::vector<float> 	Amplitudes of non-zero channels. The channel IDs are given in Channel (at the same index)
        ,vecChannel   //FV0A.channel() 		  //channel 	std::vector<uint8_t> 	Channel IDs which had non-zero amplitudes. There are at maximum 48 channels.
        ,FV0A.time() 		    //time 	float 	Time in ns
        ,FV0A.triggerMask()  //triggerMask 	uint8_t
        // ,time
        ,TFidThis
        ,bcInTF
      );
    }
    // lastRun = -1; 

    for(const auto& FT0 : FT0s){
      bc = FT0.bc_as<myBCTable>();
      run = bc.runNumber();
      GetTimingInfo(bc, lastRun, nBCsPerTF, bcSOR, time, TFidThis, bcInTF);

      std::vector<float> vecAmplitudeA(FT0.amplitudeA().begin(), FT0.amplitudeA().end());
      if( vecAmplitudeA.size() != FT0.amplitudeA().size()){ LOG(info)<<"DEBUG :: Error in vecAmplitudeA";}

      std::vector<float> vecAmplitudeC(FT0.amplitudeC().begin(), FT0.amplitudeC().end());
      if( vecAmplitudeC.size() != FT0.amplitudeC().size()){ LOG(info)<<"DEBUG :: Error in vecAmplitudeC";}

      std::vector<uint8_t> vecChannelA(FT0.channelA().begin(), FT0.channelA().end());
      if( vecChannelA.size() != FT0.channelA().size()) { LOG(info)<<"DEBUG :: Error in vecChannelA";}

      std::vector<uint8_t> vecChannelC(FT0.channelC().begin(), FT0.channelC().end());
      if( vecChannelC.size() != FT0.channelC().size()) { LOG(info)<<"DEBUG :: Error in vecChannelC";}

      DerivedFT0s(
       FT0.globalIndex()+gFT0Checker
      ,FT0.globalIndex()
      ,bc.globalBC()

      ,FT0.bcId() 	        //I 	bcId 	int32 	BC index
      ,vecAmplitudeA 		//amplitudeA 	std::vector<float> 	Amplitudes of non-zero channels on the A-side. The channel IDs are given in ChannelA (at the same index)
      ,vecChannelA 		  //channelA 	std::vector<uint8_t> 	Channel IDs on the A side which had non-zero amplitudes. There are at maximum 96 channels.
      ,vecAmplitudeC 		//amplitudeC 	std::vector<float> 	Amplitudes of non-zero channels on the C-side. The channel IDs are given in ChannelC (at the same index)
      ,vecChannelC 		  //channelC 	std::vector<uint8_t> 	Channel IDs on the C side which had non-zero amplitudes. There are at maximum 112 channels.
      ,FT0.timeA() 		    //timeA 	float 	Average A-side time
      ,FT0.timeC() 		    //timeC 	float 	Average C-side time
      ,FT0.triggerMask() 	//triggerMask 	uint8_t 	
      ,FT0.posZ() 	        //posZ     	float 	Z position calculated from timeA and timeC in cm
      ,FT0.collTime() 	    //collTime 	float 	Collision time, one need also check validation (code below) for timeA and timeC
      ,FT0.isValidTimeA() 	//isValidTimeA 	bool 	Checks if time from A side was calculated, and if is not dummy
      ,FT0.isValidTimeC() 	//isValidTimeC 	bool 	Checks if time from C side was calculated
      ,FT0.isValidTime() 	//isValidTime 	bool 	Checks if times from A and C side were calculated simultaneously
      ,FT0.sumAmpA() 	    //sumAmpA 	float 	Calculates sum of positive amplitudes from side A
      ,FT0.sumAmpC() 	    //sumAmpC 	float 	Calculates sum of positive amplitudes from side C
      ,time
      ,TFidThis
      ,bcInTF
      );
    }
  */