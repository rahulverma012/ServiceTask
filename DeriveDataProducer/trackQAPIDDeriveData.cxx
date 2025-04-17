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

#include "DerivedPIDTables.h"
#include "Common/DataModel/OccupancyTables.h"

#include <typeinfo>
#include <TRandom.h>

#include "Framework/AnalysisDataModel.h"
#include "Framework/HistogramRegistry.h"
#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsFT0/Digit.h"
#include "DataFormatsParameters/GRPLHCIFData.h"
// #include "MetadataHelper.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

int32_t nBCsPerOrbit = o2::constants::lhc::LHCMaxBunches;

int hPt_BinFor1        = -1; //= recoEvent.get<TH1>(HIST("Event/track/h_Pt"  )->FindBin(1))
int hQPt_BinFor_Plus1  = -1; //= recoEvent.get<TH1>(HIST("Event/track/h_QPt" )->FindBin(-1))
int hQPt_BinFor_Minus1 = -1; //= recoEvent.get<TH1>(HIST("Event/track/h_QPt" )->FindBin(-1))

struct trackqapidderivedata{

  Service<o2::ccdb::BasicCCDBManager> ccdb;

  // declare production of tables
  Produces<aod::DrCollisions> DerivedCollisions;
  Produces<aod::DrTracks> DerivedTracks;
  Produces<aod::DrFV0As> DerivedFV0As;
  Produces<aod::DrFT0s> DerivedFT0s;
  Produces<aod::DrV0s> DerivedV0s;
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
  Configurable<int> cfgSkipDataFrame{"cfgSkipDataFrame", 3, "cfgSkipDataFrame"};
  //check Fill condition for downsampling

  Configurable<bool> doV0andTrackMatchingCheck{"doV0andTrackMatchingCheck", true, "doV0andTrackMatchingCheck"};

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

  void init(InitContext const&){
    LOGF(info,"Starting init");

    // auto grplhcif = ccdb->getForTimeStamp<o2::parameters::GRPLHCIFData>("GLO/Config/GRPLHCIF", tsSOR);
    // bcPatternB = grplhcif->getBunchFilling().getBCPattern();
    // auto vector = ccdb->getForTimeStamp<std::vector<Long64_t>>(1731764207045, 1731850607045);

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
      recoEvent.add("Event/track/Full/h_QPt", "h_QPt", {HistType::kTH1F, {{200,-10,10}}});

      recoEvent.addClone("Event/track/Full/","Event/track/Rejected/");
      recoEvent.addClone("Event/track/Rejected/","Event/track/Accepted/");

      recoEvent.add("Event/track/isDownsampledMB" , "isDownsampledMB" , kTH1F, {{33*2,-1, 32}});  
      recoEvent.add("Event/track/isDownsampledPT" , "isDownsampledPT" , kTH1F, {{33*2,-1, 32}});  
      recoEvent.add("Event/track/isDownsampledQPT", "isDownsampledQPT", kTH1F, {{33*2,-1, 32}}); 
      recoEvent.add("Event/track/FullMask"        , "FullMask"        , kTH1F, {{33*2,-1, 32}});

      recoEvent.add("Event/track/probRatioPT" , "probRatioPT" , kTH1F, {{400,-20, 20}});  
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

      recoEvent.add("Event/percentV0Filled", "percentV0Filled", kTH1F, {{1020, -1, 102}});

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
    kAntiHypertriton
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

  int FindTrackIdInList(int64_t Key, std::vector<int64_t> trackGlobalIndexList, std::vector<int64_t> DrTrackPosList){
    int elementPos = -1;
    if (trackGlobalIndexList.size() > 0) { 
      elementPos = BinarySearchVector(Key,trackGlobalIndexList, 0, trackGlobalIndexList.size());
    } 
    if (elementPos == -1){ return -2000000000;}
    else                 {
      if ( DrTrackPosList[elementPos] < 0 ) { LOG(info) <<"DEBUG :: ERROR :: ERROR in DrTrackPosList[elementPos]";}
      return DrTrackPosList[elementPos];
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
  int BinarySearchVector(int64_t Key, std::vector<int64_t> List, int low, int high)
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
    float prob  = static_cast<float>(static_cast<double>(hist->GetBinContent(hist->FindBin(Qpt)))                                                            /static_cast<double>(hist->GetEntries()));
    float prob1 = static_cast<float>(static_cast<double>(hist->GetBinContent(hQPt_BinFor_Plus1))                                                             /static_cast<double>(hist->GetEntries()));
    // float prob1 = static_cast<float>(static_cast<double>(hist->GetBinContent(hQPt_BinFor_Plus1))+static_cast<double>(hist->GetBinContent(hQPt_BinFor_Minus1))/static_cast<double>(hist->GetEntries()));
    weight = prob;
    probRatio = prob/prob1;
    bool isTrigger = (gRandom->Rndm()*(prob/prob1))<factorDS;
    return isTrigger;
  }

  template<typename T>
  void getMotherListFromDaughterList( const int64_t &Key, const std::vector<int64_t> &dauList, std::vector<int64_t> &trigV0List, std::vector<int64_t> &newMotherList, const T& V0s){
    for(uint32_t i = 0; i < dauList.size(); i++){
      if(dauList[i] == Key) {
        // LOG(info)<<"DEBUG :: Key :: "<<Key<<" :: pos = "<<i<<" :: trigV0 = "<<trigV0List[i]<<" :: GI of V0 at it = "<<V0s.iteratorAt(trigV0List[i]).globalIndex()<<" :: posTrackId() = "<<V0s.iteratorAt(trigV0List[i]).posTrackId() ;
        if(doV0andTrackMatchingCheck){ if(V0s.iteratorAt(trigV0List[i]).posTrackId() != Key && V0s.iteratorAt(trigV0List[i]).negTrackId() != Key){
          LOG(info)<<"DEBUG :: v0 Information checking :: something is wrong :: check the errors :: Indices not mathcing";
          LOG(info)<<"Key :: "<<Key<<" :: ";
        }}
        newMotherList.push_back(trigV0List[i]);
      }
    }
  }


  // Event Filter
  // Filter eventFilter = (o2::aod::evsel::sel8 == true);
  // Filter posZFilter  = (nabs(o2::aod::collision::posZ) < cutZvertex);

  // Track Filter
  // Filter PtFilter  = (o2::aod::track::pt) > 0.15f && (o2::aod::track::pt) < 5.0f ;

  using myCollisions = soa::Join<aod::Collisions, aod::Mults, aod::MultsExtra, aod::FT0sCorrected, aod::EvSels>;
  
  //, aod::EvSels>>;//;, 
                      // aod::CentFT0Ms, aod::CentFT0Cs, aod::CentFT0As, aod::Mults>>; //aod::CentFV0As,

                      //aod::FV0As, o2::aod::FT0s These tables cant be joined to the collision tables
                      
                      // o2::aod::Mults = soa::Join<o2::aod::BarrelMults, o2::aod::FV0Mults, o2::aod::FT0Mults, o2::aod::FDDMults, o2::aod::ZDCMults> 

  using myTracks = soa::Join<aod::Tracks, o2::aod::TracksCov, aod::TracksExtra, aod::TracksDCA, aod::TracksDCACov //>; //, aod::TracksQA>;//, , aod::TrackSelection
                    ,aod::TOFSignal,  aod::pidTOFbeta, aod::pidTOFmass, aod::EvTimeTOFOnly, aod::pidTOFFlags, aod::pidEvTimeFlags
                    ,aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTPCFullEl, aod::pidTPCFullDe
                    ,aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFFullEl, aod::pidTOFFullDe>;

  using myOccsTable = soa::Join<aod::OccsBCsList, aod::OccsDet, aod::OccsTrackMult, aod::OccsMultExtra, aod::OccsRobust,
                                                  aod::OccsMeanDet, aod::OccsMeanTrkMult, aod::OccsMnMultExtra, aod::OccsMeanRobust>;
  
  using myTrackMeanOccs = soa::Join<aod::TrackMeanOccs0, aod::TrackMeanOccs1, aod::TrackMeanOccs2, aod::TrackMeanOccs3, aod::TrackMeanOccs4,
                                                         aod::TrackMeanOccs5, aod::TrackMeanOccs6, aod::TrackMeanOccs7, aod::TrackMeanOccs8>;

  using myBCTable = soa::Join<aod::BCsWithTimestamps,aod::OccIndexTable>;
  // using MyTracksQA = aod::TracksQA_000;// aod::TracksQAVersion; //aod::TracksQA
  // using MyTracksQA = aod::TracksQA_002;  //aod::TracksQA_000;// aod::TracksQAVersion; //aod::TracksQA
  using MyTracksQA = aod::TracksQA_002;


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
  int dfCount = 0;
  int debugCounter = 0;
  int CountCheck = 0;
  int FV0A_Debug = 0;

  uint64_t gCollisionChecker = 0;
  uint64_t gFV0AChecker = 0;
  uint64_t gFT0Checker  = 0;
  uint64_t gTrackChecker = 0;
  uint64_t gV0Checker    = 0;
  
  // int      dfCount = 0;
  int lastRun = -999;
  int64_t bcSOR = -1;
  int32_t nBCsPerTF = -1;
  uint64_t time = -1;
  int64_t TFidThis = -1;
  int bcInTF = -1;
  
  bool startDownsampling = false;

  void process( myBCTable const& BCs
               ,myCollisions const& collisions
               ,myTracks const& tracks
               ,MyTracksQA const& tracksQA
               ,o2::aod::Origins const& Origins
               ,o2::aod::AmbiguousTracks const& ambgTracks
               ,o2::aod::FV0As const& FV0As
               //,o2::aod::FV0Cs const& FV0Cs
               ,o2::aod::FT0s  const& FT0s
               ,o2::aod::V0Datas  const& V0s
               ,o2::aod::OccIndexTable const& occIdxTable
               ,myOccsTable const& occTables
               ,myTrackMeanOccs const& trackMeanOccs
              )
  {
    lastRun = -999;
    nBCsPerTF = -999;
    bcSOR = -999;
    dfCount++;
    // if(dfCount > 10) {return;}
    LOG(info)<<"DEBUG :: df_"<<dfCount<<" :: DF_"<<Origins.iteratorAt(0).dataframeID()<<" :: collisions.size() = "<<collisions.size()<<" :: tracks.size() = "<<tracks.size()<<" :: tracksQA.size() = "<<tracksQA.size();
    LOG(info)<<"DEBUG :: df_"<<dfCount<<" :: DF_"<<Origins.iteratorAt(0).dataframeID()
                                                                                      <<" :: FV0As.size() = "<<FV0As.size()
                                                                                      // <<" :: FV0Cs.size() = "<<FV0Cs.size()
                                                                                      <<" :: FT0s.size() = "<<FT0s.size()
                                                                                      <<" :: V0s.size() = "<<V0s.size()
                                                                                      <<" :: BCs.size() = "<<BCs.size()
                                                                                      <<" :: occIdxTable.size() = "<<occIdxTable.size()
                                                                                      <<" :: occTables.size() = "<<occTables.size()
                                                                                      <<" :: tracKMeanOccs.size() = "<<trackMeanOccs.size()
                                                                                      ;

    if( dfCount < 10 ) { 
      LOG(info)<<"DEBUG :: df_"<<dfCount<<" :: DF_"<<Origins.iteratorAt(0).dataframeID()
                                        <<" :: FT0s :: posZ = "        <<FT0s.iteratorAt(0).posZ()
                                                <<" :: CollTime = "    <<FT0s.iteratorAt(0).collTime()
                                                <<" :: IsValidTimeA = "<<FT0s.iteratorAt(0).isValidTimeA()
                                                <<" :: IsValidTimeC = "<<FT0s.iteratorAt(0).isValidTimeC()
                                                <<" :: IsValidTime = " <<FT0s.iteratorAt(0).isValidTime()
                                                <<" :: SumAmpA = "     <<FT0s.iteratorAt(0).sumAmpA()
                                                <<" :: SumAmpC = "     <<FT0s.iteratorAt(0).sumAmpC()
                                                ;
    }

    for(auto const& myocc : occTables){
      std::vector<int64_t> myVector = std::vector<int64_t>(myocc.bcsInTFList().begin(),myocc.bcsInTFList().end());
      BuildOCCS(
                  myocc.tfId(),
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
      ,collision.bcInTF()  //	int	Position of a (found) bunch crossing inside a given timeframe
      ,collision.trackOccupancyInTimeRange()  //	int	Occupancy in specified time interval by a number of tracks from nearby collisions
      ,collision.ft0cOccupancyInTimeRange()  //	float	Occupancy in specified time interval by a sum of FT0C amplitudes from nearby collisions

      ,collision.t0ACorrected()	    //float	Collision time A-side, corrected with primary vertex
      ,collision.t0CCorrected()	    //float	Collision time C-side, corrected with primary vertex
      ,collision.t0AC()	            //float	Collision time (A+C)/2
      ,collision.t0ACorrectedValid()	//bool	Was T0ACorrected computable?
      ,collision.t0CCorrectedValid()	//bool	Was T0CCorrected computable?
      ,collision.t0ACValid()	        //bool	Was T0AC computable?
      ,collision.t0resolution()	    //float	Was T0CCorrected computable?
      );
    }
    //collision Loop is over

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

    std::vector<int> vecTriggerMaskV0MB;
    std::vector<int> vecTriggerMaskV0PT;

    int nV0Triggered = 0;

    float weightDSPt[10] ; float probRatioPT[10];
    bool isDownsampledMB[10]; bool isDownsampledPT[10];

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

      if(dfCount < cfgSkipDataFrame) continue;
      weightDSPt[kV0]              = -9999.0 , probRatioPT[kV0]              = -9999.0;
      weightDSPt[kLambda]          = -9999.0 , probRatioPT[kLambda]          = -9999.0;
      weightDSPt[kAntiLambda]      = -9999.0 , probRatioPT[kAntiLambda]      = -9999.0;
      weightDSPt[kK0Short]         = -9999.0 , probRatioPT[kK0Short]         = -9999.0;
      weightDSPt[kGamma]           = -9999.0 , probRatioPT[kGamma]           = -9999.0;
      weightDSPt[kHypertriton]     = -9999.0 , probRatioPT[kHypertriton]     = -9999.0;
      weightDSPt[kAntiHypertriton] = -9999.0 , probRatioPT[kAntiHypertriton] = -9999.0;

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

      isDownsampledPT[kV0]              = checkTriggerPt(v0.pt(), factorPt, weightDSPt[kV0]              , probRatioPT[kV0]              , recoEvent.get<TH1>(HIST("Event/V0/h_Pt"              )));
      if( isLambda          ){isDownsampledPT[kLambda]          = checkTriggerPt(v0.pt(), factorPt, weightDSPt[kLambda]          , probRatioPT[kLambda]          , recoEvent.get<TH1>(HIST("Event/Lambda/h_Pt"          )));}
      if( isAntiLambda      ){isDownsampledPT[kAntiLambda]      = checkTriggerPt(v0.pt(), factorPt, weightDSPt[kAntiLambda]      , probRatioPT[kAntiLambda]      , recoEvent.get<TH1>(HIST("Event/AntiLambda/h_Pt"      )));}
      if( isK0Short         ){isDownsampledPT[kK0Short]         = checkTriggerPt(v0.pt(), factorPt, weightDSPt[kK0Short]         , probRatioPT[kK0Short]         , recoEvent.get<TH1>(HIST("Event/K0Short/h_Pt"         )));}
      // if( isGamma           ){isDownsampledPT[kGamma]           = checkTriggerPt(v0.pt(), factorPt, weightDSPt[kGamma]           , probRatioPT[kGamma]           , recoEvent.get<TH1>(HIST("Event/Gamma/h_Pt"           )));}
      // if( isHypertriton     ){isDownsampledPT[kHypertriton]     = checkTriggerPt(v0.pt(), factorPt, weightDSPt[kHypertriton]     , probRatioPT[kHypertriton]     , recoEvent.get<TH1>(HIST("Event/Hypertriton/h_Pt"     )));}
      // if( isAntiHypertriton ){isDownsampledPT[kAntiHypertriton] = checkTriggerPt(v0.pt(), factorPt, weightDSPt[kAntiHypertriton] , probRatioPT[kAntiHypertriton] , recoEvent.get<TH1>(HIST("Event/AntiHypertriton/h_Pt" )));}
      
      int triggerMaskV0MB = (1<<0)*isDownsampledMB[kV0]+(1<<1)*isDownsampledMB[kLambda]+(1<<2)*isDownsampledMB[kAntiLambda]+(1<<3)*isDownsampledMB[kK0Short];
      int triggerMaskV0PT = (1<<0)*isDownsampledPT[kV0]+(1<<1)*isDownsampledPT[kLambda]+(1<<2)*isDownsampledPT[kAntiLambda]+(1<<3)*isDownsampledPT[kK0Short];

      if(triggerMaskV0MB == 0 && triggerMaskV0PT == 0) { continue;}
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

      vecTriggerMaskV0MB.push_back(triggerMaskV0MB);
      vecTriggerMaskV0PT.push_back(triggerMaskV0PT);
    }

    recoEvent.fill(HIST("Event/percentV0Filled"), static_cast<double>(nV0Triggered)/static_cast<double>(V0s.size()));

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
    std::vector<int64_t> DrTrackPosList;
    std::vector<int64_t> v0MotherList;
    int trackQAcounter = -1;
    int trackSkippedBecauseOfMask = 0;

    // using myTracks = soa::Join<aod::Tracks, o2::aod::TracksCov, aod::TracksExtra //>; //, aod::TracksQA>;//, aod::TracksDCA, aod::TrackSelection
    //                   ,aod::TOFSignal,  aod::pidTOFbeta, aod::pidTOFmass, aod::EvTimeTOFOnly
    //                   ,aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTPCFullEl, aod::pidTPCFullDe
    //                   ,aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFFullEl, aod::pidTOFFullDe>;


    // auto tracksWithPid = soa::Attach<o2::soa::Join<aod::TracksIU, aod::TracksExtra>,
    //                                   aod::pidits::ITSNSigmaEl, aod::pidits::ITSNSigmaMu, aod::pidits::ITSNSigmaPi,
    //                                   aod::pidits::ITSNSigmaKa, aod::pidits::ITSNSigmaPr, aod::pidits::ITSNSigmaDe,
    //                                   aod::pidits::ITSNSigmaTr, aod::pidits::ITSNSigmaHe, aod::pidits::ITSNSigmaAl>(tracks);

    // using myTracks = soa::Join<aod::Tracks, o2::aod::TracksCov, aod::TracksExtra //>; //, aod::TracksQA>;//, aod::TracksDCA, aod::TrackSelection
    //                   ,aod::TOFSignal,  aod::pidTOFbeta, aod::pidTOFmass, aod::EvTimeTOFOnly
    //                   ,aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTPCFullEl, aod::pidTPCFullDe
    //                   ,aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFFullEl, aod::pidTOFFullDe>;


    auto tracksWithPid = soa::Attach<myTracks,
                              aod::pidits::ITSNSigmaEl, aod::pidits::ITSNSigmaMu, aod::pidits::ITSNSigmaPi,
                              aod::pidits::ITSNSigmaKa, aod::pidits::ITSNSigmaPr, aod::pidits::ITSNSigmaDe,
                              aod::pidits::ITSNSigmaTr, aod::pidits::ITSNSigmaHe, aod::pidits::ITSNSigmaAl>(tracks);

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

      recoEvent.fill(HIST("Event/track/Full/h_Pt" ),track.pt());
      recoEvent.fill(HIST("Event/track/Full/h_QPt"),track.signed1Pt());

      float weightDSPt = -1.0;
      float weightDSQPt = -1.0;
      int  triggerMaskDS = 0;

      bool isDownsampledMB  = false;
      bool isDownsampledPT  = false;
      bool isDownsampledQPT = false;
      bool isDownsampledV0  = false;

      float probRatioPT = -9999999;
      float probRatioQPT = -9999999;
      
      if(dfCount < cfgSkipDataFrame) continue; // what if data frames are empty // In data you did encountered empty df's

      // //check non zero count 
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
      isDownsampledPT  =  checkTriggerPt(track.pt()       , factorPt, weightDSPt , probRatioPT , recoEvent.get<TH1>(HIST("Event/track/Full/h_Pt" )));
      isDownsampledQPT = checkTriggerQPt(track.signed1Pt(), factorPt, weightDSQPt, probRatioQPT, recoEvent.get<TH1>(HIST("Event/track/Full/h_QPt")));

      // check If It is a V0 daughter and getItsV0s. 
      if(track.sign() > 0) { isDownsampledV0 = std::binary_search(posV0DauList.begin(), posV0DauList.end(), track.globalIndex());}
      else                 { isDownsampledV0 = std::binary_search(negV0DauList.begin(), negV0DauList.end(), track.globalIndex());}

      triggerMaskDS=(1<<0)*isDownsampledMB+(1<<1)*isDownsampledPT+(1<<2)*isDownsampledQPT+(1<<3)*isDownsampledV0;
      
      recoEvent.fill(HIST("Event/track/isDownsampledMB" ), isDownsampledMB );  
      recoEvent.fill(HIST("Event/track/isDownsampledPT" ), isDownsampledPT );  
      recoEvent.fill(HIST("Event/track/isDownsampledQPT"), isDownsampledQPT); 
      recoEvent.fill(HIST("Event/track/FullMask"        ), triggerMaskDS);

      recoEvent.fill(HIST("Event/track/probRatioPT" ), probRatioPT );
      recoEvent.fill(HIST("Event/track/probRatioQPT"), probRatioQPT);

      // TH1F *hist1 = (TH1F*)_file0->Get("trackqapidderivedata/recoEvent/Event/track/isDownsampledMB");
      // TH1F *hist2 = (TH1F*)_file0->Get("trackqapidderivedata/recoEvent/Event/track/isDownsampledPT");
      // TH1F *hist3 = (TH1F*)_file0->Get("trackqapidderivedata/recoEvent/Event/track/isDownsampledQPT");
      // TH1F *hist4 = (TH1F*)_file0->Get("trackqapidderivedata/recoEvent/Event/track/FullMask");

      // // // auto hist=h1;
      // template<typename T> void printNonZeroBins (const T& hist){ 
      //   cout<<"Hist Entries   : "<<hist->GetEntries()<<endl;
      //   cout<<"Hist Overflow  : "<<hist->GetBinContent(hist->GetNbinsX()+1)<<endl;
      //   cout<<"Hist Underflow : "<<hist->GetBinContent(0)<<endl;
      //   for(int i = 0; i < hist->GetNbinsX(); i++){
      //     if(hist->GetBinContent(i) != 0) {
      //       cout<<"bin i = "<<i
      //           <<" :: content = "<<std::setprecision(20)<<hist->GetBinContent(i)
      //           <<" :: ratTot  = "<<std::setprecision(7)<<100*static_cast<double>(hist->GetBinContent(i))/static_cast<double>(hist->GetEntries())
      //           <<" :: ratInt  = "<<std::setprecision(7)<<100*static_cast<double>(hist->GetBinContent(i))/static_cast<double>(hist->Integral(1,hist->GetNbinsX()))
      //           <<endl;}
      //   }
      //   cout<<endl;
      // }

      // printNonZeroBins(hist1);
      // printNonZeroBins(hist2);
      // printNonZeroBins(hist3);
      // printNonZeroBins(hist4);

      // track-q-a-p-i-d-derive-data
      // 3046153
      // 338399/3046153

      if( triggerMaskDS == 0) {
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
      DrTrackPosList.push_back(trackQAcounter);
      
      auto const& trackMeanOcc = trackMeanOccs.iteratorAt(trackQAcounter+trackSkippedBecauseOfMask); //(trackQA.globalIndex()); //because of ambgTracks, you cant use globalIndex()
      if( trackMeanOcc.trackId() != trackQA.trackId()   ){ LOG(info)<<"DEBUG :: ERROR :: TrackQA and trackMeanOcc mismatch";}

      //Store Track Data in trees
      DerivedTracks(
         trackQAcounter+gTrackChecker     //GlobalIndex
        ,track.globalIndex()                   //OrignalIndex
        ,bc.globalBC()                         //GlobalBC
        ,track.collisionId()+gCollisionChecker //GCollId

        ,track.collisionId()                   //DrCollisionId
        ,track.collisionId() 	      // Collision to which this track belongs
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
        ,track.hasITS() 	                          // Flag to check if track has a ITS match
        ,track.hasTPC() 	                          // Flag to check if track has a TPC match
        ,track.hasTRD() 	                          // Flag to check if track has a TRD match
        ,track.hasTOF() 	                          // Flag to check if track has a TOF measurement
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

        ,trackMeanOcc.meanOccPrimUnfm80() 
        ,trackMeanOcc.meanOccFV0AUnfm80() 
        ,trackMeanOcc.meanOccFV0CUnfm80() 
        ,trackMeanOcc.meanOccFT0AUnfm80() 
        ,trackMeanOcc.meanOccFT0CUnfm80() 
        ,trackMeanOcc.meanOccFDDAUnfm80() 
        ,trackMeanOcc.meanOccFDDCUnfm80() 
        ,trackMeanOcc.meanOccNTrackITSUnfm80()   
        ,trackMeanOcc.meanOccNTrackTPCUnfm80()   
        ,trackMeanOcc.meanOccNTrackTRDUnfm80()   
        ,trackMeanOcc.meanOccNTrackTOFUnfm80()   
        ,trackMeanOcc.meanOccNTrackSizeUnfm80()   
        ,trackMeanOcc.meanOccNTrackTPCAUnfm80()  
        ,trackMeanOcc.meanOccNTrackTPCCUnfm80()  
        ,trackMeanOcc.meanOccNTrackITSTPCUnfm80()
        ,trackMeanOcc.meanOccNTrackITSTPCAUnfm80()
        ,trackMeanOcc.meanOccNTrackITSTPCCUnfm80()
        ,trackMeanOcc.meanOccMultNTracksHasITSUnfm80()
        ,trackMeanOcc.meanOccMultNTracksHasTPCUnfm80()
        ,trackMeanOcc.meanOccMultNTracksHasTOFUnfm80()
        ,trackMeanOcc.meanOccMultNTracksHasTRDUnfm80()
        ,trackMeanOcc.meanOccMultNTracksITSOnlyUnfm80()
        ,trackMeanOcc.meanOccMultNTracksTPCOnlyUnfm80()
        ,trackMeanOcc.meanOccMultNTracksITSTPCUnfm80()
        ,trackMeanOcc.meanOccMultAllTracksTPCOnlyUnfm80()
        ,trackMeanOcc.meanOccRobustT0V0PrimUnfm80()   
        ,trackMeanOcc.meanOccRobustFDDT0V0PrimUnfm80()
        ,trackMeanOcc.meanOccRobustNtrackDetUnfm80()  
        ,trackMeanOcc.meanOccRobustMultExtraTableUnfm80()  
        ,trackMeanOcc.weightMeanOccPrimUnfm80() 
        ,trackMeanOcc.weightMeanOccFV0AUnfm80() 
        ,trackMeanOcc.weightMeanOccFV0CUnfm80() 
        ,trackMeanOcc.weightMeanOccFT0AUnfm80() 
        ,trackMeanOcc.weightMeanOccFT0CUnfm80() 
        ,trackMeanOcc.weightMeanOccFDDAUnfm80() 
        ,trackMeanOcc.weightMeanOccFDDCUnfm80() 
        ,trackMeanOcc.weightMeanOccNTrackITSUnfm80()   
        ,trackMeanOcc.weightMeanOccNTrackTPCUnfm80()   
        ,trackMeanOcc.weightMeanOccNTrackTRDUnfm80()   
        ,trackMeanOcc.weightMeanOccNTrackTOFUnfm80()   
        ,trackMeanOcc.weightMeanOccNTrackSizeUnfm80()   
        ,trackMeanOcc.weightMeanOccNTrackTPCAUnfm80()  
        ,trackMeanOcc.weightMeanOccNTrackTPCCUnfm80()  
        ,trackMeanOcc.weightMeanOccNTrackITSTPCUnfm80()
        ,trackMeanOcc.weightMeanOccNTrackITSTPCAUnfm80()
        ,trackMeanOcc.weightMeanOccNTrackITSTPCCUnfm80()
        ,trackMeanOcc.weightMeanOccMultNTracksHasITSUnfm80()
        ,trackMeanOcc.weightMeanOccMultNTracksHasTPCUnfm80()
        ,trackMeanOcc.weightMeanOccMultNTracksHasTOFUnfm80()
        ,trackMeanOcc.weightMeanOccMultNTracksHasTRDUnfm80()
        ,trackMeanOcc.weightMeanOccMultNTracksITSOnlyUnfm80()
        ,trackMeanOcc.weightMeanOccMultNTracksTPCOnlyUnfm80()
        ,trackMeanOcc.weightMeanOccMultNTracksITSTPCUnfm80()
        ,trackMeanOcc.weightMeanOccMultAllTracksTPCOnlyUnfm80()
        ,trackMeanOcc.weightMeanOccRobustT0V0PrimUnfm80()   
        ,trackMeanOcc.weightMeanOccRobustFDDT0V0PrimUnfm80()
        ,trackMeanOcc.weightMeanOccRobustNtrackDetUnfm80()  
        ,trackMeanOcc.weightMeanOccRobustMultExtraTableUnfm80()
        
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
        ,weightDSPt
        ,weightDSQPt
        ,triggerMaskDS
      );
      trackCounter++;
    }// end of trackQA loop

    int lastCollId = -10;
    int v0Counter = 0;

    int PosIndexPosition = -1;
    int NegIndexPosition = -1;

    std::vector<std::vector<int64_t>> v0PositionList;
    
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
      if(!std::binary_search(trigV0IndexList.begin(), trigV0IndexList.end(), v0.globalIndex())){continue;}
      PosIndexPosition = FindTrackIdInList(v0.posTrackId(),trackGlobalIndexList, DrTrackPosList);
      NegIndexPosition = FindTrackIdInList(v0.negTrackId(),trackGlobalIndexList, DrTrackPosList);
      if(PosIndexPosition < 0 && NegIndexPosition < 0) {continue;}      

      v0PositionList.push_back({v0.globalIndex(), static_cast<int64_t>(v0Counter+gV0Checker)});
      int indexPosition = BinarySearchVector(v0.globalIndex(), trigV0IndexList, 0, trigV0IndexList.size()-1);
      
      //weightDSPtV0[indexPosition];
      //weightDSPtLambda[indexPosition];
      //weightDSPtAntiLambda[indexPosition];
      //weightDSPtK0Short[indexPosition];
      //weightDSPtGamma[indexPosition];
      //weightDSPtHypertriton[indexPosition];
      //weightDSPtAntiHypertriton[indexPosition];


      DerivedV0s(
        v0Counter+gV0Checker
        ,v0.globalIndex()
        ,bc.globalBC()
        ,v0.collisionId()+gCollisionChecker
        ,v0Counter 
        ,PosIndexPosition+gTrackChecker
        ,NegIndexPosition+gTrackChecker
        ,PosIndexPosition
        ,NegIndexPosition
        ,v0.posTrackId() 	  //int 	Pointer into Tracks
        ,v0.negTrackId() 	  //int 	Pointer into Tracks
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
        
        ,weightDSPtV0[indexPosition]
        ,weightDSPtLambda[indexPosition]
        ,weightDSPtAntiLambda[indexPosition]
        ,weightDSPtK0Short[indexPosition]
        ,weightDSPtGamma[indexPosition]
        ,weightDSPtHypertriton[indexPosition]
        ,weightDSPtAntiHypertriton[indexPosition]
        ,vecTriggerMaskV0MB[indexPosition]
        ,vecTriggerMaskV0PT[indexPosition]
      );
      v0Counter++;
    }

    std::vector<int64_t> v0GIlist;
    for (const auto& track : tracks){
      v0MotherList.clear();
      v0GIlist.clear();
      if(!std::binary_search(trackGlobalIndexList.begin(), trackGlobalIndexList.end(), track.globalIndex())){continue;}
      if(track.sign() > 0) { getMotherListFromDaughterList(track.globalIndex(), posV0DauList, trigV0IndexList, v0MotherList, V0s);}
      else                 { getMotherListFromDaughterList(track.globalIndex(), negV0DauList, trigV0IndexList, v0MotherList, V0s);}
      
      for(const auto& v0Mother: v0MotherList){
        for(uint i = 0; i < v0PositionList.size(); i++){if(v0Mother == v0PositionList[i][0]){ v0GIlist.push_back(v0PositionList[i][1]);}}
      }

      DerivedDrTracksMothers(v0GIlist, v0MotherList);
    }


  //   trackGlobalIndexList
  //   //Making a pointer Table ;
  //   for (const auto& track : tracks()){
  //     if(track.sign() > 0 ) 
  //   }
  //   trackGlobalIndexList.push_back(track.globalIndex());

  //   DerivedTracks(
  //     trackQAcounter+gTrackChecker     //GlobalIndex
  //    ,track.globalIndex()                   //OrignalIndex
  //    ,bc.globalBC()                         //GlobalBC
  //    ,track.collisionId()+gCollisionChecker //GCollId

  //    ,track.collisionId()                   //DrCollisionId
  //    ,track.collisionId() 	      // Collision to which this track belongs
  //    ,track.trackType() 	        // Type of track. See enum TrackTypeEnum. This cannot be used to decide which detector has contributed to this track. Use hasITS, hasTPC, etc.
  //    ,track.x() 	                // 
  //    ,track.alpha() 	            // 
  //    ,track.y() 	                // 
  //    ,track.z() 	                // 
  //    ,track.snp() 	              // 
  //    ,track.tgl() 	              // 
  //    ,track.signed1Pt() 	        // (sign of charge)/Pt in c/GeV. Use pt() and sign() instead
  //    ,track.pt() 	              // Transverse momentum of the track in GeV/c
  //    ,track.p() 	                // Momentum in Gev/c
  //    ,track.eta() 	              // Pseudorapidity
  //    ,track.phi() 	              // Phi of the track, in radians within [0, 2pi)

  //    ,track.isWithinBeamPipe() 	  // Is the track within the beam pipe (= successfully propagated to a collision vertex)
  //    ,track.px() 	                // Momentum in x-direction in GeV/c
  //    ,track.py() 	                // Momentum in y-direction in GeV/c
  //    ,track.pz() 	                // Momentum in z-direction in GeV/c
  //    // ,drTrack::PVector         // D 	pVector             std::array<float,3> 	Momentum vector in x,y,z-directions in GeV/c
  //    ,track.energy(0.13957) 	            // Track energy, computed under the mass assumption given as input
  //    ,track.rapidity(0.13957) 	          // Track rapidity, computed under the mass assumption given as input
  //    ,track.sign() 	              // Charge: positive: 1, negative: -1

  //    // // TracksCov
  //    ,track.sigmaY() 	    //Covariance matrix
  //    ,track.sigmaZ() 	    //Covariance matrix
  //    ,track.sigmaSnp() 	  //Covariance matrix
  //    ,track.sigmaTgl() 	  //Covariance matrix
  //    ,track.sigma1Pt() 	  //Covariance matrix
  //    ,track.rhoZY() 	    //Covariance matrix in compressed form
  //    ,track.rhoSnpY() 	 	//Covariance matrix in compressed form
  //    ,track.rhoSnpZ() 	 	//Covariance matrix in compressed form
  //    ,track.rhoTglY() 	 	//Covariance matrix in compressed form
  //    ,track.rhoTglZ() 	 	//Covariance matrix in compressed form
  //    ,track.rhoTglSnp() 	//Covariance matrix in compressed form
  //    ,track.rho1PtY() 	 	//Covariance matrix in compressed form
  //    ,track.rho1PtZ() 	 	//Covariance matrix in compressed form
  //    ,track.rho1PtSnp() 	//Covariance matrix in compressed form
  //    ,track.rho1PtTgl() 	//Covariance matrix in compressed form
  //    ,track.cYY() 	  //
  //    ,track.cZY() 	  
  //    ,track.cZZ() 	  
  //    ,track.cSnpY() 	  
  //    ,track.cSnpZ() 	  
  //    ,track.cSnpSnp()
  //    ,track.cTglY() 	  
  //    ,track.cTglZ() 	  
  //    ,track.cTglSnp()
  //    ,track.cTglTgl()
  //    ,track.c1PtY() 	  
  //    ,track.c1PtZ() 	  
  //    ,track.c1PtSnp()
  //    ,track.c1PtTgl()
  //    ,track.c1Pt21Pt2()

  //    // // TracksExtra
  //    ,track.tpcInnerParam() 	                    // Momentum at inner wall of the TPC
  //    ,track.flags() 	                            // Track flags. Run 2: see TrackFlagsRun2Enum | Run 3: see TrackFlags
  //    ,track.itsClusterSizes() 	                  // Clusters sizes, four bits per a layer, starting from the innermost
  //    ,track.tpcNClsFindable() 	                  // Findable TPC clusters for this track geometry
  //    ,track.tpcNClsFindableMinusFound() 	        // TPC Clusters: Findable - Found
  //    ,track.tpcNClsFindableMinusCrossedRows() 	  // TPC Clusters: Findable - crossed rows
  //    ,track.tpcNClsShared() 	                    // Number of shared TPC clusters
  //    ,track.trdPattern() 	                      // Contributor to the track on TRD layer in bits 0-5, starting from the innermost, bit 6 indicates a potentially split tracklet, bit 7 if the track crossed a padrow
  //    ,track.itsChi2NCl() 	                      // Chi2 / cluster for the ITS track segment
  //    ,track.tpcChi2NCl() 	                      // Chi2 / cluster for the TPC track segment
  //    ,track.trdChi2() 	                          // Chi2 for the TRD track segment
  //    ,track.tofChi2() 	                          // Chi2 for the TOF track segment
  //    ,track.tpcSignal() 	                        // dE/dx signal in the TPC
  //    ,track.trdSignal() 	                        // PID signal in the TRD
  //    ,track.length() 	                          // Track length
  //    ,track.tofExpMom() 	                        // TOF expected momentum obtained in tracking, used to compute the expected times
  //    ,track.trackEtaEmcal() 	                    // 
  //    ,track.trackPhiEmcal() 	                    // 
  //    ,track.trackTime() 	                        // Estimated time of the track in ns wrt collision().bc() or ambiguoustrack.bcSlice()[0]
  //    ,track.trackTimeRes() 	                    // Resolution of the track time in ns (see TrackFlags::TrackTimeResIsRange)
     
  //    // ,track.tpcDeltaTFwd() 	                    // Delta Forward of track time in TPC time bis
  //    // ,track.tpcDeltaTBwd() 	                    // Delta Backward of track time in TPC time bis
  //    ,track.pidForTracking() 	                  // PID hypothesis used during tracking. See the constants in the class PID in PID.h
  //    ,track.isPVContributor() 	                  // Run 3: Has this track contributed to the collision vertex fit
  //    ,track.hasITS() 	                          // Flag to check if track has a ITS match
  //    ,track.hasTPC() 	                          // Flag to check if track has a TPC match
  //    ,track.hasTRD() 	                          // Flag to check if track has a TRD match
  //    ,track.hasTOF() 	                          // Flag to check if track has a TOF measurement
  //    ,track.tpcNClsFound() 	                    // Number of found TPC clusters
  //    ,track.tpcNClsCrossedRows() 	              // Number of crossed TPC Rows
  //    ,track.itsClusterMap() 	                    // ITS cluster map, one bit per a layer, starting from the innermost
  //    ,track.itsNCls() 	                          // Number of ITS clusters
  //    ,track.itsNClsInnerBarrel() 	              // Number of ITS clusters in the Inner Barrel
  //    // ,track.itsClsSizeInLayer() 	                // Size of the ITS cluster in a given layer
  //    // ,track.isITSAfterburner() 	                // If the track used the afterburner in the ITS
  //    // ,track.tofExpTime() 	                      // Expected time for the track to reach the TOF
  //    // ,track.tofExpTimePi() 	                    // Expected time for the track to reach the TOF under the pion hypothesis
  //    // ,track.tofExpTimeKa() 	                    // Expected time for the track to reach the TOF under the kaon hypothesis
  //    // ,track.tofExpTimePr() 	                    // Expected time for the track to reach the TOF under the proton hypothesis
  //    ,track.tpcCrossedRowsOverFindableCls() 	    // Ratio crossed rows over findable clusters
  //    ,track.tpcFoundOverFindableCls()            // Ratio of found over findable clusters
  //    ,track.tpcFractionSharedCls() 	            // Fraction of shared TPC clusters
  //    ,track.detectorMap() 	                      // Detector map version 1, see enum DetectorMapEnum
   
  //    // TracksQA
  //    // // ,o2::soa::Index 	GI 	globalIndex 	    int64_t 	
  //    // // ,o2::aod::trackqa::TrackId 	            //I 	    trackId 	        //int32 	track to which this QA information belongs
  //    ,tpcTime0             //tpc only time0 (mTime0 in TPC track)
  //    ,tpcdcaR              //tpc only DCAr
  //    ,tpcdcaZ              //tpc only DCAz
  //    ,tpcClusterByteMask   //tracklet bitmask - track defining 8 tracklets (152=8*19 rows) bit set if nCluster>thr (default 5)
  //    ,tpcdEdxMax0R         //TPC dEdxQMax -ROC0/dEdx
  //    ,tpcdEdxMax1R         //TPC dEdxQMax -ROC1/dEdx
  //    ,tpcdEdxMax2R         //TPC dEdxQMax -ROC2/dEdx
  //    ,tpcdEdxMax3R         //TPC dEdxQMax -ROC3/dEdx
  //    ,tpcdEdxTot0R         //TPC dEdxQtot -ROC0/dEdx
  //    ,tpcdEdxTot1R         //TPC dEdxQtot -ROC1/dEdx
  //    ,tpcdEdxTot2R         //TPC dEdxQtot -ROC2/dEdx
  //    ,tpcdEdxTot3R         //TPC dEdxQtot -ROC3/dEdx
  //    //
  //    ,track.beta()
  //    ,track.betaerror()
  //    ,track.mass()

  //    ,track.tpcNSigmaPi()
  //    ,track.tofNSigmaPi()
  //    ,track.tpcNSigmaKa()
  //    ,track.tofNSigmaKa()
  //    ,track.tpcNSigmaPr()
  //    ,track.tofNSigmaPr()
  //    ,track.tpcNSigmaEl()
  //    ,track.tofNSigmaEl()
  //    ,track.tpcNSigmaDe()
  //    ,track.tofNSigmaDe()
     
  //    ,track.tofSignal()
  //    // ,track.eventCollisionTime()
  //    ,track.goodTOFMatch()

  //    ,track.tofFlags()        // o2::aod::pidflags::TOFFlags		     //  tofFlags	uint8_t	Flag for the complementary TOF PID information for the event time
  //    ,track.isEvTimeDefined() //  D	isEvTimeDefined	bool	True if the Event Time was computed with any method i.e. there is a usable event time
  //    ,track.isEvTimeTOF()	   //  D	isEvTimeTOF	bool	True if the Event Time was computed with the TOF
  //    ,track.isEvTimeT0AC()	   //  D	isEvTimeT0AC	bool	True if the Event Time was computed with the T0AC
  //    ,track.isEvTimeTOFT0AC() //  

  //    ,time
  //    ,TFidThis
  //    ,bcInTF

  //    ,trackMeanOcc.trackId()

  //    ,trackMeanOcc.meanOccPrimUnfm80() 
  //    ,trackMeanOcc.meanOccFV0AUnfm80() 
  //    ,trackMeanOcc.meanOccFV0CUnfm80() 
  //    ,trackMeanOcc.meanOccFT0AUnfm80() 
  //    ,trackMeanOcc.meanOccFT0CUnfm80() 
  //    ,trackMeanOcc.meanOccFDDAUnfm80() 
  //    ,trackMeanOcc.meanOccFDDCUnfm80() 
  //    ,trackMeanOcc.meanOccNTrackITSUnfm80()   
  //    ,trackMeanOcc.meanOccNTrackTPCUnfm80()   
  //    ,trackMeanOcc.meanOccNTrackTRDUnfm80()   
  //    ,trackMeanOcc.meanOccNTrackTOFUnfm80()   
  //    ,trackMeanOcc.meanOccNTrackSizeUnfm80()   
  //    ,trackMeanOcc.meanOccNTrackTPCAUnfm80()  
  //    ,trackMeanOcc.meanOccNTrackTPCCUnfm80()  
  //    ,trackMeanOcc.meanOccNTrackITSTPCUnfm80()
  //    ,trackMeanOcc.meanOccNTrackITSTPCAUnfm80()
  //    ,trackMeanOcc.meanOccNTrackITSTPCCUnfm80()
  //    ,trackMeanOcc.meanOccMultNTracksHasITSUnfm80()
  //    ,trackMeanOcc.meanOccMultNTracksHasTPCUnfm80()
  //    ,trackMeanOcc.meanOccMultNTracksHasTOFUnfm80()
  //    ,trackMeanOcc.meanOccMultNTracksHasTRDUnfm80()
  //    ,trackMeanOcc.meanOccMultNTracksITSOnlyUnfm80()
  //    ,trackMeanOcc.meanOccMultNTracksTPCOnlyUnfm80()
  //    ,trackMeanOcc.meanOccMultNTracksITSTPCUnfm80()
  //    ,trackMeanOcc.meanOccMultAllTracksTPCOnlyUnfm80()
  //    ,trackMeanOcc.meanOccRobustT0V0PrimUnfm80()   
  //    ,trackMeanOcc.meanOccRobustFDDT0V0PrimUnfm80()
  //    ,trackMeanOcc.meanOccRobustNtrackDetUnfm80()  
  //    ,trackMeanOcc.meanOccRobustMultExtraTableUnfm80()  
  //    ,trackMeanOcc.weightMeanOccPrimUnfm80() 
  //    ,trackMeanOcc.weightMeanOccFV0AUnfm80() 
  //    ,trackMeanOcc.weightMeanOccFV0CUnfm80() 
  //    ,trackMeanOcc.weightMeanOccFT0AUnfm80() 
  //    ,trackMeanOcc.weightMeanOccFT0CUnfm80() 
  //    ,trackMeanOcc.weightMeanOccFDDAUnfm80() 
  //    ,trackMeanOcc.weightMeanOccFDDCUnfm80() 
  //    ,trackMeanOcc.weightMeanOccNTrackITSUnfm80()   
  //    ,trackMeanOcc.weightMeanOccNTrackTPCUnfm80()   
  //    ,trackMeanOcc.weightMeanOccNTrackTRDUnfm80()   
  //    ,trackMeanOcc.weightMeanOccNTrackTOFUnfm80()   
  //    ,trackMeanOcc.weightMeanOccNTrackSizeUnfm80()   
  //    ,trackMeanOcc.weightMeanOccNTrackTPCAUnfm80()  
  //    ,trackMeanOcc.weightMeanOccNTrackTPCCUnfm80()  
  //    ,trackMeanOcc.weightMeanOccNTrackITSTPCUnfm80()
  //    ,trackMeanOcc.weightMeanOccNTrackITSTPCAUnfm80()
  //    ,trackMeanOcc.weightMeanOccNTrackITSTPCCUnfm80()
  //    ,trackMeanOcc.weightMeanOccMultNTracksHasITSUnfm80()
  //    ,trackMeanOcc.weightMeanOccMultNTracksHasTPCUnfm80()
  //    ,trackMeanOcc.weightMeanOccMultNTracksHasTOFUnfm80()
  //    ,trackMeanOcc.weightMeanOccMultNTracksHasTRDUnfm80()
  //    ,trackMeanOcc.weightMeanOccMultNTracksITSOnlyUnfm80()
  //    ,trackMeanOcc.weightMeanOccMultNTracksTPCOnlyUnfm80()
  //    ,trackMeanOcc.weightMeanOccMultNTracksITSTPCUnfm80()
  //    ,trackMeanOcc.weightMeanOccMultAllTracksTPCOnlyUnfm80()
  //    ,trackMeanOcc.weightMeanOccRobustT0V0PrimUnfm80()   
  //    ,trackMeanOcc.weightMeanOccRobustFDDT0V0PrimUnfm80()
  //    ,trackMeanOcc.weightMeanOccRobustNtrackDetUnfm80()  
  //    ,trackMeanOcc.weightMeanOccRobustMultExtraTableUnfm80()
     
  //    ,track.dcaXY()   //	dcaXY	float	Impact parameter in XY of the track to the primary vertex
  //    ,track.dcaZ()	   //	dcaZ	float	Impact parameter in Z of the track to the primary vertex
  //    ,track.sigmaDcaXY2()  //		sigmaDcaXY2	float	Impact parameter sigma^2 in XY of the track to the primary vertex
  //    ,track.sigmaDcaZ2()   //		sigmaDcaZ2	float	Impact parameter sigma^2 in Z of the track to the primary verte

  //    ,track.usedForTOFEvTime()
  //    ,track.evTimeTOF()
  //    ,track.evTimeTOFErr()
  //    ,track.evTimeTOFMult()

  //    ,track.tofExpTimeEl() //D	tofExpTimeEl	float	Expected time for the track to reach the TOF under the electron hypothesis
  //    ,track.tofExpTimeMu() //D	tofExpTimeMu	float	Expected time for the track to reach the TOF under the muon hypothesis
  //    ,track.tofExpTimePi() //D	tofExpTimePi	float	Expected time for the track to reach the TOF under the pion hypothesis
  //    ,track.tofExpTimeKa() //D	tofExpTimeKa	float	Expected time for the track to reach the TOF under the kaon hypothesis
  //    ,track.tofExpTimePr() //D	tofExpTimePr	float	Expected time for the track to reach the TOF under the proton hypothesis
  //    ,track.tofExpTimeDe() //D	tofExpTimeDe	float	Expected time for the track to reach the TOF under the deuteron hypothesis
  //    ,track.tofExpTimeTr() //D	tofExpTimeTr	float	Expected time for the track to reach the TOF under the triton hypothesis
  //    ,track.tofExpTimeHe() //D	tofExpTimeHe	float	Expected time for the track to reach the TOF under the helium3 hypothesis
  //    ,track.tofExpTimeAl() //D	tofExpTimeAl

  //    //o2::aod::track::v001::ITSClsSizeInLayer	D	itsClsSizeInLayer	uint8_t	Size of the ITS cluster in a given layer
  //    // ,track.itsClsSizeInLayer() //	uint8_t	Size of the ITS cluster in a given layer

  //    ,track.itsNSigmaEl()
  //    ,track.itsNSigmaMu()
  //    ,track.itsNSigmaPi()
  //    ,track.itsNSigmaKa()
  //    ,track.itsNSigmaPr()
  //    ,track.itsNSigmaDe()
  //    ,track.itsNSigmaTr()
  //    ,track.itsNSigmaHe()
  //    ,track.itsNSigmaAl()

  //    ,deltaRefContParamY   
  //    ,deltaRefContParamZ   
  //    ,deltaRefContParamSnp 
  //    ,deltaRefContParamTgl 
  //    ,deltaRefContParamQ2Pt
  //    ,deltaRefGloParamY    
  //    ,deltaRefGloParamZ    
  //    ,deltaRefGloParamSnp  
  //    ,deltaRefGloParamTgl  
  //    ,deltaRefGloParamQ2Pt 
  //    ,deltaTOFdX
  //    ,deltaTOFdZ           
  //    ,weightDSPt
  //    ,weightDSQPt
  //    ,triggerMaskDS
  //  );


    LOG(info)<<"DEBUG ::";

    gCollisionChecker += collisions.size();
    gFV0AChecker      += FV0As.size();
    gFT0Checker       += FT0s.size();
    gTrackChecker     += (trackQAcounter+1);
    gV0Checker        += (v0Counter);

  }//Process function ends

};

/*
Data will be sampled to be uniform if the PDF.
We will start with the flat pt ttigger and flat q/pt trigger and MB trigger (random downsampling)

Algorithm:
* 1D histogram of the pt  for the TPC+ITS tracks  - hisPt(250,0.1,50)
* 1D hsitgram of the q/pt for the TPC+ITS tracks  - hisQPt(200,-6,6)
* 1D histogram of the pt  for the K0s  - hisPtK0(250,0.1,50)
* 1D histogram of the pt  for the Lambda  - hisPtLambda(250,0.1,50)
* 1D histogram of the pt  for the Gamma  - hisPtGamma(250,0.1,50)
* 1D histogram of the pt  for the Phi  - hisPtPhi(250,0.1,50)
* 1D histogram of the pt  for the D0  - hisPtD0(250,0.1,50)

If V0 is triggered both daugher tracks are triggers  (mask to add to the track trigger mask)
* isK0Trigger, isLambdTrigger, isGammaTriggers,isPhiTriggered


```
bool isTriggeredPt(float pt, float factorPt, float & weight){
    float prob= hisPt->GetValue(pt);
    float prob1=hisPt->GetValue(1);
    bool isTriggered = gRandom->Rndm()*prob/prob1<factorPt;
    weight=prob;
    return isTriggered;
}
```
```
bool isTriggeredQPt(float pt, float factorQPt, float & weight){
    float prob= hisQPt->GetValue(pt);
    float prob1=hisQPt->GetValue(1);
    bool isTriggered = gRandom->Rndm()*prob/prob1<factorQPt;
    weight=prob;
    return isTriggered;
}
```
```c++
bool isDownsampledMB =gRandom->Rndm()<factorMB;
bool isDownsamledHPT=isTriggeredPt(pt,factorPt);
bool isDownsampledQPT=isTriggeredQPt(pt,factorPt);
int triggermask=(1<<0)*isDownsampledMB+(1<<1)*isDownsamledHPT+(1<<2)*isDownsamledQPT;
if triggerasmk>0 track written
triggermask to be written and weight should be written

```


In the old code with analytical PDF
```
Int_t  DownsampleTsalisCharged(Double_t pt, Double_t factorPt, Double_t factor1Pt, Double_t sqrts, Double_t mass, Double_t *weight){
  Double_t prob=TsalisCharged(pt,mass,sqrts);
  Double_t probNorm=TsalisCharged(1.,mass,sqrts);
  Int_t triggerMask=0;
  (*weight)=prob/probNorm;
  if (gRandom->Rndm()*prob/probNorm<factorPt) triggerMask|=1;
  if ((gRandom->Rndm()*((prob/probNorm)*pt*pt))<factor1Pt) triggerMask|=2;
  if (gRandom->Rndm()<factorPt) triggerMask|=4;
  return triggerMask;
}
```
*/

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc){
  return WorkflowSpec{ 
                      adaptAnalysisTask<trackqapidderivedata>(cfgc)
                      // ,adaptAnalysisTask<occupancyTableConsumer>(cfgc)
  };
}
