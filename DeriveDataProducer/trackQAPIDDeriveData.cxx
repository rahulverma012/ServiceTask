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

/// https://its.cern.ch/jira/browse/O2-5095
/// https://gitlab.cern.ch/alice-tpc-offline/alice-tpc-notes/-/tree/master/JIRA/O2-4592

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
#include <TLorentzVector.h>
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

#define BITSET(mask, ithBit) ((mask) |= (1 << (ithBit)))  // avoid name bitset as std::bitset is already there
#define BITCHECK(mask, ithBit) ((mask) & (1 << (ithBit))) // bit check will return int value of (1<<ithBit), not bool, use BITCHECK != 0 in Analysis
#define BITCHECKOTHER(mask, ithBit) (((mask) & ~(1 << (ithBit))) != 0) //checks if any of the other bit execpt ithBit is non zero or not

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

//_________________________________________________________________________________________________
/*
.L $NOTES/JIRA/Issue14289/lossyCompressionUtils.C

 */

// lossyCompressionUtils.C

#include "lossyCompressionUtils.h"

#include <cmath>
#include <map>
#include <memory>
#include <random>
#include <tuple>
#include <algorithm>
#include <cassert>
#include <unordered_map>

namespace lossy {

inline int clamp(int val, int minVal, int maxVal) {
    return std::max(minVal, std::min(maxVal, val));
}

void SqrtScalingLUT::build(float sigma0_, float sigma1_, float xmin_, float xmax_,
                           int clampMin_, int clampMax_, int nBins_) {
    sigma0 = sigma0_;
    sigma1 = sigma1_;
    xmin = xmin_;
    xmax = xmax_;
    nBins = nBins_;
    clampMin = clampMin_;
    clampMax = clampMax_;

    codeToVal.resize(nBins);
    valToCode.resize(nBins);

    float step = (xmax - xmin) / (nBins - 1);
    for (int i = 0; i < nBins; ++i) {
        float val = xmin + i * step;
        float code_f = std::asinh((sigma1 * val) / sigma0) / sigma0;
        int code = clamp(static_cast<int>(std::round(code_f)), clampMin, clampMax);
        codeToVal[i] = val;
        valToCode[i] = code;
    }
}

int SqrtScalingLUT::encode(float val) const {
    if (val <= xmin) return valToCode.front();
    if (val >= xmax) return valToCode.back();
    float binWidth = (xmax - xmin) / (nBins - 1);
    int idx = static_cast<int>((val - xmin) / binWidth);
    return valToCode[idx];
}

float SqrtScalingLUT::decode(int code, bool addRndm) const {
    float code_f = static_cast<float>(code);
    if (addRndm) {
        static thread_local std::mt19937 gen(std::random_device{}());
        static thread_local std::uniform_real_distribution<float> dist(-0.5f, 0.5f);
        code_f += dist(gen);
    }
    return (sigma0 / sigma1) * std::sinh(sigma0 * code_f);
}

int SqrtScalingLUT::bitRange() const {
    return static_cast<int>(std::ceil(std::log2(clampMax - clampMin + 1)));
}

double SqrtScalingLUT::estimateEntropy() const {
    std::unordered_map<int, int> histogram;
    for (int code : valToCode) {
        ++histogram[code];
    }
    double entropy = 0.0;
    double total = static_cast<double>(valToCode.size());
    for (const auto& [code, count] : histogram) {
        double p = static_cast<double>(count) / total;
        entropy -= p * std::log2(p);
    }
    return entropy;
}

// --- Global cache for LUTs ---
static std::map<std::tuple<float, float, int, int>, std::shared_ptr<SqrtScalingLUT>> lutCache;

int codeSqrtScaling(float val, float sigma0, float sigma1,
                    int clampMin, int clampMax, QuantMode mode) {
    if (mode == QuantMode::Analytical) {
        float code_f = std::asinh((sigma1 * val) / sigma0) / sigma0;
        return clamp(static_cast<int>(std::round(code_f)), clampMin, clampMax);
    } else {
        auto key = std::make_tuple(sigma0, sigma1, clampMin, clampMax);
        auto& lutPtr = lutCache[key];
        if (!lutPtr) {
            auto lut = std::make_shared<SqrtScalingLUT>();
            float range = 1000.0f; // typical symmetric range for LUT
            lut->build(sigma0, sigma1, -range, range, clampMin, clampMax);
            lutPtr = lut;
        }
        return lutPtr->encode(val);
    }
}

float decodeSqrtScaling(int code, float sigma0, float sigma1,
                        bool addRndm, int clampMin, int clampMax, QuantMode mode) {
    if (mode == QuantMode::Analytical) {
        float code_f = static_cast<float>(code);
        if (addRndm) {
            static thread_local std::mt19937 gen(std::random_device{}());
            static thread_local std::uniform_real_distribution<float> dist(-0.5f, 0.5f);
            code_f += dist(gen);
        }
        return (sigma0 / sigma1) * std::sinh(sigma0 * code_f);
    } else {
        auto key = std::make_tuple(sigma0, sigma1, clampMin, clampMax);
        auto& lutPtr = lutCache[key];
        if (!lutPtr) {
            auto lut = std::make_shared<SqrtScalingLUT>();
            float range = 10.0f;
            lut->build(sigma0, sigma1, -range, range, clampMin, clampMax);
            lutPtr = lut;
        }
        return lutPtr->decode(code, addRndm);
    }
}

} // namespace lossy
//_________________________________________________________________________________________________
using namespace lossy;


//My functions 
inline double deltaAngleCalculator(double aA, double tA) {
  constexpr double twoPi = 6.283185307179586476925286766559; // 2*pi
  constexpr double pi = 3.1415926535897932384626433832795;  // pi

  //transform range from [-pi, pi] to [0,2pi]
  aA = (aA < 0.0) ? aA + twoPi : aA; //aA = associated Particle Angle
  tA = (tA < 0.0) ? tA + twoPi : tA; //tA = trigger Particle Angle
  double dA = aA - tA; //dA = delta Angle

  // Wrap into [-π, π]
  if (dA > pi) dA -= twoPi; //dA was in range [pi, 2pi]
  else if (dA < -pi) dA += twoPi; //dA was in range [-2pi, -pi]
  return dA; //dA in range [-pi, pi]
}

template <typename T, std::size_t N>
void sortVectorOfArray(std::vector<std::array<T, N>>& myVector, const int& myIDX)
{
  std::sort(myVector.begin(), myVector.end(), [myIDX](const std::array<T, N>& a, const std::array<T, N>& b) {
    return a[myIDX] < b[myIDX]; // sort at the required index
  });
}

template <typename T, std::size_t N>
void checkUniqueness(const std::vector<std::array<T, N>>& myVector, const int& myIDX)
{
  for (size_t i = 1; i < myVector.size(); i++) {
    if (myVector[i][myIDX] <= myVector[i - 1][myIDX]) {
      LOG(error) << "Duplicate Entries while creating Index tables :: (vec[" << i << "][" << myIDX << "]) " << myVector[i][myIDX] << " >= " << myVector[i - 1][myIDX] << " (vec[" << i - 1 << "][" << myIDX << "])";
    }
  }
}    

template<typename T>
inline bool binary_search_globalIndex(const std::vector<T>& vec, int64_t key, int64_t& idMethod, uint16_t& flag) {
    size_t left = 0;
    size_t right = vec.size();
    while (left < right) {
        size_t mid = left + (right - left) / 2;
        int64_t midVal = vec[mid].globalIndex;
        if (midVal == key) {
            idMethod = vec[mid].idMethod;
            flag = vec[mid].flag;
            return true;
        }
        else if (midVal < key) {
            left = mid + 1;
        } else {
            right = mid;
        }
    }
    return false;
}

struct trackqapidderivedata{

  Service<o2::ccdb::BasicCCDBManager> ccdb;

  // declare production of tables
  Produces<aod::DrCollisions> DerivedCollisions;
  Produces<aod::DrTracks0> DerivedTracks0;
  Produces<aod::DrTracks> DerivedTracks;
  Produces<aod::DrFV0As> DerivedFV0As;
  Produces<aod::DrFT0s> DerivedFT0s;
  Produces<aod::DrV0s> DerivedV0s;
  Produces<aod::DrPhis> DerivedPhis;
  Produces<aod::DrD0s> DerivedD0s;
  Produces<aod::DrZ0s> DerivedZ0s;
  Produces<aod::DrCosmicPairs> DerivedCosmicPairs;
  Produces<aod::OCCS> BuildOCCS;
  Produces<aod::DrTracksMothers> DerivedDrTracksMothers;

  //Histogram registry;
  HistogramRegistry recoTracks  {"recoTracks"  , {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry recoEvent   {"recoEvent"   , {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  //Configurables
  
  Configurable<int> customOrbitOffset{"customOrbitOffset", 0, "customOrbitOffset for MC"};
  Configurable<float> cfgDSMBfactorTracks{"cfgDSMBfactorTracks", 0.002, "cfgDSMBfactorTracks"};
  Configurable<float> cfgDSPtfactorTracks{"cfgDSPtfactorTracks", 0.002, "cfgDSPtfactorTracks"};

  Configurable<float> cfgDSMBfactorV0{"cfgDSMBfactorV0", 0.02, "cfgDSMBfactorV0"};
  Configurable<float> cfgDSPtfactorV0{"cfgDSPtfactorV0", 0.02, "cfgDSPtfactorV0"};

  Configurable<float> cfgDSMBfactorD0s{"cfgDSMBfactorD0s", 0.01, "cfgDSMBfactorD0s"};
  Configurable<float> cfgDSPtfactorD0s{"cfgDSPtfactorD0s", 0.01, "cfgDSPtfactorD0s"};

  Configurable<float> cfgDSMBfactorPhis{"cfgDSMBfactorPhis", 0.01, "cfgDSMBfactorPhis"};
  Configurable<float> cfgDSPtfactorPhis{"cfgDSPtfactorPhis", 0.01, "cfgDSPtfactorPhis"};


  Configurable<bool> cfgDoSkipDF{"cfgDoSkipDF", false, "cfgDoSkipDF"};
  Configurable<int > cfgNSkipDF {"cfgNSkipDF" , 3    , "cfgNSkipDF"};

  Configurable<bool> cfgDoDSOfTrack{"cfgDoDSOfTrack", false, "cfgDoDSOfTrack"};
  Configurable<bool> cfgDoDSOfV0s  {"cfgDoDSOfV0s"  , false, "cfgDoDSOfV0s"};
  Configurable<bool> cfgDoDSOfD0s  {"cfgDoDSOfD0s"  , false, "cfgDoDSOfD0s"};
  Configurable<bool> cfgDoDSOfPhi  {"cfgDoDSOfPhi"  , false, "cfgDoDSOfPhi"};

  Configurable<bool> cfgDoTableOfTrack{"cfgDoTableOfTrack", false, "cfgDoTableOfTrack"};
  Configurable<bool> cfgDoTableOfV0s  {"cfgDoTableOfV0s"  , false, "cfgDoTableOfV0s"};
  Configurable<bool> cfgDoTableOfD0s  {"cfgDoTableOfD0s"  , false, "cfgDoTableOfD0s"};
  Configurable<bool> cfgDoTableOfPhi  {"cfgDoTableOfPhi"  , false, "cfgDoTableOfPhi"};


  Configurable<bool> doV0andTrackMatchingCheck{"doV0andTrackMatchingCheck", true, "doV0andTrackMatchingCheck"};
  Configurable<bool> doD0andTrackMatchingCheck{"doD0andTrackMatchingCheck", true, "doD0andTrackMatchingCheck"};

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
  
  struct : ConfigurableGroup {
  Configurable<float> csmMuThrPt        {"csmMuThrPt"        , 0.5, "csmMuThrPt"};
  Configurable<float> csmMuThrDCAMin    {"csmMuThrDCAMin"    , 0.5, "csmMuThrDCAMin"};
  Configurable<float> csmMuThrDCAMax    {"csmMuThrDCAMax"    , 120, "csmMuThrDCAMax"};

  Configurable<float> csmMuSumPtPair    {"csmMuSumPtPair"    , 2.0 , "csmMuSumPtPair"};
  Configurable<float> csmMuSumQPtPair   {"csmMuSumQPtPair"   , 0.2 , "csmMuSumQPtPair"};
  Configurable<float> csmMuSumTglPair   {"csmMuSumTglPair"   , 0.01, "csmMuSumTglPair"};
  Configurable<float> csmMuSumDcaXY     {"csmMuSumDcaXY"     , 0.6 , "csmMuSumDcaXY"};
  Configurable<float> csmMuDiffDcaXY    {"csmMuDiffDcaXY"    , 0.6 , "csmMuDiffDcaXY"};
  Configurable<float> csmMuDiffAlphaPair{"csmMuDiffAlphaPair", 0.01, "csmMuDiffAlphaPair"};

  Configurable<bool> csmMuCheckSumDcaXY     {"csmMuCheckSumDcaXY"     , true, "csmMuCheckSumDcaXY"};
  Configurable<bool> csmMuCheckDiffDcaXY    {"csmMuCheckDiffDcaXY"    , false, "csmMuCheckDiffDcaXY"};
  Configurable<bool> csmMuCheckDiffAlphaPair{"csmMuCheckDiffAlphaPair", true, "csmMuCheckDiffAlphaPair"};
  } cfgCM;

  Configurable<bool> fillTracksFullInfo {"fillTracksFullInfo", false, "fillTracksFullInfo"};

  enum CosmicPairRejectionType{
    kPairPassed = 0,
    kFailSumPtPair,    
    kFailSumQPtPair,   
    kFailSumTglPair,   
    kFailSumDcaXY,     
    kFailDiffDcaXY,    
    kFailDiffAlphaPair
  };

  static constexpr std::string_view CMRejectionTag[]{
    "kPairPassed",
    "kFailSumPtPair",    
    "kFailSumQPtPair",   
    "kFailSumTglPair",   
    "kFailSumDcaXY",     
    "kFailDiffDcaXY",    
    "kFailDiffAlphaPair"
  };
  
  void init(InitContext const&){
    //
    LOGF(info,"Starting init");

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

      recoEvent.add("CosmicMuon/PreSel/sumPtPair"    ,"sumPtPair"   , kTH1F, {{510, -1., 50.}});
      recoEvent.add("CosmicMuon/PreSel/sumQPtPair"   ,"sumQPtPair"  , kTH1F, {{510, -1., 50.}});
      recoEvent.add("CosmicMuon/PreSel/sumTglPair"   ,"sumTglPair"  , kTH1F, {{510, -1., 50.}});
      recoEvent.add("CosmicMuon/PreSel/sumDcaXYPair" ,"sumDcaXYPair", kTH1F, {{510, -1., 50.}});
      recoEvent.add("CosmicMuon/PreSel/sumAlphaPair" ,"sumAlphaPair", kTH1F, {{510, -1., 50.}});

      recoEvent.add("CosmicMuon/PreSel/diffPtPair"   , "diffPtPair"   , kTH1F, {{510, -25., 26.}});
      recoEvent.add("CosmicMuon/PreSel/diffQPtPair"  , "diffQPtPair"  , kTH1F, {{510, -25., 26.}});
      recoEvent.add("CosmicMuon/PreSel/diffTglPair"  , "diffTglPair"  , kTH1F, {{510, -25., 26.}});
      recoEvent.add("CosmicMuon/PreSel/diffDcaXYPair", "diffDcaXYPair", kTH1F, {{510, -25., 26.}});
      recoEvent.add("CosmicMuon/PreSel/diffAlphaPair", "diffAlphaPair", kTH1F, {{510, -25., 26.}});

      recoEvent.addClone("CosmicMuon/PreSel/", "CosmicMuon/PostSel/");
      recoEvent.add("CosmicMuon/PairRejectionTrigger", "PairRejectionTrigger", kTH1F, {{11*2, -1,10}});
      for(int i = 0 ; i <= kFailDiffAlphaPair; i++){
        int iBin = recoEvent.get<TH1>(HIST("CosmicMuon/PairRejectionTrigger"))->FindBin(i);
        recoEvent.get<TH1>(HIST("CosmicMuon/PairRejectionTrigger"))->GetXaxis()->SetBinLabel(iBin, CMRejectionTag[i].data());
      }

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

  // ------------------------------  Lossy Compression [END]  ------------------------------------- //

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
    bcSOR = orbitSOR * nBCsPerOrbit + customOrbitOffset * nBCsPerOrbit; // customOrbitOffset is a Configurable
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
    // Global Track Selection
    if (!track.isGlobalTrack()) {
      return false;
    }

    // Particle Identification
    if(track.hasTOF()){ //if track has TOF signal
      if( (std::pow(track.tpcNSigmaKa(),2) + std::pow(track.tofNSigmaKa(),2)) < 4.0){ // changing nsigma cut to 2 (3 is too loose !!!)
        return true;
      } else {
        return false;
      }
    } else { // No TOF signal, select from TPC (We may need to further apply pT dependent selection for PbPb)
      // Select Kaon and Reject Pion and Proton Hypothesis
      // We may need to remove electrons as well
      if(std::abs(track.tpcNSigmaKa()) < 2.0 && std::abs(track.tpcNSigmaPi()) > 3.0 && std::abs(track.tpcNSigmaPr()) > 3.0) {
        return true;
      }else {
        return false;
      }
    }
  }

  template<typename T> 
  bool selElectron(const T& track){
    if(track.hasTOF()){ 
      if( std::sqrt(std::pow(track.tpcNSigmaEl(),2) + std::pow(track.tofNSigmaEl(),2)) < 2.0 ) {return true;}
    }
    else {
      if (std::abs(track.tpcNSigmaEl()) < 2.0){
        return true;
      }
    }
    return false; 
  }

  template<typename T> 
  bool selMuon(const T& track){
    if(track.hasTOF()){ 
      if( std::sqrt(std::pow(track.tpcNSigmaMu(),2) + std::pow(track.tofNSigmaMu(),2)) < 2.0 ) {return true;}
    }
    else {
      if (std::abs(track.tpcNSigmaMu()) < 2.0){
        return true;
      }
    }
    return false; 
  }

  using myCollisionsWithoutCent = soa::Join<aod::Collisions, aod::Mults, aod::MultsExtra, aod::FT0sCorrected, aod::EvSels>;

  using myCollisionsWithCent = soa::Join<aod::Collisions, aod::Mults, aod::MultsExtra, aod::FT0sCorrected, aod::EvSels
                                ,aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::CentFDDMs, aod::CentNTPVs>;

  using myTracks = soa::Join<aod::Tracks, aod::TrackToTracksQA, aod::TrackToAmbgTrk, o2::aod::TracksCov, aod::TracksExtra, aod::TracksDCA, aod::TracksDCACov, aod::TrackSelection //, aod::TracksQA>;//, 
                    ,aod::TOFSignal,  aod::pidTOFbeta, aod::pidTOFmass, aod::EvTimeTOFOnly, aod::pidTOFFlags, aod::pidEvTimeFlags
                    ,aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullDe
                    ,aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullDe>;

  using myOccsTable = soa::Join<aod::OccsBCsList, aod::OccsPrim, aod::OccsMeanPrim, 
                    aod::OccsT0V0, aod::OccsMeanT0V0, aod::ORT0V0Prim, aod::OMRT0V0Prim,
                    aod::OccsFDD, aod::OccsMeanFDD, aod::ORFDDT0V0Prim, aod::OMRFDDT0V0Prim,
                    aod::OccsNTrackDet, aod::OccsMeanNTrkDet, aod::ORNtrackDet, aod::OMRNtrackDet,
                    aod::OccsMultExtra, aod::OccsMnMultExtra, aod::ORMultExtra, aod::OMRMultExtra>;

  using myTrackMeanOccs = soa::Join<aod::TmoTrackIds, aod::TmoPrim , aod::TmoT0V0 , aod::TmoFDD , aod::TmoNTrackDet , aod::TmoMultExtra , aod::TmoRT0V0Prim , aod::TmoRFDDT0V0Prim, aod::TmoRNtrackDet , aod::TmoRMultExtra,
                                                    aod::TwmoPrim, aod::TwmoT0V0, aod::TwmoFDD, aod::TwmoNTrackDet, aod::TwmoMultExtra, aod::TwmoRT0V0Prim, aod::TwmoRFDDT0V0Pri, aod::TwmoRNtrackDet, aod::TwmoRMultExtra>;
                                                         
  using myBCTable = soa::Join<aod::BCsWithTimestamps,aod::OccIndexTable>;
  using MyTracksQA = aod::TracksQAVersion; //aod::TracksQA_003;

  // For manual sliceBy
  Preslice<myTracks>  TracksPerCollisionPreslice = o2::aod::track::collisionId;
  Preslice<aod::TracksQA>  trackQA_Preslice = o2::aod::trackqa::trackId;
  // Preslice<myTracks>  Tracks_PreSlice = o2::aod::track::globalIndex;
  // Preslice<myTracks>  Tracks_PreSlice = o2::soa::globalIndex;

  //Use Partition after definition of filtered Tracks
  SliceCache cache;
  Partition<myTracks> posTracks = aod::track::signed1Pt > 0.f;// track.sign() is dynamic column so use signed1Pt
  Partition<myTracks> negTracks = aod::track::signed1Pt < 0.f;


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
  uint64_t gZ0Checker =0;
  uint64_t gCosmicPairChecker = 0; 
  // int      dfCount = 0;
  int lastRun = -999;
  int64_t bcSOR = -1;
  int32_t nBCsPerTF = -1;
  uint64_t time = -1;
  int64_t TFidThis = -1;
  int bcInTF = -1;
  
  bool startDownsampling = false;

  float cfgOccSigma0 = 0.04, cfgOccSigma1 = 0.04;
  int  cfgOccClampMin = 0, cfgOccClampMax = 255;
  
  float cfgPIDSigma0 = 0.05, cfgPIDSigma1 = 0.05;
  int  cfgPIDClampMin = -128, cfgPIDClampMax = 128;

  std::chrono::high_resolution_clock::time_point Start0 = std::chrono::high_resolution_clock::now();

  template<bool processWithCent, typename B, typename C, typename T, typename A>
  void executeProcess(const B& BCs, const C& collisions, const T& tracks, const auto& tracksQA, const auto& Origins, 
      const A& ambgTracks, const auto &FV0As, const auto& FT0s, const auto &V0s, const auto& D0s, const auto& occIdxTable, 
      const auto& occTables, const auto &trackMeanOccs){
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

    float collisionCentFV0A = -999;
    float collisionCentFT0M = -999;
    float collisionCentFT0A = -999;
    float collisionCentFT0C = -999;
    float collisionCentFDDM = -999;
    float collisionCentNTPV = -999;

    std::vector<int64_t> TFIDList;
    for(const auto& collision : collisions){
      iColl++;
      const auto& bc = collision.template bc_as<myBCTable>();
      GetTimingInfo(bc, lastRun, nBCsPerTF, bcSOR, time, TFidThis, bcInTF);

      if (nBCsPerTF > nBCinTF) {
        LOG(error) << "DEBUG :: FATAL ERROR :: nBCsPerTF > nBCinTF i.e " << nBCsPerTF << " > " << nBCinTF << " will cause crash in further process";
        return;
      }
      
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

      if constexpr (processWithCent){ collisionCentFV0A = collision.centFV0A();}
      if constexpr (processWithCent){ collisionCentFT0M = collision.centFT0M();}
      if constexpr (processWithCent){ collisionCentFT0A = collision.centFT0A();}
      if constexpr (processWithCent){ collisionCentFT0C = collision.centFT0C();}
      if constexpr (processWithCent){ collisionCentFDDM = collision.centFDDM();}
      if constexpr (processWithCent){ collisionCentNTPV = collision.centNTPV();}

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

      ,collisionCentFV0A
      ,collisionCentFT0M
      ,collisionCentFT0A
      ,collisionCentFT0C
      ,collisionCentFDDM
      ,collisionCentNTPV
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

      isDownsampledMB[kV0]              = gRandom->Rndm()<cfgDSMBfactorV0;
      if( isLambda          ){isDownsampledMB[kLambda]          = gRandom->Rndm()<cfgDSMBfactorV0;}
      if( isAntiLambda      ){isDownsampledMB[kAntiLambda]      = gRandom->Rndm()<cfgDSMBfactorV0;}
      if( isK0Short         ){isDownsampledMB[kK0Short]         = gRandom->Rndm()<cfgDSMBfactorV0;}
      if( isGamma           ){isDownsampledMB[kGamma]           = gRandom->Rndm()<cfgDSMBfactorV0;}
      if( isHypertriton     ){isDownsampledMB[kHypertriton]     = gRandom->Rndm()<cfgDSMBfactorV0;}
      if( isAntiHypertriton ){isDownsampledMB[kAntiHypertriton] = gRandom->Rndm()<cfgDSMBfactorV0;}

      isDownsampledPT[kV0]              = false;
      isDownsampledPT[kLambda]          = false;
      isDownsampledPT[kAntiLambda]      = false;
      isDownsampledPT[kK0Short]         = false;
      isDownsampledPT[kGamma]           = false;
      isDownsampledPT[kHypertriton]     = false;
      isDownsampledPT[kAntiHypertriton] = false;

      isDownsampledPT[kV0]              = checkTriggerPt(v0.pt(), cfgDSPtfactorV0, weightDSPt[kV0]              , probRatioPt[kV0]              , recoEvent.get<TH1>(HIST("Event/V0/h_Pt"              )));
      if( isLambda          ){isDownsampledPT[kLambda]          = checkTriggerPt(v0.pt(), cfgDSPtfactorV0, weightDSPt[kLambda]          , probRatioPt[kLambda]          , recoEvent.get<TH1>(HIST("Event/Lambda/h_Pt"          )));}
      if( isAntiLambda      ){isDownsampledPT[kAntiLambda]      = checkTriggerPt(v0.pt(), cfgDSPtfactorV0, weightDSPt[kAntiLambda]      , probRatioPt[kAntiLambda]      , recoEvent.get<TH1>(HIST("Event/AntiLambda/h_Pt"      )));}
      if( isK0Short         ){isDownsampledPT[kK0Short]         = checkTriggerPt(v0.pt(), cfgDSPtfactorV0, weightDSPt[kK0Short]         , probRatioPt[kK0Short]         , recoEvent.get<TH1>(HIST("Event/K0Short/h_Pt"         )));}
      // if( isGamma           ){isDownsampledPT[kGamma]           = checkTriggerPt(v0.pt(), cfgDSPtfactorV0, weightDSPt[kGamma]           , probRatioPt[kGamma]           , recoEvent.get<TH1>(HIST("Event/Gamma/h_Pt"           )));}
      // if( isHypertriton     ){isDownsampledPT[kHypertriton]     = checkTriggerPt(v0.pt(), cfgDSPtfactorV0, weightDSPt[kHypertriton]     , probRatioPt[kHypertriton]     , recoEvent.get<TH1>(HIST("Event/Hypertriton/h_Pt"     )));}
      // if( isAntiHypertriton ){isDownsampledPT[kAntiHypertriton] = checkTriggerPt(v0.pt(), cfgDSPtfactorV0, weightDSPt[kAntiHypertriton] , probRatioPt[kAntiHypertriton] , recoEvent.get<TH1>(HIST("Event/AntiHypertriton/h_Pt" )));}
      
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
      if (isPhi   ) { isDownsampledMB[kPhi1020] = gRandom->Rndm()<cfgDSMBfactorPhis;}

      isDownsampledPT[kPhi1020] = false;
      if (isPhi   ) {isDownsampledPT[kPhi1020] = checkTriggerPt(phiPtlist[phiCand], cfgDSPtfactorPhis, weightDSPt[kPhi1020], probRatioPt[kPhi1020], recoEvent.get<TH1>(HIST("Event/Phi/h_Pt")));}

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
      if (isD0   ) { isDownsampledMB[kD0] = gRandom->Rndm()<cfgDSMBfactorD0s;}
      if (isD0Bar) { isDownsampledMB[kD0Bar] = gRandom->Rndm()<cfgDSMBfactorD0s;}

      isDownsampledPT[kD0] = false;
      isDownsampledPT[kD0Bar] = false;
      if (isD0   ) {isDownsampledPT[kD0] = checkTriggerPt(d0.pt(), cfgDSPtfactorD0s, weightDSPt[kD0], probRatioPt[kD0], recoEvent.get<TH1>(HIST("Event/D0/h_Pt")));}
      if (isD0Bar) {isDownsampledPT[kD0Bar] = checkTriggerPt(d0.pt(), cfgDSPtfactorD0s, weightDSPt[kD0Bar], probRatioPt[kD0Bar], recoEvent.get<TH1>(HIST("Event/D0Bar/h_Pt")));}

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

    std::vector<float> Occ_Prim_Unfm_80            ;
    std::vector<float> Occ_FV0A_Unfm_80            ;
    std::vector<float> Occ_FV0C_Unfm_80            ;
    std::vector<float> Occ_FT0A_Unfm_80            ;
    std::vector<float> Occ_FT0C_Unfm_80            ;
    std::vector<float> Occ_FDDA_Unfm_80            ;
    std::vector<float> Occ_FDDC_Unfm_80            ;

    std::vector<float> Occ_NTrack_PVC_Unfm_80      ;
    std::vector<float> Occ_NTrack_ITS_Unfm_80      ;
    std::vector<float> Occ_NTrack_TPC_Unfm_80      ;
    std::vector<float> Occ_NTrack_TRD_Unfm_80      ;
    std::vector<float> Occ_NTrack_TOF_Unfm_80      ;
    std::vector<float> Occ_NTrackSize_Unfm_80      ;
    std::vector<float> Occ_NTrackTPC_A_Unfm_80     ;
    std::vector<float> Occ_NTrackTPC_C_Unfm_80     ;
    std::vector<float> Occ_NTrackITS_TPC_Unfm_80   ;
    std::vector<float> Occ_NTrackITS_TPC_A_Unfm_80 ;
    std::vector<float> Occ_NTrackITS_TPC_C_Unfm_80 ;

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
          bc = collision.template bc_as<myBCTable>();
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
          bc = collision.template bc_as<myBCTable>();
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

      isDownsampledMB  = gRandom->Rndm()<cfgDSMBfactorTracks;
      isDownsampledPT  =  checkTriggerPt(track.pt()       , cfgDSPtfactorTracks, weightDSPt , probRatioPt , recoEvent.get<TH1>(HIST("Event/track/Full/h_Pt" )));
      isDownsampledQPT = checkTriggerQPt(track.signed1Pt(), cfgDSPtfactorTracks, weightDSQPt, probRatioQPT, recoEvent.get<TH1>(HIST("Event/track/Full/h_QPt")));

      // check If It is a V0 daughter and getItsV0s. 
      if(track.sign() > 0) { isDownsampledV0 = std::binary_search(posV0DauList.begin(), posV0DauList.end(), track.globalIndex());}
      else                 { isDownsampledV0 = std::binary_search(negV0DauList.begin(), negV0DauList.end(), track.globalIndex());}

      // // check if It is a prong of Phi
      // if(track.sign() > 0) { isDownsampledPhi = std::binary_search(posKaonPhiList.begin(), posKaonPhiList.end(), track.globalIndex());}
      // else                 { isDownsampledPhi = std::binary_search(negKaonPhiList.begin(), negKaonPhiList.end(), track.globalIndex());}

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

      if constexpr (processWithCent){ collisionCentFV0A = collision.centFV0A();}
      if constexpr (processWithCent){ collisionCentFT0M = collision.centFT0M();}
      if constexpr (processWithCent){ collisionCentFT0A = collision.centFT0A();}
      if constexpr (processWithCent){ collisionCentFT0C = collision.centFT0C();}
      if constexpr (processWithCent){ collisionCentFDDM = collision.centFDDM();}
      if constexpr (processWithCent){ collisionCentNTPV = collision.centNTPV();}
      
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

        // // // TracksCov
        // ,track.sigmaY() 	    //Covariance matrix
        // ,track.sigmaZ() 	    //Covariance matrix
        // ,track.sigmaSnp() 	  //Covariance matrix
        // ,track.sigmaTgl() 	  //Covariance matrix
        // ,track.sigma1Pt() 	  //Covariance matrix
        // ,track.rhoZY() 	    //Covariance matrix in compressed form
        // ,track.rhoSnpY() 	 	//Covariance matrix in compressed form
        // ,track.rhoSnpZ() 	 	//Covariance matrix in compressed form
        // ,track.rhoTglY() 	 	//Covariance matrix in compressed form
        // ,track.rhoTglZ() 	 	//Covariance matrix in compressed form
        // ,track.rhoTglSnp() 	//Covariance matrix in compressed form
        // ,track.rho1PtY() 	 	//Covariance matrix in compressed form
        // ,track.rho1PtZ() 	 	//Covariance matrix in compressed form
        // ,track.rho1PtSnp() 	//Covariance matrix in compressed form
        // ,track.rho1PtTgl() 	//Covariance matrix in compressed form
        // ,track.cYY() 	  //
        // ,track.cZY() 	  
        // ,track.cZZ() 	  
        // ,track.cSnpY() 	  
        // ,track.cSnpZ() 	  
        // ,track.cSnpSnp()
        // ,track.cTglY() 	  
        // ,track.cTglZ() 	  
        // ,track.cTglSnp()
        // ,track.cTglTgl()
        // ,track.c1PtY() 	  
        // ,track.c1PtZ() 	  
        // ,track.c1PtSnp()
        // ,track.c1PtTgl()
        // ,track.c1Pt21Pt2()

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

        // ,track.tpcNSigmaPi()
        // ,track.tofNSigmaPi()
        // ,track.tpcNSigmaKa()
        // ,track.tofNSigmaKa()
        // ,track.tpcNSigmaPr()
        // ,track.tofNSigmaPr()
        // ,track.tpcNSigmaEl()
        // ,track.tofNSigmaEl()
        // ,track.tpcNSigmaDe()
        // ,track.tofNSigmaDe()

        ,codeSqrtScaling(track.tpcNSigmaPi(), cfgPIDSigma0, cfgPIDSigma1, cfgPIDClampMin, cfgPIDClampMax)
        ,codeSqrtScaling(track.tofNSigmaPi(), cfgPIDSigma0, cfgPIDSigma1, cfgPIDClampMin, cfgPIDClampMax)
        ,codeSqrtScaling(track.tpcNSigmaKa(), cfgPIDSigma0, cfgPIDSigma1, cfgPIDClampMin, cfgPIDClampMax)
        ,codeSqrtScaling(track.tofNSigmaKa(), cfgPIDSigma0, cfgPIDSigma1, cfgPIDClampMin, cfgPIDClampMax)
        ,codeSqrtScaling(track.tpcNSigmaPr(), cfgPIDSigma0, cfgPIDSigma1, cfgPIDClampMin, cfgPIDClampMax)
        ,codeSqrtScaling(track.tofNSigmaPr(), cfgPIDSigma0, cfgPIDSigma1, cfgPIDClampMin, cfgPIDClampMax)
        ,codeSqrtScaling(track.tpcNSigmaEl(), cfgPIDSigma0, cfgPIDSigma1, cfgPIDClampMin, cfgPIDClampMax)
        ,codeSqrtScaling(track.tofNSigmaEl(), cfgPIDSigma0, cfgPIDSigma1, cfgPIDClampMin, cfgPIDClampMax)
        ,codeSqrtScaling(track.tpcNSigmaDe(), cfgPIDSigma0, cfgPIDSigma1, cfgPIDClampMin, cfgPIDClampMax)
        ,codeSqrtScaling(track.tofNSigmaDe(), cfgPIDSigma0, cfgPIDSigma1, cfgPIDClampMin, cfgPIDClampMax)
        
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

        // ,trackMeanOcc.tmoPrimUnfm80() 
        // ,trackMeanOcc.tmoFV0AUnfm80() 
        // ,trackMeanOcc.tmoFV0CUnfm80() 
        // ,trackMeanOcc.tmoFT0AUnfm80() 
        // ,trackMeanOcc.tmoFT0CUnfm80()

        // ,trackMeanOcc.tmoFDDAUnfm80() 
        // ,trackMeanOcc.tmoFDDCUnfm80() 
        
        // ,trackMeanOcc.tmoNTrackITSUnfm80()   
        // ,trackMeanOcc.tmoNTrackTPCUnfm80()   
        // ,trackMeanOcc.tmoNTrackTRDUnfm80()   
        // ,trackMeanOcc.tmoNTrackTOFUnfm80()   
        // ,trackMeanOcc.tmoNTrackSizeUnfm80()   
        // ,trackMeanOcc.tmoNTrackTPCAUnfm80()  
        // ,trackMeanOcc.tmoNTrackTPCCUnfm80()  
        // ,trackMeanOcc.tmoNTrackITSTPCUnfm80()
        // ,trackMeanOcc.tmoNTrackITSTPCAUnfm80()
        // ,trackMeanOcc.tmoNTrackITSTPCCUnfm80()
        
        // ,trackMeanOcc.tmoMultNTracksHasITSUnfm80()
        // ,trackMeanOcc.tmoMultNTracksHasTPCUnfm80()
        // ,trackMeanOcc.tmoMultNTracksHasTOFUnfm80()
        // ,trackMeanOcc.tmoMultNTracksHasTRDUnfm80()
        // ,trackMeanOcc.tmoMultNTracksITSOnlyUnfm80()
        // ,trackMeanOcc.tmoMultNTracksTPCOnlyUnfm80()
        // ,trackMeanOcc.tmoMultNTracksITSTPCUnfm80()
        // ,trackMeanOcc.tmoMultAllTracksTPCOnlyUnfm80()
        
        // ,trackMeanOcc.tmoRobustT0V0PrimUnfm80()   
        // ,trackMeanOcc.tmoRobustFDDT0V0PrimUnfm80()
        // ,trackMeanOcc.tmoRobustNtrackDetUnfm80()  
        // ,trackMeanOcc.tmoRobustMultExtraTableUnfm80()  
        
        // ,trackMeanOcc.twmoPrimUnfm80() 
        // ,trackMeanOcc.twmoFV0AUnfm80() 
        // ,trackMeanOcc.twmoFV0CUnfm80() 
        // ,trackMeanOcc.twmoFT0AUnfm80() 
        // ,trackMeanOcc.twmoFT0CUnfm80() 
        // ,trackMeanOcc.twmoFDDAUnfm80() 
        // ,trackMeanOcc.twmoFDDCUnfm80() 
        // ,trackMeanOcc.twmoNTrackITSUnfm80()   
        // ,trackMeanOcc.twmoNTrackTPCUnfm80()   
        // ,trackMeanOcc.twmoNTrackTRDUnfm80()   
        // ,trackMeanOcc.twmoNTrackTOFUnfm80()   
        // ,trackMeanOcc.twmoNTrackSizeUnfm80()   
        // ,trackMeanOcc.twmoNTrackTPCAUnfm80()  
        // ,trackMeanOcc.twmoNTrackTPCCUnfm80()  
        // ,trackMeanOcc.twmoNTrackITSTPCUnfm80()
        // ,trackMeanOcc.twmoNTrackITSTPCAUnfm80()
        // ,trackMeanOcc.twmoNTrackITSTPCCUnfm80()
        // ,trackMeanOcc.twmoMultNTracksHasITSUnfm80()
        // ,trackMeanOcc.twmoMultNTracksHasTPCUnfm80()
        // ,trackMeanOcc.twmoMultNTracksHasTOFUnfm80()
        // ,trackMeanOcc.twmoMultNTracksHasTRDUnfm80()
        // ,trackMeanOcc.twmoMultNTracksITSOnlyUnfm80()
        // ,trackMeanOcc.twmoMultNTracksTPCOnlyUnfm80()
        // ,trackMeanOcc.twmoMultNTracksITSTPCUnfm80()
        // ,trackMeanOcc.twmoMultAllTracksTPCOnlyUnfm80()
        // ,trackMeanOcc.twmoRobustT0V0PrimUnfm80()   
        // ,trackMeanOcc.twmoRobustFDDT0V0PrimUnfm80()
        // ,trackMeanOcc.twmoRobustNtrackDetUnfm80()  
        // ,trackMeanOcc.twmoRobustMultExtraTableUnfm80()

        // Lossy Compression [START]
        ,codeSqrtScaling(trackMeanOcc.tmoPrimUnfm80(), cfgOccSigma0, cfgOccSigma1, cfgOccClampMin, cfgOccClampMax) 
        ,codeSqrtScaling(trackMeanOcc.tmoFV0AUnfm80(), cfgOccSigma0, cfgOccSigma1, cfgOccClampMin, cfgOccClampMax) 
        ,codeSqrtScaling(trackMeanOcc.tmoFV0CUnfm80(), cfgOccSigma0, cfgOccSigma1, cfgOccClampMin, cfgOccClampMax) 
        ,codeSqrtScaling(trackMeanOcc.tmoFT0AUnfm80(), cfgOccSigma0, cfgOccSigma1, cfgOccClampMin, cfgOccClampMax) 
        ,codeSqrtScaling(trackMeanOcc.tmoFT0CUnfm80(), cfgOccSigma0, cfgOccSigma1, cfgOccClampMin, cfgOccClampMax)

        ,codeSqrtScaling(trackMeanOcc.tmoFDDAUnfm80(), cfgOccSigma0, cfgOccSigma1, cfgOccClampMin, cfgOccClampMax) 
        ,codeSqrtScaling(trackMeanOcc.tmoFDDCUnfm80(), cfgOccSigma0, cfgOccSigma1, cfgOccClampMin, cfgOccClampMax) 
        
        ,codeSqrtScaling(trackMeanOcc.tmoNTrackITSUnfm80(), cfgOccSigma0, cfgOccSigma1, cfgOccClampMin, cfgOccClampMax)   
        ,codeSqrtScaling(trackMeanOcc.tmoNTrackTPCUnfm80(), cfgOccSigma0, cfgOccSigma1, cfgOccClampMin, cfgOccClampMax)   
        ,codeSqrtScaling(trackMeanOcc.tmoNTrackTRDUnfm80(), cfgOccSigma0, cfgOccSigma1, cfgOccClampMin, cfgOccClampMax)   
        ,codeSqrtScaling(trackMeanOcc.tmoNTrackTOFUnfm80(), cfgOccSigma0, cfgOccSigma1, cfgOccClampMin, cfgOccClampMax)   
        ,codeSqrtScaling(trackMeanOcc.tmoNTrackSizeUnfm80(), cfgOccSigma0, cfgOccSigma1, cfgOccClampMin, cfgOccClampMax)   
        ,codeSqrtScaling(trackMeanOcc.tmoNTrackTPCAUnfm80(), cfgOccSigma0, cfgOccSigma1, cfgOccClampMin, cfgOccClampMax)  
        ,codeSqrtScaling(trackMeanOcc.tmoNTrackTPCCUnfm80(), cfgOccSigma0, cfgOccSigma1, cfgOccClampMin, cfgOccClampMax)  
        ,codeSqrtScaling(trackMeanOcc.tmoNTrackITSTPCUnfm80(), cfgOccSigma0, cfgOccSigma1, cfgOccClampMin, cfgOccClampMax)
        ,codeSqrtScaling(trackMeanOcc.tmoNTrackITSTPCAUnfm80(), cfgOccSigma0, cfgOccSigma1, cfgOccClampMin, cfgOccClampMax)
        ,codeSqrtScaling(trackMeanOcc.tmoNTrackITSTPCCUnfm80(), cfgOccSigma0, cfgOccSigma1, cfgOccClampMin, cfgOccClampMax)
        
        ,codeSqrtScaling(trackMeanOcc.tmoMultNTracksHasITSUnfm80(), cfgOccSigma0, cfgOccSigma1, cfgOccClampMin, cfgOccClampMax)
        ,codeSqrtScaling(trackMeanOcc.tmoMultNTracksHasTPCUnfm80(), cfgOccSigma0, cfgOccSigma1, cfgOccClampMin, cfgOccClampMax)
        ,codeSqrtScaling(trackMeanOcc.tmoMultNTracksHasTOFUnfm80(), cfgOccSigma0, cfgOccSigma1, cfgOccClampMin, cfgOccClampMax)
        ,codeSqrtScaling(trackMeanOcc.tmoMultNTracksHasTRDUnfm80(), cfgOccSigma0, cfgOccSigma1, cfgOccClampMin, cfgOccClampMax)
        ,codeSqrtScaling(trackMeanOcc.tmoMultNTracksITSOnlyUnfm80(), cfgOccSigma0, cfgOccSigma1, cfgOccClampMin, cfgOccClampMax)
        ,codeSqrtScaling(trackMeanOcc.tmoMultNTracksTPCOnlyUnfm80(), cfgOccSigma0, cfgOccSigma1, cfgOccClampMin, cfgOccClampMax)
        ,codeSqrtScaling(trackMeanOcc.tmoMultNTracksITSTPCUnfm80(), cfgOccSigma0, cfgOccSigma1, cfgOccClampMin, cfgOccClampMax)
        ,codeSqrtScaling(trackMeanOcc.tmoMultAllTracksTPCOnlyUnfm80(), cfgOccSigma0, cfgOccSigma1, cfgOccClampMin, cfgOccClampMax)
        
        ,codeSqrtScaling(trackMeanOcc.tmoRobustT0V0PrimUnfm80(), cfgOccSigma0, cfgOccSigma1, cfgOccClampMin, cfgOccClampMax)   
        ,codeSqrtScaling(trackMeanOcc.tmoRobustFDDT0V0PrimUnfm80(), cfgOccSigma0, cfgOccSigma1, cfgOccClampMin, cfgOccClampMax)
        ,codeSqrtScaling(trackMeanOcc.tmoRobustNtrackDetUnfm80(), cfgOccSigma0, cfgOccSigma1, cfgOccClampMin, cfgOccClampMax)  
        ,codeSqrtScaling(trackMeanOcc.tmoRobustMultExtraTableUnfm80(), cfgOccSigma0, cfgOccSigma1, cfgOccClampMin, cfgOccClampMax)  
        
        ,codeSqrtScaling(trackMeanOcc.twmoPrimUnfm80(), cfgOccSigma0, cfgOccSigma1, cfgOccClampMin, cfgOccClampMax) 
        ,codeSqrtScaling(trackMeanOcc.twmoFV0AUnfm80(), cfgOccSigma0, cfgOccSigma1, cfgOccClampMin, cfgOccClampMax) 
        ,codeSqrtScaling(trackMeanOcc.twmoFV0CUnfm80(), cfgOccSigma0, cfgOccSigma1, cfgOccClampMin, cfgOccClampMax) 
        ,codeSqrtScaling(trackMeanOcc.twmoFT0AUnfm80(), cfgOccSigma0, cfgOccSigma1, cfgOccClampMin, cfgOccClampMax) 
        ,codeSqrtScaling(trackMeanOcc.twmoFT0CUnfm80(), cfgOccSigma0, cfgOccSigma1, cfgOccClampMin, cfgOccClampMax) 
        ,codeSqrtScaling(trackMeanOcc.twmoFDDAUnfm80(), cfgOccSigma0, cfgOccSigma1, cfgOccClampMin, cfgOccClampMax) 
        ,codeSqrtScaling(trackMeanOcc.twmoFDDCUnfm80(), cfgOccSigma0, cfgOccSigma1, cfgOccClampMin, cfgOccClampMax) 
        ,codeSqrtScaling(trackMeanOcc.twmoNTrackITSUnfm80(), cfgOccSigma0, cfgOccSigma1, cfgOccClampMin, cfgOccClampMax)   
        ,codeSqrtScaling(trackMeanOcc.twmoNTrackTPCUnfm80(), cfgOccSigma0, cfgOccSigma1, cfgOccClampMin, cfgOccClampMax)   
        ,codeSqrtScaling(trackMeanOcc.twmoNTrackTRDUnfm80(), cfgOccSigma0, cfgOccSigma1, cfgOccClampMin, cfgOccClampMax)   
        ,codeSqrtScaling(trackMeanOcc.twmoNTrackTOFUnfm80(), cfgOccSigma0, cfgOccSigma1, cfgOccClampMin, cfgOccClampMax)   
        ,codeSqrtScaling(trackMeanOcc.twmoNTrackSizeUnfm80(), cfgOccSigma0, cfgOccSigma1, cfgOccClampMin, cfgOccClampMax)   
        ,codeSqrtScaling(trackMeanOcc.twmoNTrackTPCAUnfm80(), cfgOccSigma0, cfgOccSigma1, cfgOccClampMin, cfgOccClampMax)  
        ,codeSqrtScaling(trackMeanOcc.twmoNTrackTPCCUnfm80(), cfgOccSigma0, cfgOccSigma1, cfgOccClampMin, cfgOccClampMax)  
        ,codeSqrtScaling(trackMeanOcc.twmoNTrackITSTPCUnfm80(), cfgOccSigma0, cfgOccSigma1, cfgOccClampMin, cfgOccClampMax)
        ,codeSqrtScaling(trackMeanOcc.twmoNTrackITSTPCAUnfm80(), cfgOccSigma0, cfgOccSigma1, cfgOccClampMin, cfgOccClampMax)
        ,codeSqrtScaling(trackMeanOcc.twmoNTrackITSTPCCUnfm80(), cfgOccSigma0, cfgOccSigma1, cfgOccClampMin, cfgOccClampMax)
        ,codeSqrtScaling(trackMeanOcc.twmoMultNTracksHasITSUnfm80(), cfgOccSigma0, cfgOccSigma1, cfgOccClampMin, cfgOccClampMax)
        ,codeSqrtScaling(trackMeanOcc.twmoMultNTracksHasTPCUnfm80(), cfgOccSigma0, cfgOccSigma1, cfgOccClampMin, cfgOccClampMax)
        ,codeSqrtScaling(trackMeanOcc.twmoMultNTracksHasTOFUnfm80(), cfgOccSigma0, cfgOccSigma1, cfgOccClampMin, cfgOccClampMax)
        ,codeSqrtScaling(trackMeanOcc.twmoMultNTracksHasTRDUnfm80(), cfgOccSigma0, cfgOccSigma1, cfgOccClampMin, cfgOccClampMax)
        ,codeSqrtScaling(trackMeanOcc.twmoMultNTracksITSOnlyUnfm80(), cfgOccSigma0, cfgOccSigma1, cfgOccClampMin, cfgOccClampMax)
        ,codeSqrtScaling(trackMeanOcc.twmoMultNTracksTPCOnlyUnfm80(), cfgOccSigma0, cfgOccSigma1, cfgOccClampMin, cfgOccClampMax)
        ,codeSqrtScaling(trackMeanOcc.twmoMultNTracksITSTPCUnfm80(), cfgOccSigma0, cfgOccSigma1, cfgOccClampMin, cfgOccClampMax)
        ,codeSqrtScaling(trackMeanOcc.twmoMultAllTracksTPCOnlyUnfm80(), cfgOccSigma0, cfgOccSigma1, cfgOccClampMin, cfgOccClampMax)
        ,codeSqrtScaling(trackMeanOcc.twmoRobustT0V0PrimUnfm80(), cfgOccSigma0, cfgOccSigma1, cfgOccClampMin, cfgOccClampMax)   
        ,codeSqrtScaling(trackMeanOcc.twmoRobustFDDT0V0PrimUnfm80(), cfgOccSigma0, cfgOccSigma1, cfgOccClampMin, cfgOccClampMax)
        ,codeSqrtScaling(trackMeanOcc.twmoRobustNtrackDetUnfm80(), cfgOccSigma0, cfgOccSigma1, cfgOccClampMin, cfgOccClampMax)  
        ,codeSqrtScaling(trackMeanOcc.twmoRobustMultExtraTableUnfm80(), cfgOccSigma0, cfgOccSigma1, cfgOccClampMin, cfgOccClampMax)
        // Lossy Compression [END]
        
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

        // ,track.itsNSigmaEl()
        // ,track.itsNSigmaMu()
        // ,track.itsNSigmaPi()
        // ,track.itsNSigmaKa()
        // ,track.itsNSigmaPr()
        // ,track.itsNSigmaDe()
        // ,track.itsNSigmaTr()
        // ,track.itsNSigmaHe()
        // ,track.itsNSigmaAl()

        ,codeSqrtScaling(track.itsNSigmaEl(), cfgPIDSigma0, cfgPIDSigma1, cfgPIDClampMin, cfgPIDClampMax)
        ,codeSqrtScaling(track.itsNSigmaMu(), cfgPIDSigma0, cfgPIDSigma1, cfgPIDClampMin, cfgPIDClampMax)
        ,codeSqrtScaling(track.itsNSigmaPi(), cfgPIDSigma0, cfgPIDSigma1, cfgPIDClampMin, cfgPIDClampMax)
        ,codeSqrtScaling(track.itsNSigmaKa(), cfgPIDSigma0, cfgPIDSigma1, cfgPIDClampMin, cfgPIDClampMax)
        ,codeSqrtScaling(track.itsNSigmaPr(), cfgPIDSigma0, cfgPIDSigma1, cfgPIDClampMin, cfgPIDClampMax)
        ,codeSqrtScaling(track.itsNSigmaDe(), cfgPIDSigma0, cfgPIDSigma1, cfgPIDClampMin, cfgPIDClampMax)
        ,codeSqrtScaling(track.itsNSigmaTr(), cfgPIDSigma0, cfgPIDSigma1, cfgPIDClampMin, cfgPIDClampMax)
        ,codeSqrtScaling(track.itsNSigmaHe(), cfgPIDSigma0, cfgPIDSigma1, cfgPIDClampMin, cfgPIDClampMax)
        ,codeSqrtScaling(track.itsNSigmaAl(), cfgPIDSigma0, cfgPIDSigma1, cfgPIDClampMin, cfgPIDClampMax)
 
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

        ,collisionCentFV0A
        ,collisionCentFT0M
        ,collisionCentFT0A
        ,collisionCentFT0C
        ,collisionCentFDDM
        ,collisionCentNTPV

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
    
      if(fillTracksFullInfo){
        DerivedTracks0(      
           track.tpcNSigmaPi()
          ,track.tofNSigmaPi()
          ,track.tpcNSigmaKa()
          ,track.tofNSigmaKa()
          ,track.tpcNSigmaPr()
          ,track.tofNSigmaPr()
          ,track.tpcNSigmaEl()
          ,track.tofNSigmaEl()
          ,track.tpcNSigmaDe()
          ,track.tofNSigmaDe()
        
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

          ,track.itsNSigmaEl()
          ,track.itsNSigmaMu()
          ,track.itsNSigmaPi()
          ,track.itsNSigmaKa()
          ,track.itsNSigmaPr()
          ,track.itsNSigmaDe()
          ,track.itsNSigmaTr()
          ,track.itsNSigmaHe()
          ,track.itsNSigmaAl()
        );
      }
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
        bc = myV0Coll.template bc_as<myBCTable>();
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
        bc = myPhiColl.template bc_as<myBCTable>();
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
/*
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
*/
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
        bc = myD0Coll.template bc_as<myBCTable>();
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
      // DerivedD0s(
      //   d0Counter+gD0Checker
      //   ,d0.globalIndex()
      //   ,bc.globalBC()
      //   ,dfCount
      //   ,d0.collisionId()+gCollisionChecker
      //   ,d0Counter
      //   ,Prong0IndexPosition+gTrackChecker
      //   ,Prong1IndexPosition+gTrackChecker
      //   ,Prong0IndexPosition
      //   ,Prong1IndexPosition
      //   ,d0.prong0Id() 	      
      //   ,d0.prong1Id() 	      
      //   ,d0.collisionId()    
      //   ,d0.posX() 	          
      //   ,d0.posY() 	          
      //   ,d0.posZ() 	          
      //   ,d0.xSecondaryVertex()
      //   ,d0.ySecondaryVertex()
      //   ,d0.zSecondaryVertex()
      //   ,d0.errorDecayLength()
      //   ,d0.errorDecayLengthXY()
      //   ,d0.chi2PCA()
      //   ,d0.rSecondaryVertex() 	          //GI 		? 	
      //   ,d0.decayLength() 	
      //   ,d0.decayLengthXY() 	
      //   ,d0.decayLengthNormalised() 	
      //   ,d0.decayLengthXYNormalised() 	
      //   ,d0.impactParameterNormalised0() 	//GI 		? 	
      //   ,d0.ptProng0() 	
      //   ,d0.pt2Prong0() 	
      //   // ,d0.pVectorProng0()
      //   ,d0.impactParameterNormalised1() 	//GI 		? 	
      //   ,d0.ptProng1() 	
      //   ,d0.pt2Prong1() 	
      //   // ,d0.pVectorProng1()
      //   ,d0.pxProng0() 	
      //   ,d0.pyProng0() 	
      //   ,d0.pzProng0() 	
      //   ,d0.pxProng1() 	
      //   ,d0.pyProng1() 	
      //   ,d0.pzProng1() 	
      //   ,d0.impactParameter0() 	
      //   ,d0.impactParameter1() 	
      //   ,d0.errorImpactParameter0() 	
      //   ,d0.errorImpactParameter1() 	
      //   ,d0.impactParameterZ0() 	
      //   ,d0.impactParameterZ1() 	
      //   ,d0.errorImpactParameterZ0() 	
      //   ,d0.errorImpactParameterZ1() 	
      //   ,d0.nProngsContributorsPV()
      //   ,d0.hfflag()
      //   // ,d0.m() 	                          //    GI 		? 	
      //   // ,d0.m2() 	
      //   ,d0.impactParameterProduct() 	
      //   // ,d0.cosThetaStar() 	
      //   ,d0.impactParameterProngSqSum() 	
      //   ,d0.pt() 	                                //GI 		? 	
      //   ,d0.pt2() 	
      //   ,d0.p() 	
      //   ,d0.p2() 	
      //   // // // ,d0.pVector()
      //   ,d0.cpa() 	
      //   ,d0.cpaXY() 	
      //   // // // ,d0.ct() 	
      //   ,d0.impactParameterXY() 	
      //   ,d0.maxNormalisedDeltaIP()
      //   ,d0.px()
      //   ,d0.py()
      //   ,d0.pz()
      //   ,d0.eta() 	
      //   ,d0.phi() 	
      //   ,d0.y(1.86484) 	//Mass in GeV/c2
      //   ,d0.e(1.86484)  //Mass in GeV/c2

      //   ,time
      //   ,TFidThis
      //   ,bcInTF

      //   ,weightDSPtD0[indexPosition]
      //   ,weightDSPtD0Bar[indexPosition]
      //   ,vecTriggerMaskD0[indexPosition]

      //   ,hEntriesD0Pt             
      //   ,hEntriesD0BarPt         

      //   ,histPDFPol1Evaluation(recoEvent.get<TH1>(HIST("Event/D0/h_Pt" ))  , d0.pt() , 2)
      //   ,histPDFPol1Evaluation(recoEvent.get<TH1>(HIST("Event/D0Bar/h_Pt")), d0.pt() , 2)
      // );

      d0Counter++;
    }
    PrintTime(StartD0Table, Form("DEBUG :: df_%ld :: D0 Table Time :: ",dfCount));

    //
    //Phi reconstruction table
    for(const auto& collision: collisions){
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

          // Do Kaon selection
          // We need to select the tracks carefully to have a Phi peak
          // We add global track selection to see if we get it
          // Otherwise we need to put some cut on collision as well for this (We'll do this later)


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

    /////////////////////////////////////////////////////////////////////////////////////////
    auto StartZ0Candidate = std::chrono::high_resolution_clock::now();
    bool isZ0Candidate = false;
    uint8_t  decayBit = 0;
    bool isZ0CandidateFromEl =false;
    bool isZ0CandidateFromMu =false;
    TLorentzVector dauEl1, dauEl2, motherFromEl;
    TLorentzVector dauMu1, dauMu2, motherFromMu;

    int z0Counter = 0;
    lastCollId = -999;
    int posDuaIndexPosition = -1;
    int negDuaIndexPosition = -1;
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
          
          isZ0CandidateFromEl =false;
          isZ0CandidateFromMu =false;
          decayBit = 0;
          // Do Electron and Muon selection
          if(selElectron(posTrack) && selElectron(negTrack)){isZ0CandidateFromEl = true;}
          if(selMuon(posTrack) && selMuon(negTrack) ){isZ0CandidateFromMu = true;}

          if(!isZ0CandidateFromEl && !isZ0CandidateFromMu){ continue; }

          if(isZ0CandidateFromEl){
            //LorentzVecotrApproach for mass
            dauEl1.SetXYZM(posTrack.px(), posTrack.py(), posTrack.pz(), MassElectron);
            dauEl2.SetXYZM(negTrack.px(), negTrack.py(), negTrack.pz(), MassElectron);
            motherFromEl = dauEl1 + dauEl2;

            // Mll​=\sqrt(2 * p_{T_1} *​ p_{T_2}​(cosh(Δη)−cos(Δϕ))
            float mass = motherFromEl.M();
            if( 80 < mass && mass < 100 ) {
              BITSET(decayBit, 0);
            }
          }

          if(isZ0CandidateFromMu){
            //LorentzVecotrApproach for mass
            dauMu1.SetXYZM(posTrack.px(), posTrack.py(), posTrack.pz(), MassMuonMinus);
            dauMu2.SetXYZM(negTrack.px(), negTrack.py(), negTrack.pz(), MassMuonMinus);
            motherFromMu = dauMu1 + dauMu2;
            float mass = motherFromMu.M();
            if( 80 < mass && mass < 100 ) {
              BITSET(decayBit, 1);
            }
          }          

          if (lastCollId != collision.globalIndex()) {
            lastCollId = collision.globalIndex();
            bc = collision.template bc_as<myBCTable>();
            GetTimingInfo(bc, lastRun, nBCsPerTF, bcSOR, time, TFidThis, bcInTF);
          }

          posDuaIndexPosition = FindTrackIdInList(posTrack.globalIndex(), trackGlobalIndexList, DrTrackPositionIndexList, trackGlobalIndexList.size());
          negDuaIndexPosition = FindTrackIdInList(negTrack.globalIndex(), trackGlobalIndexList, DrTrackPositionIndexList, trackGlobalIndexList.size());

          
          if(BITCHECK(decayBit, 0) != 0){
            DerivedZ0s(
              z0Counter+gZ0Checker
              ,z0Counter
              ,bc.globalBC()
              ,dfCount
              ,collision.globalIndex()+gCollisionChecker
              ,z0Counter+gZ0Checker

              ,posDuaIndexPosition+gTrackChecker //global Id in merged DFs in derived Tracks
              ,negDuaIndexPosition+gTrackChecker //global Id in merged DFs in derived Tracks

              ,posDuaIndexPosition //global Id in single DF in derived Tracks
              ,negDuaIndexPosition //global Id in single DF in derived Tracks

              ,posTrack.globalIndex() //original Id in single DF in original Tracks
              ,negTrack.globalIndex() //original Id in single DF in original Tracks

              ,decayBit

              ,motherFromEl.Px()    // phiPx
              ,motherFromEl.Py()    // phiPy
              ,motherFromEl.Pz()    // phiPz
              ,motherFromEl.Pt()
              ,motherFromEl.Eta()
              ,motherFromEl.Phi()
              ,motherFromEl.Rapidity()
              ,motherFromEl.M()
              ,motherFromEl.E()

              ,posTrack.px()
              ,posTrack.py()
              ,posTrack.pz()
              ,posTrack.pt()
              ,posTrack.eta()

              ,negTrack.px()
              ,negTrack.py()
              ,negTrack.pz()
              ,negTrack.pt()
              ,negTrack.eta()

              ,time
              ,TFidThis
              ,bcInTF
            );
            z0Counter++;
          }

          if(BITCHECK(decayBit, 1) != 0){
            DerivedZ0s(
              z0Counter+gZ0Checker
              ,z0Counter
              ,bc.globalBC()
              ,dfCount
              ,collision.globalIndex()+gCollisionChecker
              ,z0Counter+gZ0Checker

              ,posDuaIndexPosition+gTrackChecker //global Id in merged DFs in derived Tracks
              ,negDuaIndexPosition+gTrackChecker //global Id in merged DFs in derived Tracks

              ,posDuaIndexPosition //global Id in single DF in derived Tracks
              ,negDuaIndexPosition //global Id in single DF in derived Tracks

              ,posTrack.globalIndex() //original Id in single DF in original Tracks
              ,negTrack.globalIndex() //original Id in single DF in original Tracks

              ,decayBit

              ,motherFromMu.Px()    // phiPx
              ,motherFromMu.Py()    // phiPy
              ,motherFromMu.Pz()    // phiPz
              ,motherFromMu.Pt()
              ,motherFromMu.Eta()
              ,motherFromMu.Phi()
              ,motherFromMu.Rapidity()
              ,motherFromMu.M()
              ,motherFromMu.E()

              ,posTrack.px()
              ,posTrack.py()
              ,posTrack.pz()
              ,posTrack.pt()
              ,posTrack.eta()

              ,negTrack.px()
              ,negTrack.py()
              ,negTrack.pz()
              ,negTrack.pt()
              ,negTrack.eta()

              ,time
              ,TFidThis
              ,bcInTF
            );
            z0Counter++;
          }
        }
      }
    }

    PrintTime(StartZ0Candidate, Form("DEBUG :: df_%ld :: Z0 Table Time :: ",dfCount));

    auto StartCosmicMuon = std::chrono::high_resolution_clock::now();

    int cosmicPairCounter = 0;
    int upperTrkIndexPosition ; 
    int lowerTrkIndexPosition ; 

    float fSumPtPair;    
    float fSumQPtPair;   
    float fSumTglPair;   
    float fSumDcaXYPair;   
    float fSumAlphaPair;

    float fDiffPtPair;   
    float fDiffQPtPair;  
    float fDiffTglPair;  
    float fDiffDcaXYPair;  
    float fDiffAlphaPair;

    auto time1  = time;
    auto TFidThis1  = TFidThis;
    auto bcInTF1  = bcInTF;

    auto time2  = time;
    auto TFidThis2  = TFidThis;
    auto bcInTF2  = bcInTF;

    int lastCollId1 = -999;
    int lastCollId2 = -999;

    float upperTrkTpcTime0           = -999 ; float lowerTrkTpcTime0           = -999 ;
    float upperTrkTpcdcaR            = -999 ; float lowerTrkTpcdcaR            = -999 ;
    float upperTrkTpcdcaZ            = -999 ; float lowerTrkTpcdcaZ            = -999 ;
    float upperTrkTpcClusterByteMask = -999 ; float lowerTrkTpcClusterByteMask = -999 ;
    float upperTrkTpcdEdxMax0R       = -999 ; float lowerTrkTpcdEdxMax0R       = -999 ;
    float upperTrkTpcdEdxMax1R       = -999 ; float lowerTrkTpcdEdxMax1R       = -999 ;
    float upperTrkTpcdEdxMax2R       = -999 ; float lowerTrkTpcdEdxMax2R       = -999 ;
    float upperTrkTpcdEdxMax3R       = -999 ; float lowerTrkTpcdEdxMax3R       = -999 ;
    float upperTrkTpcdEdxTot0R       = -999 ; float lowerTrkTpcdEdxTot0R       = -999 ;
    float upperTrkTpcdEdxTot1R       = -999 ; float lowerTrkTpcdEdxTot1R       = -999 ;
    float upperTrkTpcdEdxTot2R       = -999 ; float lowerTrkTpcdEdxTot2R       = -999 ;
    float upperTrkTpcdEdxTot3R       = -999 ; float lowerTrkTpcdEdxTot3R       = -999 ;

    float upperTrkDeltaRefContParamY    = -999; float lowerTrkDeltaRefContParamY    = -999;
    float upperTrkDeltaRefContParamZ    = -999; float lowerTrkDeltaRefContParamZ    = -999;
    float upperTrkDeltaRefContParamSnp  = -999; float lowerTrkDeltaRefContParamSnp  = -999;
    float upperTrkDeltaRefContParamTgl  = -999; float lowerTrkDeltaRefContParamTgl  = -999;
    float upperTrkDeltaRefContParamQ2Pt = -999; float lowerTrkDeltaRefContParamQ2Pt = -999;
    float upperTrkDeltaRefGloParamY     = -999; float lowerTrkDeltaRefGloParamY     = -999;
    float upperTrkDeltaRefGloParamZ     = -999; float lowerTrkDeltaRefGloParamZ     = -999;
    float upperTrkDeltaRefGloParamSnp   = -999; float lowerTrkDeltaRefGloParamSnp   = -999;
    float upperTrkDeltaRefGloParamTgl   = -999; float lowerTrkDeltaRefGloParamTgl   = -999;
    float upperTrkDeltaRefGloParamQ2Pt  = -999; float lowerTrkDeltaRefGloParamQ2Pt  = -999;
    float upperTrkDeltaTOFdX            = -999; float lowerTrkDeltaTOFdX            = -999;
    float upperTrkDeltaTOFdZ            = -999; float lowerTrkDeltaTOFdZ            = -999;

    uint8_t upperTrkDetectorMask = 0;
    uint8_t lowerTrkDetectorMask = 0;

    int trackCounterDebugger = 0;
    std::vector<int64_t> upperTrackTime;
    std::vector<int64_t> upperTrackTFidThis;
    std::vector<int64_t> upperTrackBcInTF;

    std::vector<int64_t> lowerTrackTime;
    std::vector<int64_t> lowerTrackTFidThis;
    std::vector<int64_t> lowerTrackBcInTF;

    auto uColl = collisions.begin();
    auto lColl = collisions.begin();

    auto pairRejectionHist = recoEvent.get<TH1>(HIST("CosmicMuon/PairRejectionTrigger"));
    
    for(const auto& upperTrk : tracks){
      // if(upperTrk.y() < 0 ) {continue;}
      // if ( !(0 <= upperTrk.alpha() && upperTrk.alpha() < TMath::Pi())) {continue;}
      if(upperTrk.pt() < cfgCM.csmMuThrPt ) {continue;}
      if(upperTrk.dcaXY() < cfgCM.csmMuThrDCAMin ) {continue;}
      if(upperTrk.dcaXY() > cfgCM.csmMuThrDCAMax ) {continue;}

      if(upperTrk.collisionId()< 0 ) {
        upperTrackTime.clear();
        upperTrackTFidThis.clear();
        upperTrackBcInTF.clear();
        auto bcs = upperTrk.template ambgTrack_as<A>().template bc_as<B>();
        for (const auto& bc: bcs){
          GetTimingInfo(bc, lastRun, nBCsPerTF, bcSOR, time, TFidThis, bcInTF);
          upperTrackTime.push_back(time);
          upperTrackTFidThis.push_back(TFidThis);
          upperTrackBcInTF.push_back(bcInTF);
        }
      }
      else {
        uColl = upperTrk.template collision_as<C>();
        if (lastCollId1 != uColl.globalIndex()) {
          upperTrackTime.clear();
          upperTrackTFidThis.clear();
          upperTrackBcInTF.clear();
          lastCollId1 = uColl.globalIndex();
          bc = uColl.template bc_as<myBCTable>();
          GetTimingInfo(bc, lastRun, nBCsPerTF, bcSOR, time, TFidThis, bcInTF);
          upperTrackTime.push_back(time);
          upperTrackTFidThis.push_back(TFidThis);
          upperTrackBcInTF.push_back(bcInTF);
        }
      }

      upperTrkTpcTime0           = -999 ; lowerTrkTpcTime0           = -999 ;
      upperTrkTpcdcaR            = -999 ; lowerTrkTpcdcaR            = -999 ;
      upperTrkTpcdcaZ            = -999 ; lowerTrkTpcdcaZ            = -999 ;
      upperTrkTpcClusterByteMask = -999 ; lowerTrkTpcClusterByteMask = -999 ;
      upperTrkTpcdEdxMax0R       = -999 ; lowerTrkTpcdEdxMax0R       = -999 ;
      upperTrkTpcdEdxMax1R       = -999 ; lowerTrkTpcdEdxMax1R       = -999 ;
      upperTrkTpcdEdxMax2R       = -999 ; lowerTrkTpcdEdxMax2R       = -999 ;
      upperTrkTpcdEdxMax3R       = -999 ; lowerTrkTpcdEdxMax3R       = -999 ;
      upperTrkTpcdEdxTot0R       = -999 ; lowerTrkTpcdEdxTot0R       = -999 ;
      upperTrkTpcdEdxTot1R       = -999 ; lowerTrkTpcdEdxTot1R       = -999 ;
      upperTrkTpcdEdxTot2R       = -999 ; lowerTrkTpcdEdxTot2R       = -999 ;
      upperTrkTpcdEdxTot3R       = -999 ; lowerTrkTpcdEdxTot3R       = -999 ;

      upperTrkDeltaRefContParamY    = -999; lowerTrkDeltaRefContParamY    = -999;
      upperTrkDeltaRefContParamZ    = -999; lowerTrkDeltaRefContParamZ    = -999;
      upperTrkDeltaRefContParamSnp  = -999; lowerTrkDeltaRefContParamSnp  = -999;
      upperTrkDeltaRefContParamTgl  = -999; lowerTrkDeltaRefContParamTgl  = -999;
      upperTrkDeltaRefContParamQ2Pt = -999; lowerTrkDeltaRefContParamQ2Pt = -999;
      upperTrkDeltaRefGloParamY     = -999; lowerTrkDeltaRefGloParamY     = -999;
      upperTrkDeltaRefGloParamZ     = -999; lowerTrkDeltaRefGloParamZ     = -999;
      upperTrkDeltaRefGloParamSnp   = -999; lowerTrkDeltaRefGloParamSnp   = -999;
      upperTrkDeltaRefGloParamTgl   = -999; lowerTrkDeltaRefGloParamTgl   = -999;
      upperTrkDeltaRefGloParamQ2Pt  = -999; lowerTrkDeltaRefGloParamQ2Pt  = -999;
      upperTrkDeltaTOFdX            = -999; lowerTrkDeltaTOFdX            = -999;
      upperTrkDeltaTOFdZ            = -999; lowerTrkDeltaTOFdZ            = -999;

      if(upperTrk.trackQAId() != -1){
        auto trackQA_From_Tracks = upperTrk.template trackQA_as<MyTracksQA>();//Dereferencing aod::Tracks ==> aod::TracksQAVersion
        if(upperTrk.globalIndex() != trackQA_From_Tracks.trackId()){
          LOG(info)<<"DEBUG :: ERROR ERROR ERROR :: upperTrk.globalIndex() != trackQA_From_Tracks.trackId() :: "<<upperTrk.globalIndex()<<" !=  "<<trackQA_From_Tracks.trackId(); //this is also failing 
        }

        upperTrkTpcTime0           = trackQA_From_Tracks.tpcTime0          ();
        upperTrkTpcdcaR            = trackQA_From_Tracks.tpcdcaR           ();
        upperTrkTpcdcaZ            = trackQA_From_Tracks.tpcdcaZ           ();
        upperTrkTpcClusterByteMask = trackQA_From_Tracks.tpcClusterByteMask();
        upperTrkTpcdEdxMax0R       = trackQA_From_Tracks.tpcdEdxMax0R      ();
        upperTrkTpcdEdxMax1R       = trackQA_From_Tracks.tpcdEdxMax1R      ();
        upperTrkTpcdEdxMax2R       = trackQA_From_Tracks.tpcdEdxMax2R      ();
        upperTrkTpcdEdxMax3R       = trackQA_From_Tracks.tpcdEdxMax3R      ();
        upperTrkTpcdEdxTot0R       = trackQA_From_Tracks.tpcdEdxTot0R      ();
        upperTrkTpcdEdxTot1R       = trackQA_From_Tracks.tpcdEdxTot1R      ();
        upperTrkTpcdEdxTot2R       = trackQA_From_Tracks.tpcdEdxTot2R      ();
        upperTrkTpcdEdxTot3R       = trackQA_From_Tracks.tpcdEdxTot3R      ();

        upperTrkDeltaRefContParamY    = trackQA_From_Tracks.deltaRefContParamY(); 
        upperTrkDeltaRefContParamZ    = trackQA_From_Tracks.deltaRefITSParamZ();
        upperTrkDeltaRefContParamSnp  = trackQA_From_Tracks.deltaRefContParamSnp();
        upperTrkDeltaRefContParamTgl  = trackQA_From_Tracks.deltaRefContParamTgl();
        upperTrkDeltaRefContParamQ2Pt = trackQA_From_Tracks.deltaRefContParamQ2Pt();
        upperTrkDeltaRefGloParamY     = trackQA_From_Tracks.deltaRefGloParamY();
        upperTrkDeltaRefGloParamZ     = trackQA_From_Tracks.deltaRefGloParamZ();
        upperTrkDeltaRefGloParamSnp   = trackQA_From_Tracks.deltaRefGloParamSnp();
        upperTrkDeltaRefGloParamTgl   = trackQA_From_Tracks.deltaRefGloParamTgl();
        upperTrkDeltaRefGloParamQ2Pt  = trackQA_From_Tracks.deltaRefGloParamQ2Pt();
        upperTrkDeltaTOFdX            = trackQA_From_Tracks.deltaTOFdX();
        upperTrkDeltaTOFdZ            = trackQA_From_Tracks.deltaTOFdZ();
      }

      for(const auto& lowerTrk : tracks){
        // if (lowerTrk.y() > 0 ) {continue;}
        // if ( !(-TMath::Pi() <= upperTrk.alpha() && upperTrk.alpha() < 0 )) {continue;}
        if(lowerTrk.pt() < cfgCM.csmMuThrPt ) {continue;}
        if(lowerTrk.dcaXY() < cfgCM.csmMuThrDCAMin ) {continue;}
        if(lowerTrk.dcaXY() > cfgCM.csmMuThrDCAMax ) {continue;}

        if( lowerTrk.globalIndex() <= upperTrk.globalIndex()) {continue;} //Strictly upper index policy for fast grouping
        if((upperTrk.signed1Pt() * lowerTrk.signed1Pt()) < 0) {continue;} //pair track should be of same sign

        fSumPtPair     = std::abs(upperTrk.pt()        + lowerTrk.pt());
        fSumQPtPair    = std::abs(upperTrk.signed1Pt() + lowerTrk.signed1Pt());
        fSumTglPair    = std::abs(upperTrk.tgl()       + lowerTrk.tgl());
        // fSumDcaXYPair  = std::abs(upperTrk.dcaXY())    + std::abs(lowerTrk.dcaXY());
        if( upperTrk.dcaXY() * lowerTrk.dcaXY() < 0 ){
          fSumDcaXYPair  = std::abs(upperTrk.dcaXY())    + std::abs(lowerTrk.dcaXY());
        }else {
          fSumDcaXYPair  = std::abs(upperTrk.dcaXY() + lowerTrk.dcaXY());
        }
        fSumAlphaPair  = std::abs(upperTrk.alpha()     + lowerTrk.alpha());
        fDiffPtPair    = std::abs(upperTrk.pt()        - lowerTrk.pt());
        fDiffQPtPair   = std::abs(upperTrk.signed1Pt() - lowerTrk.signed1Pt());
        fDiffTglPair   = std::abs(upperTrk.tgl()       - lowerTrk.tgl());
        fDiffDcaXYPair = std::abs(upperTrk.dcaXY()     - lowerTrk.dcaXY());
        fDiffAlphaPair = std::abs(deltaAngleCalculator(upperTrk.alpha(),lowerTrk.alpha()));

        recoEvent.fill(HIST("CosmicMuon/PreSel/sumPtPair"    ), fSumPtPair     );
        recoEvent.fill(HIST("CosmicMuon/PreSel/sumQPtPair"   ), fSumQPtPair    );
        recoEvent.fill(HIST("CosmicMuon/PreSel/sumTglPair"   ), fSumTglPair    );
        recoEvent.fill(HIST("CosmicMuon/PreSel/sumDcaXYPair" ), fSumDcaXYPair  );
        recoEvent.fill(HIST("CosmicMuon/PreSel/sumAlphaPair" ), fSumAlphaPair  );

        recoEvent.fill(HIST("CosmicMuon/PreSel/diffPtPair"   ), fDiffPtPair    );
        recoEvent.fill(HIST("CosmicMuon/PreSel/diffQPtPair"  ), fDiffQPtPair   );
        recoEvent.fill(HIST("CosmicMuon/PreSel/diffTglPair"  ), fDiffTglPair   );
        recoEvent.fill(HIST("CosmicMuon/PreSel/diffDcaXYPair"), fDiffDcaXYPair );
        recoEvent.fill(HIST("CosmicMuon/PreSel/diffAlphaPair"), fDiffAlphaPair );

        if ( fSumPtPair < cfgCM.csmMuSumPtPair                                           ) { pairRejectionHist->Fill(kFailSumPtPair    ); continue;}
        if ( fSumQPtPair > cfgCM.csmMuSumQPtPair                                         ) { pairRejectionHist->Fill(kFailSumQPtPair   ); }//continue;}
        if ( fSumTglPair > cfgCM.csmMuSumTglPair                                         ) { pairRejectionHist->Fill(kFailSumTglPair   ); }//continue;}
        if (cfgCM.csmMuCheckSumDcaXY && fSumDcaXYPair > cfgCM.csmMuSumDcaXY              ) { pairRejectionHist->Fill(kFailSumDcaXY     ); }//continue;}
        if (cfgCM.csmMuCheckDiffDcaXY && fDiffDcaXYPair > cfgCM.csmMuDiffDcaXY           ) { pairRejectionHist->Fill(kFailDiffDcaXY    ); }//continue;}
        if (cfgCM.csmMuCheckDiffAlphaPair && fDiffAlphaPair > cfgCM.csmMuDiffAlphaPair   ) { pairRejectionHist->Fill(kFailDiffAlphaPair); }//continue;}

        recoEvent.fill(HIST("CosmicMuon/PostSel/sumPtPair"    ), fSumPtPair     );
        recoEvent.fill(HIST("CosmicMuon/PostSel/sumQPtPair"   ), fSumQPtPair    );
        recoEvent.fill(HIST("CosmicMuon/PostSel/sumTglPair"   ), fSumTglPair    );
        recoEvent.fill(HIST("CosmicMuon/PostSel/sumDcaXYPair" ), fSumDcaXYPair  );
        recoEvent.fill(HIST("CosmicMuon/PostSel/sumAlphaPair" ), fSumAlphaPair  );

        recoEvent.fill(HIST("CosmicMuon/PostSel/diffPtPair"   ), fDiffPtPair    );
        recoEvent.fill(HIST("CosmicMuon/PostSel/diffQPtPair"  ), fDiffQPtPair   );
        recoEvent.fill(HIST("CosmicMuon/PostSel/diffTglPair"  ), fDiffTglPair   );
        recoEvent.fill(HIST("CosmicMuon/PostSel/diffDcaXYPair"), fDiffDcaXYPair );
        recoEvent.fill(HIST("CosmicMuon/PostSel/diffAlphaPair"), fDiffAlphaPair );

        if(lowerTrk.collisionId()< 0 ) {
          lowerTrackTime.clear();
          lowerTrackTFidThis.clear();
          lowerTrackBcInTF.clear();
          auto bcs = lowerTrk.template ambgTrack_as<A>().template bc_as<B>();
          for (const auto& bc: bcs){
            GetTimingInfo(bc, lastRun, nBCsPerTF, bcSOR, time, TFidThis, bcInTF);
            lowerTrackTime.push_back(time);
            lowerTrackTFidThis.push_back(TFidThis);
            lowerTrackBcInTF.push_back(bcInTF);
          }
        }
        else {
          lColl = lowerTrk.template collision_as<C>();
          if (lastCollId1 != uColl.globalIndex()) {
            lowerTrackTime.clear();
            lowerTrackTFidThis.clear();
            lowerTrackBcInTF.clear();
            lastCollId1 = uColl.globalIndex();
            bc = uColl.template bc_as<myBCTable>();
            GetTimingInfo(bc, lastRun, nBCsPerTF, bcSOR, time, TFidThis, bcInTF);
            lowerTrackTime.push_back(time);
            lowerTrackTFidThis.push_back(TFidThis);
            lowerTrackBcInTF.push_back(bcInTF);
          }
        }

        lowerTrkTpcTime0           = -999 ; lowerTrkTpcTime0           = -999 ;
        lowerTrkTpcdcaR            = -999 ; lowerTrkTpcdcaR            = -999 ;
        lowerTrkTpcdcaZ            = -999 ; lowerTrkTpcdcaZ            = -999 ;
        lowerTrkTpcClusterByteMask = -999 ; lowerTrkTpcClusterByteMask = -999 ;
        lowerTrkTpcdEdxMax0R       = -999 ; lowerTrkTpcdEdxMax0R       = -999 ;
        lowerTrkTpcdEdxMax1R       = -999 ; lowerTrkTpcdEdxMax1R       = -999 ;
        lowerTrkTpcdEdxMax2R       = -999 ; lowerTrkTpcdEdxMax2R       = -999 ;
        lowerTrkTpcdEdxMax3R       = -999 ; lowerTrkTpcdEdxMax3R       = -999 ;
        lowerTrkTpcdEdxTot0R       = -999 ; lowerTrkTpcdEdxTot0R       = -999 ;
        lowerTrkTpcdEdxTot1R       = -999 ; lowerTrkTpcdEdxTot1R       = -999 ;
        lowerTrkTpcdEdxTot2R       = -999 ; lowerTrkTpcdEdxTot2R       = -999 ;
        lowerTrkTpcdEdxTot3R       = -999 ; lowerTrkTpcdEdxTot3R       = -999 ;

        lowerTrkDeltaRefContParamY    = -999; lowerTrkDeltaRefContParamY    = -999;
        lowerTrkDeltaRefContParamZ    = -999; lowerTrkDeltaRefContParamZ    = -999;
        lowerTrkDeltaRefContParamSnp  = -999; lowerTrkDeltaRefContParamSnp  = -999;
        lowerTrkDeltaRefContParamTgl  = -999; lowerTrkDeltaRefContParamTgl  = -999;
        lowerTrkDeltaRefContParamQ2Pt = -999; lowerTrkDeltaRefContParamQ2Pt = -999;
        lowerTrkDeltaRefGloParamY     = -999; lowerTrkDeltaRefGloParamY     = -999;
        lowerTrkDeltaRefGloParamZ     = -999; lowerTrkDeltaRefGloParamZ     = -999;
        lowerTrkDeltaRefGloParamSnp   = -999; lowerTrkDeltaRefGloParamSnp   = -999;
        lowerTrkDeltaRefGloParamTgl   = -999; lowerTrkDeltaRefGloParamTgl   = -999;
        lowerTrkDeltaRefGloParamQ2Pt  = -999; lowerTrkDeltaRefGloParamQ2Pt  = -999;
        lowerTrkDeltaTOFdX            = -999; lowerTrkDeltaTOFdX            = -999;
        lowerTrkDeltaTOFdZ            = -999; lowerTrkDeltaTOFdZ            = -999;

        if(lowerTrk.trackQAId() != -1){
          auto trackQA_From_Tracks = lowerTrk.template trackQA_as<MyTracksQA>();//Dereferencing aod::Tracks ==> aod::TracksQAVersion
          if(lowerTrk.globalIndex() != trackQA_From_Tracks.trackId()){
            LOG(info)<<"DEBUG :: ERROR ERROR ERROR :: lowerTrk.globalIndex() != trackQA_From_Tracks.trackId() :: "<<lowerTrk.globalIndex()<<" !=  "<<trackQA_From_Tracks.trackId(); //this is also failing 
          }

          lowerTrkTpcTime0           = trackQA_From_Tracks.tpcTime0          ();
          lowerTrkTpcdcaR            = trackQA_From_Tracks.tpcdcaR           ();
          lowerTrkTpcdcaZ            = trackQA_From_Tracks.tpcdcaZ           ();
          lowerTrkTpcClusterByteMask = trackQA_From_Tracks.tpcClusterByteMask();
          lowerTrkTpcdEdxMax0R       = trackQA_From_Tracks.tpcdEdxMax0R      ();
          lowerTrkTpcdEdxMax1R       = trackQA_From_Tracks.tpcdEdxMax1R      ();
          lowerTrkTpcdEdxMax2R       = trackQA_From_Tracks.tpcdEdxMax2R      ();
          lowerTrkTpcdEdxMax3R       = trackQA_From_Tracks.tpcdEdxMax3R      ();
          lowerTrkTpcdEdxTot0R       = trackQA_From_Tracks.tpcdEdxTot0R      ();
          lowerTrkTpcdEdxTot1R       = trackQA_From_Tracks.tpcdEdxTot1R      ();
          lowerTrkTpcdEdxTot2R       = trackQA_From_Tracks.tpcdEdxTot2R      ();
          lowerTrkTpcdEdxTot3R       = trackQA_From_Tracks.tpcdEdxTot3R      ();

          lowerTrkDeltaRefContParamY    = trackQA_From_Tracks.deltaRefContParamY(); 
          lowerTrkDeltaRefContParamZ    = trackQA_From_Tracks.deltaRefITSParamZ();
          lowerTrkDeltaRefContParamSnp  = trackQA_From_Tracks.deltaRefContParamSnp();
          lowerTrkDeltaRefContParamTgl  = trackQA_From_Tracks.deltaRefContParamTgl();
          lowerTrkDeltaRefContParamQ2Pt = trackQA_From_Tracks.deltaRefContParamQ2Pt();
          lowerTrkDeltaRefGloParamY     = trackQA_From_Tracks.deltaRefGloParamY();
          lowerTrkDeltaRefGloParamZ     = trackQA_From_Tracks.deltaRefGloParamZ();
          lowerTrkDeltaRefGloParamSnp   = trackQA_From_Tracks.deltaRefGloParamSnp();
          lowerTrkDeltaRefGloParamTgl   = trackQA_From_Tracks.deltaRefGloParamTgl();
          lowerTrkDeltaRefGloParamQ2Pt  = trackQA_From_Tracks.deltaRefGloParamQ2Pt();
          lowerTrkDeltaTOFdX            = trackQA_From_Tracks.deltaTOFdX();
          lowerTrkDeltaTOFdZ            = trackQA_From_Tracks.deltaTOFdZ();
        }

        upperTrkIndexPosition = FindTrackIdInList(upperTrk.globalIndex(), trackGlobalIndexList, DrTrackPositionIndexList, trackGlobalIndexList.size());
        lowerTrkIndexPosition = FindTrackIdInList(lowerTrk.globalIndex(), trackGlobalIndexList, DrTrackPositionIndexList, trackGlobalIndexList.size());

        upperTrkDetectorMask = (1<<0)*upperTrk.hasITS() + (1<<1)*upperTrk.hasTPC() + (1<<2)*upperTrk.hasTRD() + (1<<3)*upperTrk.hasTOF();
        lowerTrkDetectorMask = (1<<0)*lowerTrk.hasITS() + (1<<1)*lowerTrk.hasTPC() + (1<<2)*lowerTrk.hasTRD() + (1<<3)*lowerTrk.hasTOF();

        DerivedCosmicPairs(
          cosmicPairCounter+gCosmicPairChecker
          ,cosmicPairCounter
          ,dfCount

          ,uColl.template bc_as<B>().globalBC()
          ,uColl.globalIndex()+gCollisionChecker
          ,upperTrkIndexPosition+gTrackChecker //global Id in merged DFs in derived Tracks
          ,upperTrkIndexPosition //global Id in single DF in derived Tracks
          ,upperTrk.globalIndex() //original Id in single DF in original Tracks

          ,lColl.template bc_as<B>().globalBC()
          ,lColl.globalIndex()+gCollisionChecker
          ,lowerTrkIndexPosition+gTrackChecker //global Id in merged DFs in derived Tracks
          ,lowerTrkIndexPosition //global Id in single DF in derived Tracks
          ,lowerTrk.globalIndex() //original Id in single DF in original Tracks


          ,upperTrk.x() 	                // 
          ,upperTrk.alpha() 	            // 
          ,upperTrk.y() 	                // 
          ,upperTrk.z() 	                // 
          ,upperTrk.snp() 	              // 
          ,upperTrk.tgl() 	              // 
          ,upperTrk.signed1Pt() 	        // (sign of charge)/Pt in c/GeV. Use pt() and sign() instead
          ,upperTrk.pt() 	              // Transverse momentum of the track in GeV/c
          ,upperTrk.p() 	                // Momentum in Gev/c
          ,upperTrk.eta() 	              // Pseudorapidity
          ,upperTrk.phi() 	              // Phi of the track, in radians within [0, 2pi)

          ,upperTrk.isWithinBeamPipe() 	  // Is the track within the beam pipe (= successfully propagated to a collision vertex)
          ,upperTrk.px() 	                // Momentum in x-direction in GeV/c
          ,upperTrk.py() 	                // Momentum in y-direction in GeV/c
          ,upperTrk.pz() 	                // Momentum in z-direction in GeV/c

          // // TracksCov
          ,upperTrk.sigmaY() 	    //Covariance matrix
          ,upperTrk.sigmaZ() 	    //Covariance matrix
          ,upperTrk.sigmaSnp() 	  //Covariance matrix
          ,upperTrk.sigmaTgl() 	  //Covariance matrix
          ,upperTrk.sigma1Pt() 	  //Covariance matrix
          ,upperTrk.rhoZY() 	    //Covariance matrix in compressed form
          ,upperTrk.rhoSnpY() 	 	//Covariance matrix in compressed form
          ,upperTrk.rhoSnpZ() 	 	//Covariance matrix in compressed form
          ,upperTrk.rhoTglY() 	 	//Covariance matrix in compressed form
          ,upperTrk.rhoTglZ() 	 	//Covariance matrix in compressed form
          ,upperTrk.rhoTglSnp() 	//Covariance matrix in compressed form
          ,upperTrk.rho1PtY() 	 	//Covariance matrix in compressed form
          ,upperTrk.rho1PtZ() 	 	//Covariance matrix in compressed form
          ,upperTrk.rho1PtSnp() 	//Covariance matrix in compressed form
          ,upperTrk.rho1PtTgl() 	//Covariance matrix in compressed form
          ,upperTrk.cYY() 	  //
          ,upperTrk.cZY() 	  
          ,upperTrk.cZZ() 	  
          ,upperTrk.cSnpY() 	  
          ,upperTrk.cSnpZ() 	  
          ,upperTrk.cSnpSnp()
          ,upperTrk.cTglY() 	  
          ,upperTrk.cTglZ() 	  
          ,upperTrk.cTglSnp()
          ,upperTrk.cTglTgl()
          ,upperTrk.c1PtY() 	  
          ,upperTrk.c1PtZ() 	  
          ,upperTrk.c1PtSnp()
          ,upperTrk.c1PtTgl()
          ,upperTrk.c1Pt21Pt2()

          // // TracksExtra
          ,upperTrk.tpcInnerParam() 	                    // Momentum at inner wall of the TPC
          ,upperTrk.flags() 	                            // Track flags. Run 2: see TrackFlagsRun2Enum | Run 3: see TrackFlags
          ,upperTrk.itsClusterSizes() 	                  // Clusters sizes, four bits per a layer, starting from the innermost
          ,upperTrk.tpcNClsFindable() 	                  // Findable TPC clusters for this track geometry
          ,upperTrk.tpcNClsFindableMinusFound() 	        // TPC Clusters: Findable - Found
          ,upperTrk.tpcNClsFindableMinusCrossedRows() 	  // TPC Clusters: Findable - crossed rows
          ,upperTrk.tpcNClsShared() 	                    // Number of shared TPC clusters
          ,upperTrk.trdPattern() 	                      // Contributor to the track on TRD layer in bits 0-5, starting from the innermost, bit 6 indicates a potentially split tracklet, bit 7 if the track crossed a padrow
          ,upperTrk.itsChi2NCl() 	                      // Chi2 / cluster for the ITS track segment
          ,upperTrk.tpcChi2NCl() 	                      // Chi2 / cluster for the TPC track segment
          ,upperTrk.trdChi2() 	                          // Chi2 for the TRD track segment
          ,upperTrk.tofChi2() 	                          // Chi2 for the TOF track segment
          ,upperTrk.tpcSignal() 	                        // dE/dx signal in the TPC
          ,upperTrk.trdSignal() 	                        // PID signal in the TRD
          ,upperTrk.length() 	                          // Track length
          ,upperTrk.tofExpMom() 	                        // TOF expected momentum obtained in tracking, used to compute the expected times
          ,upperTrk.trackEtaEmcal() 	                    // 
          ,upperTrk.trackPhiEmcal() 	                    // 
          ,upperTrk.trackTime() 	                        // Estimated time of the track in ns wrt collision().bc() or ambiguousupperTrk.bcSlice()[0]
          ,upperTrk.trackTimeRes() 	                    // Resolution of the track time in ns (see TrackFlags::TrackTimeResIsRange)
          
          // ,upperTrk.tpcDeltaTFwd() 	                    // Delta Forward of track time in TPC time bis
          // ,upperTrk.tpcDeltaTBwd() 	                    // Delta Backward of track time in TPC time bis
          ,upperTrk.pidForTracking() 	                  // PID hypothesis used during tracking. See the constants in the class PID in PID.h
          ,upperTrk.isPVContributor() 	                  // Run 3: Has this track contributed to the collision vertex fit

          ,upperTrkDetectorMask
          // ,upperTrk.hasITS() 	                          // Flag to check if track has a ITS match
          // ,upperTrk.hasTPC() 	                          // Flag to check if track has a TPC match
          // ,upperTrk.hasTRD() 	                          // Flag to check if track has a TRD match
          // ,upperTrk.hasTOF() 	                          // Flag to check if track has a TOF measurement

          ,upperTrk.tpcNClsFound() 	                    // Number of found TPC clusters
          ,upperTrk.tpcNClsCrossedRows() 	              // Number of crossed TPC Rows
          ,upperTrk.itsClusterMap() 	                    // ITS cluster map, one bit per a layer, starting from the innermost
          ,upperTrk.itsNCls() 	                          // Number of ITS clusters
          ,upperTrk.itsNClsInnerBarrel() 	              // Number of ITS clusters in the Inner Barrel
          // ,upperTrk.itsClsSizeInLayer() 	                // Size of the ITS cluster in a given layer
          // ,upperTrk.isITSAfterburner() 	                // If the track used the afterburner in the ITS
          // ,upperTrk.tofExpTime() 	                      // Expected time for the track to reach the TOF
          // ,upperTrk.tofExpTimePi() 	                    // Expected time for the track to reach the TOF under the pion hypothesis
          // ,upperTrk.tofExpTimeKa() 	                    // Expected time for the track to reach the TOF under the kaon hypothesis
          // ,upperTrk.tofExpTimePr() 	                    // Expected time for the track to reach the TOF under the proton hypothesis
          ,upperTrk.tpcCrossedRowsOverFindableCls() 	    // Ratio crossed rows over findable clusters
          ,upperTrk.tpcFoundOverFindableCls()            // Ratio of found over findable clusters
          ,upperTrk.tpcFractionSharedCls() 	            // Fraction of shared TPC clusters
          ,upperTrk.detectorMap() 	                      // Detector map version 1, see enum DetectorMapEnum
        
          // TracksQA
          // // ,o2::soa::Index 	GI 	globalIndex 	    int64_t 	
          // // ,o2::aod::trackqa::TrackId 	            //I 	    trackId 	        //int32 	track to which this QA information belongs
          ,upperTrkTpcTime0             //tpc only time0 (mTime0 in TPC track)
          ,upperTrkTpcdcaR              //tpc only DCAr
          ,upperTrkTpcdcaZ              //tpc only DCAz
          ,upperTrkTpcClusterByteMask   //tracklet bitmask - track defining 8 tracklets (152=8*19 rows) bit set if nCluster>thr (default 5)
          ,upperTrkTpcdEdxMax0R         //TPC dEdxQMax -ROC0/dEdx
          ,upperTrkTpcdEdxMax1R         //TPC dEdxQMax -ROC1/dEdx
          ,upperTrkTpcdEdxMax2R         //TPC dEdxQMax -ROC2/dEdx
          ,upperTrkTpcdEdxMax3R         //TPC dEdxQMax -ROC3/dEdx
          ,upperTrkTpcdEdxTot0R         //TPC dEdxQtot -ROC0/dEdx
          ,upperTrkTpcdEdxTot1R         //TPC dEdxQtot -ROC1/dEdx
          ,upperTrkTpcdEdxTot2R         //TPC dEdxQtot -ROC2/dEdx
          ,upperTrkTpcdEdxTot3R         //TPC dEdxQtot -ROC3/dEdx
          //
          ,upperTrk.beta()
          ,upperTrk.betaerror()
          ,upperTrk.mass()

          ,codeSqrtScaling(upperTrk.tpcNSigmaPi(), cfgPIDSigma0, cfgPIDSigma1, cfgPIDClampMin, cfgPIDClampMax)
          ,codeSqrtScaling(upperTrk.tofNSigmaPi(), cfgPIDSigma0, cfgPIDSigma1, cfgPIDClampMin, cfgPIDClampMax)
          ,codeSqrtScaling(upperTrk.tpcNSigmaKa(), cfgPIDSigma0, cfgPIDSigma1, cfgPIDClampMin, cfgPIDClampMax)
          ,codeSqrtScaling(upperTrk.tofNSigmaKa(), cfgPIDSigma0, cfgPIDSigma1, cfgPIDClampMin, cfgPIDClampMax)
          ,codeSqrtScaling(upperTrk.tpcNSigmaPr(), cfgPIDSigma0, cfgPIDSigma1, cfgPIDClampMin, cfgPIDClampMax)
          ,codeSqrtScaling(upperTrk.tofNSigmaPr(), cfgPIDSigma0, cfgPIDSigma1, cfgPIDClampMin, cfgPIDClampMax)
          ,codeSqrtScaling(upperTrk.tpcNSigmaEl(), cfgPIDSigma0, cfgPIDSigma1, cfgPIDClampMin, cfgPIDClampMax)
          ,codeSqrtScaling(upperTrk.tofNSigmaEl(), cfgPIDSigma0, cfgPIDSigma1, cfgPIDClampMin, cfgPIDClampMax)
          ,codeSqrtScaling(upperTrk.tpcNSigmaDe(), cfgPIDSigma0, cfgPIDSigma1, cfgPIDClampMin, cfgPIDClampMax)
          ,codeSqrtScaling(upperTrk.tofNSigmaDe(), cfgPIDSigma0, cfgPIDSigma1, cfgPIDClampMin, cfgPIDClampMax)
          
          ,upperTrk.tofSignal()

          ,upperTrk.dcaXY()   //	dcaXY	float	Impact parameter in XY of the track to the primary vertex
          ,upperTrk.dcaZ()	   //	dcaZ	float	Impact parameter in Z of the track to the primary vertex
          ,upperTrk.sigmaDcaXY2()  //		sigmaDcaXY2	float	Impact parameter sigma^2 in XY of the track to the primary vertex
          ,upperTrk.sigmaDcaZ2()   //		sigmaDcaZ2	float	Impact parameter sigma^2 in Z of the track to the primary verte

          ,upperTrkDeltaRefContParamY   
          ,upperTrkDeltaRefContParamZ   
          ,upperTrkDeltaRefContParamSnp 
          ,upperTrkDeltaRefContParamTgl 
          ,upperTrkDeltaRefContParamQ2Pt
          ,upperTrkDeltaRefGloParamY    
          ,upperTrkDeltaRefGloParamZ    
          ,upperTrkDeltaRefGloParamSnp  
          ,upperTrkDeltaRefGloParamTgl  
          ,upperTrkDeltaRefGloParamQ2Pt 
          ,upperTrkDeltaTOFdX
          ,upperTrkDeltaTOFdZ           

        ////////////////Lower Track Variables ////////////////////////////////////////////////////////////////////////////////////////////////
          ,lowerTrk.x() 	                // 
          ,lowerTrk.alpha() 	            // 
          ,lowerTrk.y() 	                // 
          ,lowerTrk.z() 	                // 
          ,lowerTrk.snp() 	              // 
          ,lowerTrk.tgl() 	              // 
          ,lowerTrk.signed1Pt() 	        // (sign of charge)/Pt in c/GeV. Use pt() and sign() instead
          ,lowerTrk.pt() 	              // Transverse momentum of the track in GeV/c
          ,lowerTrk.p() 	                // Momentum in Gev/c
          ,lowerTrk.eta() 	              // Pseudorapidity
          ,lowerTrk.phi() 	              // Phi of the track, in radians within [0, 2pi)

          ,lowerTrk.isWithinBeamPipe() 	  // Is the track within the beam pipe (= successfully propagated to a collision vertex)
          ,lowerTrk.px() 	                // Momentum in x-direction in GeV/c
          ,lowerTrk.py() 	                // Momentum in y-direction in GeV/c
          ,lowerTrk.pz() 	                // Momentum in z-direction in GeV/c

          // // TracksCov
          ,lowerTrk.sigmaY() 	    //Covariance matrix
          ,lowerTrk.sigmaZ() 	    //Covariance matrix
          ,lowerTrk.sigmaSnp() 	  //Covariance matrix
          ,lowerTrk.sigmaTgl() 	  //Covariance matrix
          ,lowerTrk.sigma1Pt() 	  //Covariance matrix
          ,lowerTrk.rhoZY() 	    //Covariance matrix in compressed form
          ,lowerTrk.rhoSnpY() 	 	//Covariance matrix in compressed form
          ,lowerTrk.rhoSnpZ() 	 	//Covariance matrix in compressed form
          ,lowerTrk.rhoTglY() 	 	//Covariance matrix in compressed form
          ,lowerTrk.rhoTglZ() 	 	//Covariance matrix in compressed form
          ,lowerTrk.rhoTglSnp() 	//Covariance matrix in compressed form
          ,lowerTrk.rho1PtY() 	 	//Covariance matrix in compressed form
          ,lowerTrk.rho1PtZ() 	 	//Covariance matrix in compressed form
          ,lowerTrk.rho1PtSnp() 	//Covariance matrix in compressed form
          ,lowerTrk.rho1PtTgl() 	//Covariance matrix in compressed form
          ,lowerTrk.cYY() 	  //
          ,lowerTrk.cZY() 	  
          ,lowerTrk.cZZ() 	  
          ,lowerTrk.cSnpY() 	  
          ,lowerTrk.cSnpZ() 	  
          ,lowerTrk.cSnpSnp()
          ,lowerTrk.cTglY() 	  
          ,lowerTrk.cTglZ() 	  
          ,lowerTrk.cTglSnp()
          ,lowerTrk.cTglTgl()
          ,lowerTrk.c1PtY() 	  
          ,lowerTrk.c1PtZ() 	  
          ,lowerTrk.c1PtSnp()
          ,lowerTrk.c1PtTgl()
          ,lowerTrk.c1Pt21Pt2()

          // // TracksExtra
          ,lowerTrk.tpcInnerParam() 	                    // Momentum at inner wall of the TPC
          ,lowerTrk.flags() 	                            // Track flags. Run 2: see TrackFlagsRun2Enum | Run 3: see TrackFlags
          ,lowerTrk.itsClusterSizes() 	                  // Clusters sizes, four bits per a layer, starting from the innermost
          ,lowerTrk.tpcNClsFindable() 	                  // Findable TPC clusters for this track geometry
          ,lowerTrk.tpcNClsFindableMinusFound() 	        // TPC Clusters: Findable - Found
          ,lowerTrk.tpcNClsFindableMinusCrossedRows() 	  // TPC Clusters: Findable - crossed rows
          ,lowerTrk.tpcNClsShared() 	                    // Number of shared TPC clusters
          ,lowerTrk.trdPattern() 	                      // Contributor to the track on TRD layer in bits 0-5, starting from the innermost, bit 6 indicates a potentially split tracklet, bit 7 if the track crossed a padrow
          ,lowerTrk.itsChi2NCl() 	                      // Chi2 / cluster for the ITS track segment
          ,lowerTrk.tpcChi2NCl() 	                      // Chi2 / cluster for the TPC track segment
          ,lowerTrk.trdChi2() 	                          // Chi2 for the TRD track segment
          ,lowerTrk.tofChi2() 	                          // Chi2 for the TOF track segment
          ,lowerTrk.tpcSignal() 	                        // dE/dx signal in the TPC
          ,lowerTrk.trdSignal() 	                        // PID signal in the TRD
          ,lowerTrk.length() 	                          // Track length
          ,lowerTrk.tofExpMom() 	                        // TOF expected momentum obtained in tracking, used to compute the expected times
          ,lowerTrk.trackEtaEmcal() 	                    // 
          ,lowerTrk.trackPhiEmcal() 	                    // 
          ,lowerTrk.trackTime() 	                        // Estimated time of the track in ns wrt collision().bc() or ambiguouslowerTrk.bcSlice()[0]
          ,lowerTrk.trackTimeRes() 	                    // Resolution of the track time in ns (see TrackFlags::TrackTimeResIsRange)
          
          // ,lowerTrk.tpcDeltaTFwd() 	                    // Delta Forward of track time in TPC time bis
          // ,lowerTrk.tpcDeltaTBwd() 	                    // Delta Backward of track time in TPC time bis
          ,lowerTrk.pidForTracking() 	                  // PID hypothesis used during tracking. See the constants in the class PID in PID.h
          ,lowerTrk.isPVContributor() 	                  // Run 3: Has this track contributed to the collision vertex fit

          ,lowerTrkDetectorMask
          // ,lowerTrk.hasITS() 	                          // Flag to check if track has a ITS match
          // ,lowerTrk.hasTPC() 	                          // Flag to check if track has a TPC match
          // ,lowerTrk.hasTRD() 	                          // Flag to check if track has a TRD match
          // ,lowerTrk.hasTOF() 	                          // Flag to check if track has a TOF measurement

          ,lowerTrk.tpcNClsFound() 	                    // Number of found TPC clusters
          ,lowerTrk.tpcNClsCrossedRows() 	              // Number of crossed TPC Rows
          ,lowerTrk.itsClusterMap() 	                    // ITS cluster map, one bit per a layer, starting from the innermost
          ,lowerTrk.itsNCls() 	                          // Number of ITS clusters
          ,lowerTrk.itsNClsInnerBarrel() 	              // Number of ITS clusters in the Inner Barrel
          // ,lowerTrk.itsClsSizeInLayer() 	                // Size of the ITS cluster in a given layer
          // ,lowerTrk.isITSAfterburner() 	                // If the track used the afterburner in the ITS
          // ,lowerTrk.tofExpTime() 	                      // Expected time for the track to reach the TOF
          // ,lowerTrk.tofExpTimePi() 	                    // Expected time for the track to reach the TOF under the pion hypothesis
          // ,lowerTrk.tofExpTimeKa() 	                    // Expected time for the track to reach the TOF under the kaon hypothesis
          // ,lowerTrk.tofExpTimePr() 	                    // Expected time for the track to reach the TOF under the proton hypothesis
          ,lowerTrk.tpcCrossedRowsOverFindableCls() 	    // Ratio crossed rows over findable clusters
          ,lowerTrk.tpcFoundOverFindableCls()            // Ratio of found over findable clusters
          ,lowerTrk.tpcFractionSharedCls() 	            // Fraction of shared TPC clusters
          ,lowerTrk.detectorMap() 	                      // Detector map version 1, see enum DetectorMapEnum
        
          // TracksQA
          // // ,o2::soa::Index 	GI 	globalIndex 	    int64_t 	
          // // ,o2::aod::trackqa::TrackId 	            //I 	    trackId 	        //int32 	track to which this QA information belongs
          ,lowerTrkTpcTime0             //tpc only time0 (mTime0 in TPC track)
          ,lowerTrkTpcdcaR              //tpc only DCAr
          ,lowerTrkTpcdcaZ              //tpc only DCAz
          ,lowerTrkTpcClusterByteMask   //tracklet bitmask - track defining 8 tracklets (152=8*19 rows) bit set if nCluster>thr (default 5)
          ,lowerTrkTpcdEdxMax0R         //TPC dEdxQMax -ROC0/dEdx
          ,lowerTrkTpcdEdxMax1R         //TPC dEdxQMax -ROC1/dEdx
          ,lowerTrkTpcdEdxMax2R         //TPC dEdxQMax -ROC2/dEdx
          ,lowerTrkTpcdEdxMax3R         //TPC dEdxQMax -ROC3/dEdx
          ,lowerTrkTpcdEdxTot0R         //TPC dEdxQtot -ROC0/dEdx
          ,lowerTrkTpcdEdxTot1R         //TPC dEdxQtot -ROC1/dEdx
          ,lowerTrkTpcdEdxTot2R         //TPC dEdxQtot -ROC2/dEdx
          ,lowerTrkTpcdEdxTot3R         //TPC dEdxQtot -ROC3/dEdx
          //
          ,lowerTrk.beta()
          ,lowerTrk.betaerror()
          ,lowerTrk.mass()

          ,codeSqrtScaling(lowerTrk.tpcNSigmaPi(), cfgPIDSigma0, cfgPIDSigma1, cfgPIDClampMin, cfgPIDClampMax)
          ,codeSqrtScaling(lowerTrk.tofNSigmaPi(), cfgPIDSigma0, cfgPIDSigma1, cfgPIDClampMin, cfgPIDClampMax)
          ,codeSqrtScaling(lowerTrk.tpcNSigmaKa(), cfgPIDSigma0, cfgPIDSigma1, cfgPIDClampMin, cfgPIDClampMax)
          ,codeSqrtScaling(lowerTrk.tofNSigmaKa(), cfgPIDSigma0, cfgPIDSigma1, cfgPIDClampMin, cfgPIDClampMax)
          ,codeSqrtScaling(lowerTrk.tpcNSigmaPr(), cfgPIDSigma0, cfgPIDSigma1, cfgPIDClampMin, cfgPIDClampMax)
          ,codeSqrtScaling(lowerTrk.tofNSigmaPr(), cfgPIDSigma0, cfgPIDSigma1, cfgPIDClampMin, cfgPIDClampMax)
          ,codeSqrtScaling(lowerTrk.tpcNSigmaEl(), cfgPIDSigma0, cfgPIDSigma1, cfgPIDClampMin, cfgPIDClampMax)
          ,codeSqrtScaling(lowerTrk.tofNSigmaEl(), cfgPIDSigma0, cfgPIDSigma1, cfgPIDClampMin, cfgPIDClampMax)
          ,codeSqrtScaling(lowerTrk.tpcNSigmaDe(), cfgPIDSigma0, cfgPIDSigma1, cfgPIDClampMin, cfgPIDClampMax)
          ,codeSqrtScaling(lowerTrk.tofNSigmaDe(), cfgPIDSigma0, cfgPIDSigma1, cfgPIDClampMin, cfgPIDClampMax)
          
          ,lowerTrk.tofSignal()

          ,lowerTrk.dcaXY()   //	dcaXY	float	Impact parameter in XY of the track to the primary vertex
          ,lowerTrk.dcaZ()	   //	dcaZ	float	Impact parameter in Z of the track to the primary vertex
          ,lowerTrk.sigmaDcaXY2()  //		sigmaDcaXY2	float	Impact parameter sigma^2 in XY of the track to the primary vertex
          ,lowerTrk.sigmaDcaZ2()   //		sigmaDcaZ2	float	Impact parameter sigma^2 in Z of the track to the primary verte

          ,lowerTrkDeltaRefContParamY   
          ,lowerTrkDeltaRefContParamZ   
          ,lowerTrkDeltaRefContParamSnp 
          ,lowerTrkDeltaRefContParamTgl 
          ,lowerTrkDeltaRefContParamQ2Pt
          ,lowerTrkDeltaRefGloParamY    
          ,lowerTrkDeltaRefGloParamZ    
          ,lowerTrkDeltaRefGloParamSnp  
          ,lowerTrkDeltaRefGloParamTgl  
          ,lowerTrkDeltaRefGloParamQ2Pt 
          ,lowerTrkDeltaTOFdX
          ,lowerTrkDeltaTOFdZ                     

          //Sum and difference variables
          ,fSumPtPair     
          ,fSumQPtPair    
          ,fSumTglPair    
          ,fSumDcaXYPair
          ,fSumAlphaPair
          
          ,fDiffPtPair    
          ,fDiffQPtPair   
          ,fDiffTglPair   
          ,fDiffDcaXYPair 
          ,fDiffAlphaPair

          ,time1
          ,TFidThis1
          ,bcInTF1

          ,time2
          ,TFidThis2
          ,bcInTF2
        );
        cosmicPairCounter++;
      }//lower track loop ends
    }//upper track loop ends 

    LOG(info)<<"DEBUG ::  trackCounterDebugger = "<<trackCounterDebugger<<" :: ambgTracks = "<<ambgTracks.size();
    
    LOG(info)<<"DEBUG :: cosmicPairsCounted = "<<cosmicPairCounter<<" :: gCosmicPairChecker = "<<gCosmicPairChecker;
    PrintTime(StartCosmicMuon, Form("DEBUG :: df_%ld :: Cosmic Muon Time :: ",dfCount));
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
    gZ0Checker        += (z0Counter);
    gCosmicPairChecker+= (cosmicPairCounter);
/**/
    PrintTime(Start1, Form("DEBUG :: df_%ld :: DF End    :: DF Read Time :: ",dfCount));
    PrintTime(Start0, Form("DEBUG :: df_%ld :: DF End    :: Elapsed Time :: ",dfCount));
    LOG(info)<<"DEBUG ::";

  }//Process function ends

  void processNothing(aod::Collisions const& ){
    return;
  }
  PROCESS_SWITCH(trackqapidderivedata, processNothing, "processNothing", true);

  void processWithoutCentrality(myBCTable const& BCs
               ,myCollisionsWithoutCent const& collisions
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
              ){
    executeProcess<false>(BCs, collisions, tracks, tracksQA, Origins, ambgTracks, FV0As, FT0s, V0s, D0s, occIdxTable, occTables, trackMeanOccs);
  }

  PROCESS_SWITCH(trackqapidderivedata, processWithoutCentrality, "processWithoutCentrality", false);

  void processWithCentrality( myBCTable const& BCs
               ,myCollisionsWithCent const& collisions
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
              ){
    executeProcess<true>(BCs, collisions, tracks, tracksQA, Origins, ambgTracks, FV0As, FT0s, V0s, D0s, occIdxTable, occTables, trackMeanOccs);
  }
  PROCESS_SWITCH(trackqapidderivedata, processWithCentrality, "processWithCentrality", false);
};


WorkflowSpec defineDataProcessing(ConfigContext const& cfgc){
  return WorkflowSpec{ 
                      adaptAnalysisTask<trackqapidderivedata>(cfgc)
  };
}
