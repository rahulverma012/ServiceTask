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
/// \file   DerivedPIDTables.h
/// \brief  HeaderFiles for derived data creation
/// \author Rahul Verma (rahul.verma@iitb.ac.in) :: Marian I Ivanov (marian.ivanov@cern.ch) :: Yash Patley (yash.patley@cern.ch)

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/FT0Corrected.h"
#include "Common/DataModel/OccupancyTables.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/TableProducer/PID/pidTOFBase.h"

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"

#include <vector>

namespace o2::aod
{
namespace drcolmn
{
DECLARE_SOA_COLUMN(GlobalIndex, gIDX, uint64_t); //o2-linter: disable=name/o2-column ()
DECLARE_SOA_COLUMN(OriginalIndex, oIDX, uint64_t);
DECLARE_SOA_COLUMN(GBCid, gbcId, uint64_t);
DECLARE_SOA_COLUMN(GCollId, gCollId, uint64_t);
DECLARE_SOA_COLUMN(GDFId, gDfId, int64_t);

} // namespace drcolmn
} // namespace o2::aod

// Defining Collision Tables
namespace o2::aod
{
namespace drcolmn
{
DECLARE_SOA_COLUMN(TfId, tfId, int64_t);
DECLARE_SOA_COLUMN(BCinTF, bcInTF, int);
DECLARE_SOA_COLUMN(BCinTFFromEvSel, bcInTFFromEvSel, int);
DECLARE_SOA_COLUMN(PosZ, posZ, float);                // D 	posZ     	float
DECLARE_SOA_COLUMN(CollTime, collTime, float);        // D 	collTime 	float
DECLARE_SOA_COLUMN(IsValidTimeA, isValidTimeA, bool); // D 	isValidTimeA 	bool
DECLARE_SOA_COLUMN(IsValidTimeC, isValidTimeC, bool); // D 	isValidTimeC 	bool
DECLARE_SOA_COLUMN(IsValidTime, isValidTime, bool);   // D 	isValidTime 	bool
DECLARE_SOA_COLUMN(SumAmpA, sumAmpA, float);          // D 	sumAmpA 	float
DECLARE_SOA_COLUMN(SumAmpC, sumAmpC, float);          // D 	sumAmpC 	float

DECLARE_SOA_COLUMN(NTrack_PVC, nTrack_PVC, int);
DECLARE_SOA_COLUMN(NTrack_ITS, nTrack_ITS, int);
DECLARE_SOA_COLUMN(NTrack_TPC, nTrack_TPC, int);
DECLARE_SOA_COLUMN(NTrack_TRD, nTrack_TRD, int);
DECLARE_SOA_COLUMN(NTrack_TOF, nTrack_TOF, int);
DECLARE_SOA_COLUMN(NTrackTPC_A, nTrackTPC_A, int);
DECLARE_SOA_COLUMN(NTrackTPC_C, nTrackTPC_C, int);
DECLARE_SOA_COLUMN(NTrackITS_TPC, nTrackITS_TPC, int);
DECLARE_SOA_COLUMN(NTrackITS_TPC_A, nTrackITS_TPC_A, int);
DECLARE_SOA_COLUMN(NTrackITS_TPC_C, nTrackITS_TPC_C, int);
DECLARE_SOA_COLUMN(NTrackSize, nTrackSize, int);

DECLARE_SOA_COLUMN(V0Radius, v0radius, float);                 // D 	v0radius 	float 	V0 decay radius (2D, centered at zero)
DECLARE_SOA_COLUMN(MLambda, mLambda, float);                   // D 	mLambda 	float 	mass under lambda hypothesis
DECLARE_SOA_COLUMN(MAntiLambda, mAntiLambda, float);           // D 	mAntiLambda 	float 	mass under antilambda hypothesis
DECLARE_SOA_COLUMN(MK0Short, mK0Short, float);                 // D 	mK0Short 	float 	mass under K0short hypothesis
DECLARE_SOA_COLUMN(MGamma, mGamma, float);                     // D 	mGamma 	float 	mass under gamma hypothesis
DECLARE_SOA_COLUMN(MHypertriton, mHypertriton, float);         // D 	mHypertriton 	float 	mass under hypertriton hypothesis
DECLARE_SOA_COLUMN(MAntiHypertriton, mAntiHypertriton, float); // D 	mAntiHypertriton 	float 	mass under antihypertriton hypothesis
DECLARE_SOA_COLUMN(NegativePt, negativept, float);             // D 	negativept 	float 	negative daughter pT
DECLARE_SOA_COLUMN(PositivePt, positivept, float);             // D 	positivept 	float 	positive daughter pT
DECLARE_SOA_COLUMN(NegativeEta, negativeeta, float);           // D 	negativeeta 	float 	negative daughter eta
DECLARE_SOA_COLUMN(NegativePhi, negativephi, float);           // D 	negativephi 	float 	negative daughter phi
DECLARE_SOA_COLUMN(PositiveEta, positiveeta, float);           // D 	positiveeta 	float 	positive daughter eta
DECLARE_SOA_COLUMN(PositivePhi, positivephi, float);           // D 	positivephi 	float 	positive daughter phi
DECLARE_SOA_COLUMN(IsStandardV0, isStandardV0, bool);          // D 	isStandardV0 	bool 	is standard V0
DECLARE_SOA_COLUMN(IsPhotonTPConly, isPhotonTPConly, bool);    // D 	isPhotonTPConly 	bool 	is tpc-only photon V0

DECLARE_SOA_COLUMN(Alpha, alpha, float);
DECLARE_SOA_COLUMN(QtArm, qtarm, float);
DECLARE_SOA_COLUMN(DistOverTotMom, distovertotmom, float);
DECLARE_SOA_COLUMN(PsiPair, psipair, float);
DECLARE_SOA_COLUMN(PFracPos, pfracpos, float);
DECLARE_SOA_COLUMN(PFracNeg, pfracneg, float);

DECLARE_SOA_COLUMN(T0AC, t0AC, float);                          //	Collision time (A+C)/2
DECLARE_SOA_COLUMN(T0ACorrectedValid, t0ACorrectedValid, bool); //	Was T0ACorrected computable?
DECLARE_SOA_COLUMN(T0CCorrectedValid, t0CCorrectedValid, bool); //	Was T0CCorrected computable?
DECLARE_SOA_COLUMN(T0ACValid, t0ACValid, bool);                 //	Was T0AC computable?
DECLARE_SOA_COLUMN(T0resolution, t0resolution, float);          //	Was T0CCorrected computable?

DECLARE_SOA_COLUMN(WeightDSPtV0, weightDSPtV0, float);
DECLARE_SOA_COLUMN(WeightDSPtLambda, weightDSPtLambda, float);
DECLARE_SOA_COLUMN(WeightDSPtAntiLambda, weightDSPtAntiLambda, float);
DECLARE_SOA_COLUMN(WeightDSPtK0Short, weightDSPtK0Short, float);
DECLARE_SOA_COLUMN(WeightDSPtGamma, weightDSPtGamma, float);
DECLARE_SOA_COLUMN(WeightDSPtHypertriton, weightDSPtHypertriton, float);
DECLARE_SOA_COLUMN(WeightDSPtAntiHypertriton, weightDSPtAntiHypertriton, float);

DECLARE_SOA_COLUMN(TriggerMaskV0, triggerMaskV0, int);
DECLARE_SOA_COLUMN(TriggerMaskV0PT, triggerMaskV0PT, int);

DECLARE_SOA_COLUMN(RSecondaryVertex, rSecondaryVertex, float);
DECLARE_SOA_COLUMN(DecayLength, decayLength, float);
DECLARE_SOA_COLUMN(DecayLengthXY, decayLengthXY, float);
DECLARE_SOA_COLUMN(DecayLengthNormalised, decayLengthNormalised, float);
DECLARE_SOA_COLUMN(DecayLengthXYNormalised, decayLengthXYNormalised, float);
DECLARE_SOA_COLUMN(ImpactParameterNormalised0, impactParameterNormalised0, float);
DECLARE_SOA_COLUMN(PtProng0, ptProng0, float);
DECLARE_SOA_COLUMN(Pt2Prong0, pt2Prong0, float);
// DECLARE_SOA_COLUMN(PVectorProng0, pVectorProng0, 	std::array<float,3>);
DECLARE_SOA_COLUMN(ImpactParameterNormalised1, impactParameterNormalised1, float);
DECLARE_SOA_COLUMN(PtProng1, ptProng1, float);
DECLARE_SOA_COLUMN(Pt2Prong1, pt2Prong1, float);

DECLARE_SOA_COLUMN(ImpactParameterProduct, impactParameterProduct, float);
// // // ,o2::aod::hf_cand_2prong::CosThetaStar 	                //D 	cosThetaStar 	float
DECLARE_SOA_COLUMN(ImpactParameterProngSqSum, impactParameterProngSqSum, float);
// // ,o2::aod::hf_cand_2prong::ImpactParameterProngSqSum 	  //D 	impactParameterProngSqSum 	float
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Pt2, pt2, float);
DECLARE_SOA_COLUMN(P, p, float);
DECLARE_SOA_COLUMN(P2, p2, float);
// // // // ,o2::aod::hf_cand::PVector 	                            //D 	pVector 	std::array<float,3>
DECLARE_SOA_COLUMN(CPA, cpa, float);
DECLARE_SOA_COLUMN(CPAXY, cpaXY, float);
// // // ,o2::aod::hf_cand::Ct 	                                //D 	ct 	float
DECLARE_SOA_COLUMN(ImpactParameterXY, impactParameterXY, float);
DECLARE_SOA_COLUMN(MaxNormalisedDeltaIP, maxNormalisedDeltaIP, float);
DECLARE_SOA_COLUMN(Px, px, float);
DECLARE_SOA_COLUMN(Py, py, float);
DECLARE_SOA_COLUMN(Pz, pz, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(Y, y, float);
DECLARE_SOA_COLUMN(E, e, float);

DECLARE_SOA_COLUMN(HEntriesTrkPt, hEntriesTrkPt, uint64_t);
DECLARE_SOA_COLUMN(HEntriesTrkQPt, hEntriesTrkQPt, uint64_t);

DECLARE_SOA_COLUMN(HEntriesV0Pt, hEntriesV0Pt, uint64_t);
DECLARE_SOA_COLUMN(HEntriesLambdaPt, hEntriesLambdaPt, uint64_t);
DECLARE_SOA_COLUMN(HEntriesAntiLambdaPt, hEntriesAntiLambdaPt, uint64_t);
DECLARE_SOA_COLUMN(HEntriesK0ShortPt, hEntriesK0ShortPt, uint64_t);

DECLARE_SOA_COLUMN(HEntriesGammaPt, hEntriesGammaPt, uint64_t);
DECLARE_SOA_COLUMN(HEntriesHypertritonPt, hEntriesHypertritonPt, uint64_t);
DECLARE_SOA_COLUMN(HEntriesAntiHypertritonPt, hEntriesAntiHypertritonPt, uint64_t);

DECLARE_SOA_COLUMN(MPhi, mPhi, float);

DECLARE_SOA_COLUMN(PosKaonPx, posKaonPx, float);
DECLARE_SOA_COLUMN(PosKaonPy, posKaonPy, float);
DECLARE_SOA_COLUMN(PosKaonPz, posKaonPz, float);
DECLARE_SOA_COLUMN(PosKaonPt, posKaonPt, float);
DECLARE_SOA_COLUMN(PosKaonEta, posKaonEta, float);
DECLARE_SOA_COLUMN(PosKaonY, posKaonY, float);
DECLARE_SOA_COLUMN(PosKaonIsPVContributor, posKaonIsPVContributor, bool);
DECLARE_SOA_COLUMN(PosKaonIsGlobalTrack, posKaonIsGlobalTrack, bool);

DECLARE_SOA_COLUMN(NegKaonPx, negKaonPx, float);
DECLARE_SOA_COLUMN(NegKaonPy, negKaonPy, float);
DECLARE_SOA_COLUMN(NegKaonPz, negKaonPz, float);
DECLARE_SOA_COLUMN(NegKaonPt, negKaonPt, float);
DECLARE_SOA_COLUMN(NegKaonEta, negKaonEta, float);
DECLARE_SOA_COLUMN(NegKaonY, negKaonY, float);
DECLARE_SOA_COLUMN(NegKaonIsPVContributor, negKaonIsPVContributor, bool);
DECLARE_SOA_COLUMN(NegKaonIsGlobalTrack, negKaonIsGlobalTrack, bool);

DECLARE_SOA_COLUMN(WeightDSPtPhi, weightDSPtPhi, float);
DECLARE_SOA_COLUMN(TriggerMaskPhi, triggerMaskPhi, int);
DECLARE_SOA_COLUMN(HEntriesPhiPt, hEntriesPhiPt, uint64_t);

DECLARE_SOA_COLUMN(WeightDSPtD0, weightDSPtD0, float);
DECLARE_SOA_COLUMN(WeightDSPtD0Bar, weightDSPtD0Bar, float);
DECLARE_SOA_COLUMN(TriggerMaskD0, triggerMaskD0, int);

DECLARE_SOA_COLUMN(HEntriesD0Pt, hEntriesD0Pt, uint64_t);
DECLARE_SOA_COLUMN(HEntriesD0BarPt, hEntriesD0BarPt, uint64_t);

DECLARE_SOA_COLUMN(Pol1FitTrackPt, pol1FitTrackPt, float);
DECLARE_SOA_COLUMN(Pol1FitTrackQPt, pol1FitTrackQPt, float);

DECLARE_SOA_COLUMN(Pol1FitV0Pt, pol1FitV0Pt, float);
DECLARE_SOA_COLUMN(Pol1FitLambdaPt, pol1FitLambdaPt, float);
DECLARE_SOA_COLUMN(Pol1FitAntiLambdaPt, pol1FitAntiLambdaPt, float);
DECLARE_SOA_COLUMN(Pol1FitK0ShortPt, pol1FitK0ShortPt, float);

DECLARE_SOA_COLUMN(Pol1FitPhiPt, pol1FitPhiPt, float);

DECLARE_SOA_COLUMN(Pol1FitD0Pt, pol1FitD0Pt, float);
DECLARE_SOA_COLUMN(Pol1FitD0BarPt, pol1FitD0BarPt, float);

DECLARE_SOA_COLUMN(GZ0Id, gZ0Id, int64_t);

DECLARE_SOA_COLUMN(GPosDauId, gPosDauId, int64_t);
DECLARE_SOA_COLUMN(GNegDauId, gNegDauId, int64_t);
DECLARE_SOA_COLUMN(DrPosDauId, drPosDauId, int64_t);
DECLARE_SOA_COLUMN(DrNegDauId, drNegDauId, int64_t);

DECLARE_SOA_COLUMN(PosTrackId, posTrackId, int64_t); // original Id in single DF in original Tracks
DECLARE_SOA_COLUMN(NegTrackId, negTrackId, int64_t); // original Id in single DF in original Tracks

DECLARE_SOA_COLUMN(DecayFlag, decayFlag, uint8_t); // muon decay or electron decay
DECLARE_SOA_COLUMN(MZ0, mZ0, float);
DECLARE_SOA_COLUMN(PosDauPx, posDauPx, float);
DECLARE_SOA_COLUMN(PosDauPy, posDauPy, float);
DECLARE_SOA_COLUMN(PosDauPz, posDauPz, float);
DECLARE_SOA_COLUMN(PosDauPt, posDauPt, float);
DECLARE_SOA_COLUMN(PosDauEta, posDauEta, float);
DECLARE_SOA_COLUMN(NegDauPx, negDauPx, float);
DECLARE_SOA_COLUMN(NegDauPy, negDauPy, float);
DECLARE_SOA_COLUMN(NegDauPz, negDauPz, float);
DECLARE_SOA_COLUMN(NegDauPt, negDauPt, float);
DECLARE_SOA_COLUMN(NegDauEta, negDauEta, float);

DECLARE_SOA_COLUMN(GUBCid, gubcId, uint64_t);
DECLARE_SOA_COLUMN(GUCollId, guCollId, uint64_t);

DECLARE_SOA_COLUMN(GLBCid, glbcId, uint64_t);
DECLARE_SOA_COLUMN(GLCollId, glCollId, uint64_t);

DECLARE_SOA_COLUMN(UTimestamp, uTimestamp, int64_t);//std::vector<int64_t>);
DECLARE_SOA_COLUMN(LTimestamp, lTimestamp, int64_t);//std::vector<int64_t>);

DECLARE_SOA_COLUMN(UTfId, uTfId, int64_t);// std::vector<int64_t>);
DECLARE_SOA_COLUMN(LTfId, lTfId, int64_t);// std::vector<int64_t>);

DECLARE_SOA_COLUMN(UBCinTF, uBCinTF, int64_t);// std::vector<int>);
DECLARE_SOA_COLUMN(LBCinTF, lBCinTF, int64_t);// std::vector<int>);

DECLARE_SOA_COLUMN(GCosmicPairId, gCosmicPairId, int64_t);

DECLARE_SOA_COLUMN(GUpperTrackId, gUpperTrackId, int64_t);
DECLARE_SOA_COLUMN(GLowerTrackId, gLowerTrackId, int64_t);

DECLARE_SOA_COLUMN(DrUpperTrackId, drUpperTrackId, int64_t);
DECLARE_SOA_COLUMN(DrLowerTrackId, drLowerTrackId, int64_t);

DECLARE_SOA_COLUMN(UpperTrackId, upperTrackId, int64_t);
DECLARE_SOA_COLUMN(LowerTrackId, lowerTrackId, int64_t);

// DECLARE_SOA_COLUMN(Amplitude   , amplitude 	  ,std::vector<float>); 	Amplitudes of non-zero channels. The channel IDs are given in Channel (at the same index)
} // namespace drcolmn

DECLARE_SOA_TABLE(DrCollisions, "AOD", "DRCOLLISION", o2::soa::Index<>,
                  o2::aod::drcolmn::GlobalIndex,
                  o2::aod::drcolmn::OriginalIndex,
                  o2::aod::drcolmn::GBCid,
                  o2::aod::drcolmn::GDFId,
                  o2::aod::collision::BCId,
                  o2::aod::collision::PosX,
                  o2::aod::collision::PosY,
                  o2::aod::collision::PosZ,
                  o2::aod::collision::CovXX,
                  o2::aod::collision::CovXY,
                  o2::aod::collision::CovXZ,
                  o2::aod::collision::CovYY,
                  o2::aod::collision::CovYZ,
                  o2::aod::collision::CovZZ,
                  o2::aod::collision::Flags,
                  o2::aod::collision::Chi2,
                  o2::aod::collision::NumContrib,
                  o2::aod::collision::CollisionTime,
                  o2::aod::collision::CollisionTimeRes,
                  o2::aod::timestamp::Timestamp,
                  o2::aod::drcolmn::TfId,
                  o2::aod::drcolmn::BCinTF,
                  o2::aod::mult::MultFV0A, //		multFV0A 	float
                  o2::aod::mult::MultFV0C, //		multFV0C 	float
                  // o2::aod::mult::MultFV0M, //	D 	multFV0M 	float
                  o2::aod::mult::MultFT0A, //		multFT0A 	float
                  o2::aod::mult::MultFT0C, //		multFT0C 	float
                  // o2::aod::mult::MultFT0M, //	D 	multFT0M 	float
                  o2::aod::mult::MultFDDA, //		multFDDA 	float
                  o2::aod::mult::MultFDDC, //		multFDDC 	float
                  // o2::aod::mult::MultFDDM, 	D 	multFDDM 	float
                  o2::aod::mult::MultZNA,  //		multZNA 	float
                  o2::aod::mult::MultZNC,  //		multZNC 	float
                  o2::aod::mult::MultZEM1, // 		multZEM1 	float
                  o2::aod::mult::MultZEM2, // 		multZEM2 	float
                  o2::aod::mult::MultZPA,  //		multZPA 	float
                  o2::aod::mult::MultZPC,  //		multZPC 	float
                  drcolmn::NTrack_PVC,
                  drcolmn::NTrack_ITS,
                  drcolmn::NTrack_TPC,
                  drcolmn::NTrack_TRD,
                  drcolmn::NTrack_TOF,
                  drcolmn::NTrackTPC_A,
                  drcolmn::NTrackTPC_C,
                  drcolmn::NTrackITS_TPC,
                  drcolmn::NTrackITS_TPC_A,
                  drcolmn::NTrackITS_TPC_C,
                  drcolmn::NTrackSize,

                  aod::mult::MultNTracksHasITS,
                  aod::mult::MultNTracksHasTPC,
                  aod::mult::MultNTracksHasTOF,
                  aod::mult::MultNTracksHasTRD,
                  aod::mult::MultNTracksITSOnly,
                  aod::mult::MultNTracksTPCOnly,
                  aod::mult::MultNTracksITSTPC,
                  aod::mult::MultAllTracksTPCOnly,

                  aod::evsel::Sel7,              //  sel7	bool	Event selection decision based on V0A & V0C
                  aod::evsel::Sel8,              //  sel8	bool	Event selection decision based on TVX
                  aod::evsel::FoundBCId,         // I	foundBCId	int	BC entry index in BCs table (-1 if doesn't exist)
                  aod::evsel::FoundFT0Id,        // I	foundFT0Id	int	FT0 entry index in FT0s table (-1 if doesn't exist)
                  aod::evsel::FoundFV0Id,        // I	foundFV0Id	int	FV0 entry index in FV0As table (-1 if doesn't exist)
                  aod::evsel::FoundFDDId,        // I	foundFDDId	int	FDD entry index in FDDs table (-1 if doesn't exist)
                  aod::evsel::FoundZDCId,        // I	foundZDCId	int	ZDC entry index in ZDCs table (-1 if doesn't exist)
                  aod::drcolmn::BCinTFFromEvSel, // the next line  aod::evsel::BcInTF was in conflict with the aod::drcolmn::BCinTF hence using this column
                  // aod::evsel::BcInTF,		                //  bcInTF	int	Position of a (found) bunch crossing inside a given timeframe
                  aod::evsel::NumTracksInTimeRange,  // trackOccupancyInTimeRange	int	Occupancy in specified time interval by a number of tracks from nearby collisions
                  aod::evsel::SumAmpFT0CInTimeRange, // ft0cOccupancyInTimeRange	float	Occupancy in specified time interval by a sum of FT0C amplitudes from nearby collisions

                  aod::ft0::T0ACorrected,     // t0ACorrected	    //float	Collision time A-side, corrected with primary vertex
                  aod::ft0::T0CCorrected,     // t0CCorrected	    //float	Collision time C-side, corrected with primary vertex
                  drcolmn::T0AC,              // t0AC	            //float	Collision time (A+C)/2
                  drcolmn::T0ACorrectedValid, // t0ACorrectedValid	//bool	Was T0ACorrected computable?
                  drcolmn::T0CCorrectedValid, // t0CCorrectedValid	//bool	Was T0CCorrected computable?
                  drcolmn::T0ACValid,         // t0ACValid	        //bool	Was T0AC computable?
                  drcolmn::T0resolution,      // t0resolution	    //float	Was T0CCorrected computable?

                  o2::aod::cent::CentFV0A,
                  o2::aod::cent::CentFT0M,
                  o2::aod::cent::CentFT0A,
                  o2::aod::cent::CentFT0C,
                  o2::aod::cent::CentFDDM,
                  o2::aod::cent::CentNTPV

);

using DrCollision = DrCollisions::iterator;
// //Defining Collision Tables
DECLARE_SOA_TABLE(DrFV0As, "AOD", "DRFV0As", o2::soa::Index<>,
                  o2::aod::drcolmn::GlobalIndex,
                  o2::aod::drcolmn::OriginalIndex,
                  o2::aod::drcolmn::GBCid,

                  o2::aod::fv0a::BCId,        // bcId 	int32 	BC index
                  o2::aod::fv0a::Amplitude,   // amplitude 	std::vector<float> 	Amplitudes of non-zero channels. The channel IDs are given in Channel (at the same index)
                  o2::aod::fv0a::Channel,     // channel 	std::vector<uint8_t> 	Channel IDs which had non-zero amplitudes. There are at maximum 48 channels.
                  o2::aod::fv0a::Time,        // time 	float 	Time in ns
                  o2::aod::fv0a::TriggerMask, // triggerMask 	uint8_t
                  // o2::aod::timestamp::Timestamp,
                  o2::aod::drcolmn::TfId,
                  o2::aod::drcolmn::BCinTF);

DECLARE_SOA_TABLE(DrFT0s, "AOD", "DRFT0s", o2::soa::Index<>,
                  o2::aod::drcolmn::GlobalIndex,
                  o2::aod::drcolmn::OriginalIndex,
                  o2::aod::drcolmn::GBCid,

                  o2::aod::ft0::BCId,             // I 	bcId 	int32 	BC index
                  o2::aod::ft0::AmplitudeA,       // amplitudeA 	std::vector<float> 	Amplitudes of non-zero channels on the A-side. The channel IDs are given in ChannelA (at the same index)
                  o2::aod::ft0::ChannelA,         // channelA 	std::vector<uint8_t> 	Channel IDs on the A side which had non-zero amplitudes. There are at maximum 96 channels.
                  o2::aod::ft0::AmplitudeC,       // amplitudeC 	std::vector<float> 	Amplitudes of non-zero channels on the C-side. The channel IDs are given in ChannelC (at the same index)
                  o2::aod::ft0::ChannelC,         // channelC 	std::vector<uint8_t> 	Channel IDs on the C side which had non-zero amplitudes. There are at maximum 112 channels.
                  o2::aod::ft0::TimeA,            // timeA 	float 	Average A-side time
                  o2::aod::ft0::TimeC,            // timeC 	float 	Average C-side time
                  o2::aod::ft0::TriggerMask,      // triggerMask 	uint8_t
                  o2::aod::drcolmn::PosZ,         // D 	posZ     	float 	Z position calculated from timeA and timeC in cm
                  o2::aod::drcolmn::CollTime,     // D 	collTime 	float 	Collision time, one need also check validation (code below) for timeA and timeC
                  o2::aod::drcolmn::IsValidTimeA, // D 	isValidTimeA 	bool 	Checks if time from A side was calculated, and if is not dummy
                  o2::aod::drcolmn::IsValidTimeC, // D 	isValidTimeC 	bool 	Checks if time from C side was calculated
                  o2::aod::drcolmn::IsValidTime,  // D 	isValidTime 	bool 	Checks if times from A and C side were calculated simultaneously
                  o2::aod::drcolmn::SumAmpA,      // D 	sumAmpA 	float 	Calculates sum of positive amplitudes from side A
                  o2::aod::drcolmn::SumAmpC,      // D 	sumAmpC 	float 	Calculates sum of positive amplitudes from side C
                  o2::aod::timestamp::Timestamp,
                  o2::aod::drcolmn::TfId,
                  o2::aod::drcolmn::BCinTF);

// Defining Track Tables
namespace drTrack
{
DECLARE_SOA_INDEX_COLUMN(DrCollision, drCollision);           // To index it to derived collisions
DECLARE_SOA_COLUMN(IsWithinBeamPipe, isWithinBeamPipe, bool); // D 	isWithinBeamPipe 	bool 	Is the track within the beam pipe (= successfully propagated to a collision vertex)
DECLARE_SOA_COLUMN(Px, px, float);                            // D 	px 	                float 	Momentum in x-direction in GeV/c
DECLARE_SOA_COLUMN(Py, py, float);                            // D 	py 	                float 	Momentum in y-direction in GeV/c
DECLARE_SOA_COLUMN(Pz, pz, float);                            // D 	pz 	                float 	Momentum in z-direction in GeV/c
// DECLARE_SOA_COLUMN(PVector         ,pVector         ,std::array<float,3>);// D 	pVector             std::array<float,3> 	Momentum vector in x,y,z-directions in GeV/c
DECLARE_SOA_COLUMN(Energy, energy, float);     // D 	energy 	            float 	Track energy, computed under the mass assumption given as input
DECLARE_SOA_COLUMN(Rapidity, rapidity, float); // D 	rapidity 	        float 	Track rapidity, computed under the mass assumption given as input
DECLARE_SOA_COLUMN(Sign, sign, short);         // D 	sign 	            short 	Charge: positive: 1, negative: -1

// Tracks Extra
DECLARE_SOA_COLUMN(TPCDeltaTFwd, tpcDeltaTFwd, float);        // D 	tpcDeltaTFwd 	                    float 	Delta Forward of track time in TPC time bis
DECLARE_SOA_COLUMN(TPCDeltaTBwd, tpcDeltaTBwd, float);        // D 	tpcDeltaTBwd 	                    float 	Delta Backward of track time in TPC time bis
DECLARE_SOA_COLUMN(PIDForTracking, pidForTracking, uint32_t); // D 	pidForTracking 	                    uint32_t 	PID hypothesis used during tracking. See the constants in the class PID in PID.h
DECLARE_SOA_COLUMN(IsPVContributor, isPVContributor, bool);   // D 	isPVContributor 	                bool 	Run 3: Has this track contributed to the collision vertex fit

DECLARE_SOA_COLUMN(DetectorMask, detectorMask, uint8_t); // D 	hasITS 	                            bool 	Flag to check if track has a ITS match

DECLARE_SOA_COLUMN(HasITS, hasITS, bool); // D 	hasITS 	                            bool 	Flag to check if track has a ITS match
DECLARE_SOA_COLUMN(HasTPC, hasTPC, bool); // D 	hasTPC 	                            bool 	Flag to check if track has a TPC match
DECLARE_SOA_COLUMN(HasTRD, hasTRD, bool); // D 	hasTRD 	                            bool 	Flag to check if track has a TRD match
DECLARE_SOA_COLUMN(HasTOF, hasTOF, bool); // D 	hasTOF 	                            bool 	Flag to check if track has a TOF measurement

DECLARE_SOA_COLUMN(TPCNClsFound, tpcNClsFound, int16_t);             // D 	tpcNClsFound 	                    int16_t 	Number of found TPC clusters
DECLARE_SOA_COLUMN(TPCNClsCrossedRows, tpcNClsCrossedRows, int16_t); // D 	tpcNClsCrossedRows 	                int16_t 	Number of crossed TPC Rows
DECLARE_SOA_COLUMN(ITSClusterMap, itsClusterMap, uint8_t);           // D 	itsClusterMap 	                    uint8_t 	ITS cluster map, one bit per a layer, starting from the innermost
DECLARE_SOA_COLUMN(ITSNCls, itsNCls, uint8_t);                       // D 	itsNCls 	                        uint8_t 	Number of ITS clusters
DECLARE_SOA_COLUMN(ITSNClsInnerBarrel, itsNClsInnerBarrel, uint8_t); // D 	itsNClsInnerBarrel 	                uint8_t 	Number of ITS clusters in the Inner Barrel
DECLARE_SOA_COLUMN(ITSClsSizeInLayer, itsClsSizeInLayer, uint8_t);   // D 	itsClsSizeInLayer 	                uint8_t 	Size of the ITS cluster in a given layer
DECLARE_SOA_COLUMN(IsITSAfterburner, isITSAfterburner, bool);        // D 	isITSAfterburner 	                bool 	If the track used the afterburner in the ITS
// DECLARE_SOA_COLUMN(TOFExpTime 	                 ,tofExpTime 	               ,float   );//D 	tofExpTime 	                        float 	Expected time for the track to reach the TOF
// DECLARE_SOA_COLUMN(TOFExpTimePi 	             ,tofExpTimePi 	               ,float   );//D 	tofExpTimePi 	                    float 	Expected time for the track to reach the TOF under the pion hypothesis
// DECLARE_SOA_COLUMN(TOFExpTimeKa 	             ,tofExpTimeKa 	               ,float   );//D 	tofExpTimeKa 	                    float 	Expected time for the track to reach the TOF under the kaon hypothesis
// DECLARE_SOA_COLUMN(TOFExpTimePr 	             ,tofExpTimePr 	               ,float   );//D 	tofExpTimePr 	                    float 	Expected time for the track to reach the TOF under the proton hypothesis
DECLARE_SOA_COLUMN(TPCCrossedRowsOverFindableCls, tpcCrossedRowsOverFindableCls, float); // D 	tpcCrossedRowsOverFindableCls 	    float 	Ratio crossed rows over findable clusters
DECLARE_SOA_COLUMN(TPCFoundOverFindableCls, tpcFoundOverFindableCls, float);             // D 	tpcFoundOverFindableCls             float 	Ratio of found over findable clusters
DECLARE_SOA_COLUMN(TPCFractionSharedCls, tpcFractionSharedCls, float);                   // D 	tpcFractionSharedCls 	            float 	Fraction of shared TPC clusters
DECLARE_SOA_COLUMN(DetectorMap, detectorMap, uint8_t);                                   // E 	detectorMap 	                    uint8_t 	Detector map version 1, see enum DetectorMapEnum

DECLARE_SOA_COLUMN(Beta, beta, float);
DECLARE_SOA_COLUMN(Mass, mass, float);

DECLARE_SOA_COLUMN(TPCNSigmaPi, tpcNSigmaPi, float);
DECLARE_SOA_COLUMN(TOFNSigmaPi, tofNSigmaPi, float);
DECLARE_SOA_COLUMN(TPCNSigmaKa, tpcNSigmaKa, float);
DECLARE_SOA_COLUMN(TOFNSigmaKa, tofNSigmaKa, float);
DECLARE_SOA_COLUMN(TPCNSigmaPr, tpcNSigmaPr, float);
DECLARE_SOA_COLUMN(TOFNSigmaPr, tofNSigmaPr, float);
DECLARE_SOA_COLUMN(TPCNSigmaEl, tpcNSigmaEl, float);
DECLARE_SOA_COLUMN(TOFNSigmaEl, tofNSigmaEl, float);
DECLARE_SOA_COLUMN(TPCNSigmaDe, tpcNSigmaDe, float);
DECLARE_SOA_COLUMN(TOFNSigmaDe, tofNSigmaDe, float);

DECLARE_SOA_COLUMN(CTPCNSigmaPi, ctpcNSigmaPi, uint8_t);
DECLARE_SOA_COLUMN(CTOFNSigmaPi, ctofNSigmaPi, uint8_t);
DECLARE_SOA_COLUMN(CTPCNSigmaKa, ctpcNSigmaKa, uint8_t);
DECLARE_SOA_COLUMN(CTOFNSigmaKa, ctofNSigmaKa, uint8_t);
DECLARE_SOA_COLUMN(CTPCNSigmaPr, ctpcNSigmaPr, uint8_t);
DECLARE_SOA_COLUMN(CTOFNSigmaPr, ctofNSigmaPr, uint8_t);
DECLARE_SOA_COLUMN(CTPCNSigmaEl, ctpcNSigmaEl, uint8_t);
DECLARE_SOA_COLUMN(CTOFNSigmaEl, ctofNSigmaEl, uint8_t);
DECLARE_SOA_COLUMN(CTPCNSigmaDe, ctpcNSigmaDe, uint8_t);
DECLARE_SOA_COLUMN(CTOFNSigmaDe, ctofNSigmaDe, uint8_t);

DECLARE_SOA_COLUMN(TOFExpTimeEl, tofExpTimeEl, float); // D	tofExpTimeEl	float	Expected time for the track to reach the TOF under the electron hypothesis
DECLARE_SOA_COLUMN(TOFExpTimeMu, tofExpTimeMu, float); // D	tofExpTimeMu	float	Expected time for the track to reach the TOF under the muon hypothesis
DECLARE_SOA_COLUMN(TOFExpTimePi, tofExpTimePi, float); // D	tofExpTimePi	float	Expected time for the track to reach the TOF under the pion hypothesis
DECLARE_SOA_COLUMN(TOFExpTimeKa, tofExpTimeKa, float); // D	tofExpTimeKa	float	Expected time for the track to reach the TOF under the kaon hypothesis
DECLARE_SOA_COLUMN(TOFExpTimePr, tofExpTimePr, float); // D	tofExpTimePr	float	Expected time for the track to reach the TOF under the proton hypothesis
DECLARE_SOA_COLUMN(TOFExpTimeDe, tofExpTimeDe, float); // D	tofExpTimeDe	float	Expected time for the track to reach the TOF under the deuteron hypothesis
DECLARE_SOA_COLUMN(TOFExpTimeTr, tofExpTimeTr, float); // D	tofExpTimeTr	float	Expected time for the track to reach the TOF under the triton hypothesis
DECLARE_SOA_COLUMN(TOFExpTimeHe, tofExpTimeHe, float); // D	tofExpTimeHe	float	Expected time for the track to reach the TOF under the helium3 hypothesis
DECLARE_SOA_COLUMN(TOFExpTimeAl, tofExpTimeAl, float); // D	tofExpTimeAl

DECLARE_SOA_COLUMN(ITSNSigmaEl, itsNSigmaEl, float);
DECLARE_SOA_COLUMN(ITSNSigmaMu, itsNSigmaMu, float);
DECLARE_SOA_COLUMN(ITSNSigmaPi, itsNSigmaPi, float);
DECLARE_SOA_COLUMN(ITSNSigmaKa, itsNSigmaKa, float);
DECLARE_SOA_COLUMN(ITSNSigmaPr, itsNSigmaPr, float);
DECLARE_SOA_COLUMN(ITSNSigmaDe, itsNSigmaDe, float);
DECLARE_SOA_COLUMN(ITSNSigmaTr, itsNSigmaTr, float);
DECLARE_SOA_COLUMN(ITSNSigmaHe, itsNSigmaHe, float);
DECLARE_SOA_COLUMN(ITSNSigmaAl, itsNSigmaAl, float);

DECLARE_SOA_COLUMN(CITSNSigmaEl, citsNSigmaEl, uint8_t);
DECLARE_SOA_COLUMN(CITSNSigmaMu, citsNSigmaMu, uint8_t);
DECLARE_SOA_COLUMN(CITSNSigmaPi, citsNSigmaPi, uint8_t);
DECLARE_SOA_COLUMN(CITSNSigmaKa, citsNSigmaKa, uint8_t);
DECLARE_SOA_COLUMN(CITSNSigmaPr, citsNSigmaPr, uint8_t);
DECLARE_SOA_COLUMN(CITSNSigmaDe, citsNSigmaDe, uint8_t);
DECLARE_SOA_COLUMN(CITSNSigmaTr, citsNSigmaTr, uint8_t);
DECLARE_SOA_COLUMN(CITSNSigmaHe, citsNSigmaHe, uint8_t);
DECLARE_SOA_COLUMN(CITSNSigmaAl, citsNSigmaAl, uint8_t);

DECLARE_SOA_COLUMN(EventCollisionTime, eventCollisionTime, float);
DECLARE_SOA_COLUMN(IsEvTimeDefined, isEvTimeDefined, bool); //  D	isEvTimeDefined	bool	True if the Event Time was computed with any method i.e. there is a usable event time
DECLARE_SOA_COLUMN(IsEvTimeTOF, isEvTimeTOF, bool);         //  D	isEvTimeTOF	bool	True if the Event Time was computed with the TOF
DECLARE_SOA_COLUMN(IsEvTimeT0AC, isEvTimeT0AC, bool);       //  D	isEvTimeT0AC	bool	True if the Event Time was computed with the T0AC
DECLARE_SOA_COLUMN(IsEvTimeTOFT0AC, isEvTimeTOFT0AC, bool); //

DECLARE_SOA_COLUMN(WeightDS, weightDS, float);
DECLARE_SOA_COLUMN(WeightDSPt, weightDSPt, float);
DECLARE_SOA_COLUMN(WeightDSQPt, weightDSQPt, float);

DECLARE_SOA_COLUMN(TriggerMaskDS, triggerMaskDS, int);

DECLARE_SOA_COLUMN(MotherV0GIList, motherV0GIList, std::vector<int64_t>);
DECLARE_SOA_COLUMN(MotherV0OIList, motherV0OIList, std::vector<int64_t>);
DECLARE_SOA_COLUMN(MotherD0GIList, motherD0GIList, std::vector<int64_t>);
DECLARE_SOA_COLUMN(MotherD0OIList, motherD0OIList, std::vector<int64_t>);
} // namespace drTrack

namespace trackmeanocc
{
// Lossy Compression
DECLARE_SOA_COLUMN(CmoPrimUnfm80, cmoPrimUnfm80, uint8_t);
DECLARE_SOA_COLUMN(CmoFV0AUnfm80, cmoFV0AUnfm80, uint8_t);
DECLARE_SOA_COLUMN(CmoFV0CUnfm80, cmoFV0CUnfm80, uint8_t);
DECLARE_SOA_COLUMN(CmoFT0AUnfm80, cmoFT0AUnfm80, uint8_t);
DECLARE_SOA_COLUMN(CmoFT0CUnfm80, cmoFT0CUnfm80, uint8_t);

DECLARE_SOA_COLUMN(CmoFDDAUnfm80, cmoFDDAUnfm80, uint8_t);
DECLARE_SOA_COLUMN(CmoFDDCUnfm80, cmoFDDCUnfm80, uint8_t);

DECLARE_SOA_COLUMN(CmoNTrackITSUnfm80, cmoNTrackITSUnfm80, uint8_t);
DECLARE_SOA_COLUMN(CmoNTrackTPCUnfm80, cmoNTrackTPCUnfm80, uint8_t);
DECLARE_SOA_COLUMN(CmoNTrackTRDUnfm80, cmoNTrackTRDUnfm80, uint8_t);
DECLARE_SOA_COLUMN(CmoNTrackTOFUnfm80, cmoNTrackTOFUnfm80, uint8_t);
DECLARE_SOA_COLUMN(CmoNTrackSizeUnfm80, cmoNTrackSizeUnfm80, uint8_t);
DECLARE_SOA_COLUMN(CmoNTrackTPCAUnfm80, cmoNTrackTPCAUnfm80, uint8_t);
DECLARE_SOA_COLUMN(CmoNTrackTPCCUnfm80, cmoNTrackTPCCUnfm80, uint8_t);
DECLARE_SOA_COLUMN(CmoNTrackITSTPCUnfm80, cmoNTrackITSTPCUnfm80, uint8_t);
DECLARE_SOA_COLUMN(CmoNTrackITSTPCAUnfm80, cmoNTrackITSTPCAUnfm80, uint8_t);
DECLARE_SOA_COLUMN(CmoNTrackITSTPCCUnfm80, cmoNTrackITSTPCCUnfm80, uint8_t);

DECLARE_SOA_COLUMN(CmoMultNTracksHasITSUnfm80, cmoMultNTracksHasITSUnfm80, uint8_t);
DECLARE_SOA_COLUMN(CmoMultNTracksHasTPCUnfm80, cmoMultNTracksHasTPCUnfm80, uint8_t);
DECLARE_SOA_COLUMN(CmoMultNTracksHasTOFUnfm80, cmoMultNTracksHasTOFUnfm80, uint8_t);
DECLARE_SOA_COLUMN(CmoMultNTracksHasTRDUnfm80, cmoMultNTracksHasTRDUnfm80, uint8_t);
DECLARE_SOA_COLUMN(CmoMultNTracksITSOnlyUnfm80, cmoMultNTracksITSOnlyUnfm80, uint8_t);
DECLARE_SOA_COLUMN(CmoMultNTracksTPCOnlyUnfm80, cmoMultNTracksTPCOnlyUnfm80, uint8_t);
DECLARE_SOA_COLUMN(CmoMultNTracksITSTPCUnfm80, cmoMultNTracksITSTPCUnfm80, uint8_t);
DECLARE_SOA_COLUMN(CmoMultAllTracksTPCOnlyUnfm80, cmoMultAllTracksTPCOnlyUnfm80, uint8_t);

DECLARE_SOA_COLUMN(CmoRobustT0V0PrimUnfm80, cmoRobustT0V0PrimUnfm80, uint8_t);
DECLARE_SOA_COLUMN(CmoRobustFDDT0V0PrimUnfm80, cmoRobustFDDT0V0PrimUnfm80, uint8_t);
DECLARE_SOA_COLUMN(CmoRobustNtrackDetUnfm80, cmoRobustNtrackDetUnfm80, uint8_t);
DECLARE_SOA_COLUMN(CmoRobustMultExtraTableUnfm80, cmoRobustMultExtraTableUnfm80, uint8_t);

DECLARE_SOA_COLUMN(CwmoPrimUnfm80, cwmoPrimUnfm80, uint8_t);
DECLARE_SOA_COLUMN(CwmoFV0AUnfm80, cwmoFV0AUnfm80, uint8_t);
DECLARE_SOA_COLUMN(CwmoFV0CUnfm80, cwmoFV0CUnfm80, uint8_t);
DECLARE_SOA_COLUMN(CwmoFT0AUnfm80, cwmoFT0AUnfm80, uint8_t);
DECLARE_SOA_COLUMN(CwmoFT0CUnfm80, cwmoFT0CUnfm80, uint8_t);

DECLARE_SOA_COLUMN(CwmoFDDAUnfm80, cwmoFDDAUnfm80, uint8_t);
DECLARE_SOA_COLUMN(CwmoFDDCUnfm80, cwmoFDDCUnfm80, uint8_t);

DECLARE_SOA_COLUMN(CwmoNTrackITSUnfm80, cwmoNTrackITSUnfm80, uint8_t);
DECLARE_SOA_COLUMN(CwmoNTrackTPCUnfm80, cwmoNTrackTPCUnfm80, uint8_t);
DECLARE_SOA_COLUMN(CwmoNTrackTRDUnfm80, cwmoNTrackTRDUnfm80, uint8_t);
DECLARE_SOA_COLUMN(CwmoNTrackTOFUnfm80, cwmoNTrackTOFUnfm80, uint8_t);
DECLARE_SOA_COLUMN(CwmoNTrackSizeUnfm80, cwmoNTrackSizeUnfm80, uint8_t);
DECLARE_SOA_COLUMN(CwmoNTrackTPCAUnfm80, cwmoNTrackTPCAUnfm80, uint8_t);
DECLARE_SOA_COLUMN(CwmoNTrackTPCCUnfm80, cwmoNTrackTPCCUnfm80, uint8_t);
DECLARE_SOA_COLUMN(CwmoNTrackITSTPCUnfm80, cwmoNTrackITSTPCUnfm80, uint8_t);
DECLARE_SOA_COLUMN(CwmoNTrackITSTPCAUnfm80, cwmoNTrackITSTPCAUnfm80, uint8_t);
DECLARE_SOA_COLUMN(CwmoNTrackITSTPCCUnfm80, cwmoNTrackITSTPCCUnfm80, uint8_t);

DECLARE_SOA_COLUMN(CwmoMultNTracksHasITSUnfm80, cwmoMultNTracksHasITSUnfm80, uint8_t);
DECLARE_SOA_COLUMN(CwmoMultNTracksHasTPCUnfm80, cwmoMultNTracksHasTPCUnfm80, uint8_t);
DECLARE_SOA_COLUMN(CwmoMultNTracksHasTOFUnfm80, cwmoMultNTracksHasTOFUnfm80, uint8_t);
DECLARE_SOA_COLUMN(CwmoMultNTracksHasTRDUnfm80, cwmoMultNTracksHasTRDUnfm80, uint8_t);
DECLARE_SOA_COLUMN(CwmoMultNTracksITSOnlyUnfm80, cwmoMultNTracksITSOnlyUnfm80, uint8_t);
DECLARE_SOA_COLUMN(CwmoMultNTracksTPCOnlyUnfm80, cwmoMultNTracksTPCOnlyUnfm80, uint8_t);
DECLARE_SOA_COLUMN(CwmoMultNTracksITSTPCUnfm80, cwmoMultNTracksITSTPCUnfm80, uint8_t);
DECLARE_SOA_COLUMN(CwmoMultAllTracksTPCOnlyUnfm80, cwmoMultAllTracksTPCOnlyUnfm80, uint8_t);

DECLARE_SOA_COLUMN(CwmoRobustT0V0PrimUnfm80, cwmoRobustT0V0PrimUnfm80, uint8_t);
DECLARE_SOA_COLUMN(CwmoRobustFDDT0V0PrimUnfm80, cwmoRobustFDDT0V0PrimUnfm80, uint8_t);
DECLARE_SOA_COLUMN(CwmoRobustNtrackDetUnfm80, cwmoRobustNtrackDetUnfm80, uint8_t);
DECLARE_SOA_COLUMN(CwmoRobustMultExtraTableUnfm80, cwmoRobustMultExtraTableUnfm80, uint8_t);

} // namespace trackmeanocc

// // } // namespace drTrack
DECLARE_SOA_TABLE(DrTracksMothers, "AOD", "DRTRACKSMOTHER", o2::soa::Index<>,
                  // o2::aod::drcolmn::GlobalIndex,
                  // o2::aod::drcolmn::OriginalIndex,
                  // o2::aod::drcolmn::GBCid,
                  // o2::aod::drcolmn::GCollId,
                  drTrack::MotherV0GIList,
                  drTrack::MotherV0OIList,
                  drTrack::MotherD0GIList,
                  drTrack::MotherD0OIList);

DECLARE_SOA_TABLE(DrTracks0, "AOD", "DRTRACK0", o2::soa::Index<>,
                  drTrack::TPCNSigmaPi,
                  drTrack::TOFNSigmaPi,
                  drTrack::TPCNSigmaKa,
                  drTrack::TOFNSigmaKa,
                  drTrack::TPCNSigmaPr,
                  drTrack::TOFNSigmaPr,
                  drTrack::TPCNSigmaEl,
                  drTrack::TOFNSigmaEl,
                  drTrack::TPCNSigmaDe,
                  drTrack::TOFNSigmaDe,

                  o2::aod::trackmeanocc::TmoPrimUnfm80,
                  o2::aod::trackmeanocc::TmoFV0AUnfm80,
                  o2::aod::trackmeanocc::TmoFV0CUnfm80,
                  o2::aod::trackmeanocc::TmoFT0AUnfm80,
                  o2::aod::trackmeanocc::TmoFT0CUnfm80,
                  o2::aod::trackmeanocc::TmoFDDAUnfm80,
                  o2::aod::trackmeanocc::TmoFDDCUnfm80,

                  o2::aod::trackmeanocc::TmoNTrackITSUnfm80,
                  o2::aod::trackmeanocc::TmoNTrackTPCUnfm80,
                  o2::aod::trackmeanocc::TmoNTrackTRDUnfm80,
                  o2::aod::trackmeanocc::TmoNTrackTOFUnfm80,
                  o2::aod::trackmeanocc::TmoNTrackSizeUnfm80,
                  o2::aod::trackmeanocc::TmoNTrackTPCAUnfm80,
                  o2::aod::trackmeanocc::TmoNTrackTPCCUnfm80,
                  o2::aod::trackmeanocc::TmoNTrackITSTPCUnfm80,
                  o2::aod::trackmeanocc::TmoNTrackITSTPCAUnfm80,
                  o2::aod::trackmeanocc::TmoNTrackITSTPCCUnfm80,

                  o2::aod::trackmeanocc::TmoMultNTracksHasITSUnfm80,
                  o2::aod::trackmeanocc::TmoMultNTracksHasTPCUnfm80,
                  o2::aod::trackmeanocc::TmoMultNTracksHasTOFUnfm80,
                  o2::aod::trackmeanocc::TmoMultNTracksHasTRDUnfm80,
                  o2::aod::trackmeanocc::TmoMultNTracksITSOnlyUnfm80,
                  o2::aod::trackmeanocc::TmoMultNTracksTPCOnlyUnfm80,
                  o2::aod::trackmeanocc::TmoMultNTracksITSTPCUnfm80,
                  o2::aod::trackmeanocc::TmoMultAllTracksTPCOnlyUnfm80,

                  o2::aod::trackmeanocc::TmoRobustT0V0PrimUnfm80,
                  o2::aod::trackmeanocc::TmoRobustFDDT0V0PrimUnfm80,
                  o2::aod::trackmeanocc::TmoRobustNtrackDetUnfm80,
                  o2::aod::trackmeanocc::TmoRobustMultExtraTableUnfm80,

                  o2::aod::trackmeanocc::TwmoPrimUnfm80,
                  o2::aod::trackmeanocc::TwmoFV0AUnfm80,
                  o2::aod::trackmeanocc::TwmoFV0CUnfm80,
                  o2::aod::trackmeanocc::TwmoFT0AUnfm80,
                  o2::aod::trackmeanocc::TwmoFT0CUnfm80,
                  o2::aod::trackmeanocc::TwmoFDDAUnfm80,
                  o2::aod::trackmeanocc::TwmoFDDCUnfm80,

                  o2::aod::trackmeanocc::TwmoNTrackITSUnfm80,
                  o2::aod::trackmeanocc::TwmoNTrackTPCUnfm80,
                  o2::aod::trackmeanocc::TwmoNTrackTRDUnfm80,
                  o2::aod::trackmeanocc::TwmoNTrackTOFUnfm80,
                  o2::aod::trackmeanocc::TwmoNTrackSizeUnfm80,
                  o2::aod::trackmeanocc::TwmoNTrackTPCAUnfm80,
                  o2::aod::trackmeanocc::TwmoNTrackTPCCUnfm80,
                  o2::aod::trackmeanocc::TwmoNTrackITSTPCUnfm80,
                  o2::aod::trackmeanocc::TwmoNTrackITSTPCAUnfm80,
                  o2::aod::trackmeanocc::TwmoNTrackITSTPCCUnfm80,

                  o2::aod::trackmeanocc::TwmoMultNTracksHasITSUnfm80,
                  o2::aod::trackmeanocc::TwmoMultNTracksHasTPCUnfm80,
                  o2::aod::trackmeanocc::TwmoMultNTracksHasTOFUnfm80,
                  o2::aod::trackmeanocc::TwmoMultNTracksHasTRDUnfm80,
                  o2::aod::trackmeanocc::TwmoMultNTracksITSOnlyUnfm80,
                  o2::aod::trackmeanocc::TwmoMultNTracksTPCOnlyUnfm80,
                  o2::aod::trackmeanocc::TwmoMultNTracksITSTPCUnfm80,
                  o2::aod::trackmeanocc::TwmoMultAllTracksTPCOnlyUnfm80,

                  o2::aod::trackmeanocc::TwmoRobustT0V0PrimUnfm80,
                  o2::aod::trackmeanocc::TwmoRobustFDDT0V0PrimUnfm80,
                  o2::aod::trackmeanocc::TwmoRobustNtrackDetUnfm80,
                  o2::aod::trackmeanocc::TwmoRobustMultExtraTableUnfm80,

                  drTrack::ITSNSigmaEl,
                  drTrack::ITSNSigmaMu,
                  drTrack::ITSNSigmaPi,
                  drTrack::ITSNSigmaKa,
                  drTrack::ITSNSigmaPr,
                  drTrack::ITSNSigmaDe,
                  drTrack::ITSNSigmaTr,
                  drTrack::ITSNSigmaHe,
                  drTrack::ITSNSigmaAl);

DECLARE_SOA_TABLE(DrTracks, "AOD", "DRTRACK", o2::soa::Index<>,
                  o2::aod::drcolmn::GlobalIndex,
                  o2::aod::drcolmn::OriginalIndex,
                  o2::aod::drcolmn::GBCid,
                  o2::aod::drcolmn::GCollId,
                  o2::aod::drcolmn::GDFId,

                  drTrack::DrCollisionId, // Index Column
                  // Tracks
                  o2::aod::track::CollisionId, //      collisionId 	    int32 	Collision to which this track belongs
                  o2::aod::track::TrackType,   //      trackType 	        uint8_t 	Type of track. See enum TrackTypeEnum. This cannot be used to decide which detector has contributed to this track. Use hasITShasTPCetc.
                  o2::aod::track::X,           //      x 	                float
                  o2::aod::track::Alpha,       //      alpha 	            float
                  o2::aod::track::Y,           //      y 	                float
                  o2::aod::track::Z,           //      z 	                float
                  o2::aod::track::Snp,         //      snp 	            float
                  o2::aod::track::Tgl,         //      tgl 	            float
                  o2::aod::track::Signed1Pt,   //      signed1Pt 	        float 	(sign of charge)/Pt in c/GeV. Use pt and sign instead
                  // o2::aod::track::IsWithinBeamPipe, // D 	isWithinBeamPipe 	bool 	Is the track within the beam pipe (= successfully propagated to a collision vertex)
                  // o2::aod::track::Px,               // D 	px 	                float 	Momentum in x-direction in GeV/c
                  // o2::aod::track::Py,               // D 	py 	                float 	Momentum in y-direction in GeV/c
                  // o2::aod::track::Pz,               // D 	pz 	                float 	Momentum in z-direction in GeV/c
                  // o2::aod::track::PVector,       //    D 	pVector             std::array<float3> 	Momentum vector in xyz-directions in GeV/c
                  // o2::aod::track::Energy,        //    D 	energy 	            float 	Track energycomputed under the mass assumption given as input
                  // o2::aod::track::Rapidity,      //    D 	rapidity 	        float 	Track rapiditycomputed under the mass assumption given as input
                  // o2::aod::track::Sign,             // D 	sign 	            short 	Charge: positive: 1negative: -1
                  o2::aod::track::Pt,  // E 	pt 	                float 	Transverse momentum of the track in GeV/c
                  o2::aod::track::P,   // E 	p 	                float 	Momentum in Gev/c
                  o2::aod::track::Eta, // E 	eta 	            float 	Pseudorapidity
                  o2::aod::track::Phi, // E 	phi 	            float 	Phi of the trackin radians within [02pi)

                  drTrack::IsWithinBeamPipe, // D 	isWithinBeamPipe 	bool 	Is the track within the beam pipe (= successfully propagated to a collision vertex)
                  drTrack::Px,               // D 	px 	                float 	Momentum in x-direction in GeV/c
                  drTrack::Py,               // D 	py 	                float 	Momentum in y-direction in GeV/c
                  drTrack::Pz,               // D 	pz 	                float 	Momentum in z-direction in GeV/c
                  // drTrack::PVector,                 // D 	pVector             std::array<float3> 	Momentum vector in xyz-directions in GeV/c
                  drTrack::Energy,   // D 	energy 	            float 	Track energycomputed under the mass assumption given as input
                  drTrack::Rapidity, // D 	rapidity 	        float 	Track rapiditycomputed under the mass assumption given as input
                  drTrack::Sign,     // D 	sign 	            short 	Charge: positive: 1negative: -1

                  // // // TracksCov
                  // o2::aod::track::SigmaY, 	//	sigmaY 	    float 	Covariance matrix
                  // o2::aod::track::SigmaZ, 	//	sigmaZ 	    float 	Covariance matrix
                  // o2::aod::track::SigmaSnp, 	//	sigmaSnp 	float 	Covariance matrix
                  // o2::aod::track::SigmaTgl, 	//	sigmaTgl 	float 	Covariance matrix
                  // o2::aod::track::Sigma1Pt, 	//	sigma1Pt 	float 	Covariance matrix
                  // o2::aod::track::RhoZY, 		//    rhoZY 	    int8_t 	Covariance matrix in compressed form
                  // o2::aod::track::RhoSnpY, 	//	rhoSnpY 	int8_t 	Covariance matrix in compressed form
                  // o2::aod::track::RhoSnpZ, 	//	rhoSnpZ 	int8_t 	Covariance matrix in compressed form
                  // o2::aod::track::RhoTglY, 	//	rhoTglY 	int8_t 	Covariance matrix in compressed form
                  // o2::aod::track::RhoTglZ, 	//	rhoTglZ 	int8_t 	Covariance matrix in compressed form
                  // o2::aod::track::RhoTglSnp, 	//	rhoTglSnp 	int8_t 	Covariance matrix in compressed form
                  // o2::aod::track::Rho1PtY, 	//	rho1PtY 	int8_t 	Covariance matrix in compressed form
                  // o2::aod::track::Rho1PtZ, 	//	rho1PtZ 	int8_t 	Covariance matrix in compressed form
                  // o2::aod::track::Rho1PtSnp, 	//	rho1PtSnp 	int8_t 	Covariance matrix in compressed form
                  // o2::aod::track::Rho1PtTgl, 	//	rho1PtTgl 	int8_t 	Covariance matrix in compressed form
                  // o2::aod::track::CYY, 	    //E 	cYY 	    float
                  // o2::aod::track::CZY, 	    //E 	cZY 	    float
                  // o2::aod::track::CZZ, 	    //E 	cZZ 	    float
                  // o2::aod::track::CSnpY, 	    //E 	cSnpY 	    float
                  // o2::aod::track::CSnpZ, 	    //E 	cSnpZ 	    float
                  // o2::aod::track::CSnpSnp, 	//E 	cSnpSnp 	float
                  // o2::aod::track::CTglY, 	    //E 	cTglY 	    float
                  // o2::aod::track::CTglZ, 	    //E 	cTglZ 	    float
                  // o2::aod::track::CTglSnp, 	//E 	cTglSnp 	float
                  // o2::aod::track::CTglTgl, 	//E 	cTglTgl 	float
                  // o2::aod::track::C1PtY, 	    //E 	c1PtY 	    float
                  // o2::aod::track::C1PtZ, 	    //E 	c1PtZ 	    float
                  // o2::aod::track::C1PtSnp, 	//E 	c1PtSnp 	float
                  // o2::aod::track::C1PtTgl, 	//E 	c1PtTgl 	float
                  // o2::aod::track::C1Pt21Pt2, 	//E 	c1Pt21Pt2 	float

                  // TracksExtra
                  o2::aod::track::TPCInnerParam,                   //    tpcInnerParam 	                    float 	Momentum at inner wall of the TPC
                  o2::aod::track::Flags,                           //    flags 	                            uint32_t 	Track flags. Run 2: see TrackFlagsRun2Enum | Run 3: see TrackFlags
                  o2::aod::track::ITSClusterSizes,                 //    itsClusterSizes 	                uint32_t 	Clusters sizesfour bits per a layerstarting from the innermost
                  o2::aod::track::TPCNClsFindable,                 //    tpcNClsFindable 	                uint8_t 	Findable TPC clusters for this track geometry
                  o2::aod::track::TPCNClsFindableMinusFound,       //    tpcNClsFindableMinusFound 	        int8_t 	TPC Clusters: Findable - Found
                  o2::aod::track::TPCNClsFindableMinusCrossedRows, //    tpcNClsFindableMinusCrossedRows 	int8_t 	TPC Clusters: Findable - crossed rows
                  o2::aod::track::TPCNClsShared,                   //    tpcNClsShared 	                    uint8_t 	Number of shared TPC clusters
                  // o2::aod::track::v001::extensions::TPCDeltaTFwd,  	//D 	tpcDeltaTFwd 	                    float 	Delta Forward of track time in TPC time bis
                  // o2::aod::track::v001::extensions::TPCDeltaTBwd,  	//D 	tpcDeltaTBwd 	                    float 	Delta Backward of track time in TPC time bis
                  o2::aod::track::TRDPattern, //    trdPattern 	                        uint8_t 	Contributor to the track on TRD layer in bits 0-5starting from the innermostbit 6 indicates a potentially split trackletbit 7 if the track crossed a padrow
                  o2::aod::track::ITSChi2NCl, //    itsChi2NCl 	                        float 	Chi2 / cluster for the ITS track segment
                  o2::aod::track::TPCChi2NCl, //    tpcChi2NCl 	                        float 	Chi2 / cluster for the TPC track segment
                  o2::aod::track::TRDChi2,    //    trdChi2 	                        float 	Chi2 for the TRD track segment
                  o2::aod::track::TOFChi2,    //    tofChi2 	                        float 	Chi2 for the TOF track segment
                  o2::aod::track::TPCSignal,  //    tpcSignal 	                        float 	dE/dx signal in the TPC
                  o2::aod::track::TRDSignal,  //    trdSignal 	                        float 	PID signal in the TRD
                  o2::aod::track::Length,     //    length 	                            float 	Track length
                  o2::aod::track::TOFExpMom,  //    tofExpMom 	                        float 	TOF expected momentum obtained in trackingused to compute the expected times
                  // o2::aod::track::PIDForTracking, 	                //D 	pidForTracking 	                    uint32_t 	PID hypothesis used during tracking. See the constants in the class PID in PID.h
                  // o2::aod::track::IsPVContributor, 	            //D 	isPVContributor 	                bool 	Run 3: Has this track contributed to the collision vertex fit
                  // o2::aod::track::HasITS, 	                        //D 	hasITS 	                            bool 	Flag to check if track has a ITS match
                  // o2::aod::track::HasTPC, 	                        //D 	hasTPC 	                            bool 	Flag to check if track has a TPC match
                  // o2::aod::track::HasTRD, 	                        //D 	hasTRD 	                            bool 	Flag to check if track has a TRD match
                  // o2::aod::track::HasTOF, 	                        //D 	hasTOF 	                            bool 	Flag to check if track has a TOF measurement
                  // o2::aod::track::TPCNClsFound, 	                //D 	tpcNClsFound 	                    int16_t 	Number of found TPC clusters
                  // o2::aod::track::TPCNClsCrossedRows, 	            //D 	tpcNClsCrossedRows 	                int16_t 	Number of crossed TPC Rows
                  // o2::aod::track::v001,::ITSClusterMap 	        //D 	itsClusterMap 	                    uint8_t 	ITS cluster mapone bit per a layerstarting from the innermost
                  // o2::aod::track::v001,::ITSNCls 	                //D 	itsNCls 	                        uint8_t 	Number of ITS clusters
                  // o2::aod::track::v001,::ITSNClsInnerBarrel 	    //D 	itsNClsInnerBarrel 	                uint8_t 	Number of ITS clusters in the Inner Barrel
                  // o2::aod::track::v001,::ITSClsSizeInLayer 	    //D 	itsClsSizeInLayer 	                uint8_t 	Size of the ITS cluster in a given layer
                  // o2::aod::track::v001,::IsITSAfterburner 	        //D 	isITSAfterburner 	                bool 	If the track used the afterburner in the ITS
                  // o2::aod::track::TOFExpTime, 	                    //D 	tofExpTime 	                        float 	Expected time for the track to reach the TOF
                  // o2::aod::track::TOFExpTimePi, 	                //D 	tofExpTimePi 	                    float 	Expected time for the track to reach the TOF under the pion hypothesis
                  // o2::aod::track::TOFExpTimeKa, 	                //D 	tofExpTimeKa 	                    float 	Expected time for the track to reach the TOF under the kaon hypothesis
                  // o2::aod::track::TOFExpTimePr, 	                //D 	tofExpTimePr 	                    float 	Expected time for the track to reach the TOF under the proton hypothesis
                  // o2::aod::track::TPCCrossedRowsOverFindableCls, 	//D 	tpcCrossedRowsOverFindableCls 	    float 	Ratio crossed rows over findable clusters
                  // o2::aod::track::TPCFoundOverFindableCls,         //D 	tpcFoundOverFindableCls             float 	Ratio of found over findable clusters
                  // o2::aod::track::TPCFractionSharedCls, 	        //D 	tpcFractionSharedCls 	            float 	Fraction of shared TPC clusters
                  o2::aod::track::TrackEtaEMCAL, //    trackEtaEmcal 	                    float
                  o2::aod::track::TrackPhiEMCAL, //    trackPhiEmcal 	                    float
                  o2::aod::track::TrackTime,     //    trackTime 	                        float 	Estimated time of the track in ns wrt collision.bc or ambiguoustrack.bcSlice[0]
                  o2::aod::track::TrackTimeRes,  //    trackTimeRes 	                    float 	Resolution of the track time in ns (see TrackFlags::TrackTimeResIsRange)
                  // o2::aod::track::v001::DetectorMap, 	            //E 	detectorMap 	                    uint8_t 	Detector map version 1see enum DetectorMapEnum

                  // drTrack::TPCDeltaTFwd,                          //D 	tpcDeltaTFwd 	                    float 	Delta Forward of track time in TPC time bis
                  // drTrack::TPCDeltaTBwd,                          //D 	tpcDeltaTBwd 	                    float 	Delta Backward of track time in TPC time bis
                  drTrack::PIDForTracking, // D 	pidForTracking 	                    uint32_t 	PID hypothesis used during tracking. See the constants in the class PID in PID.h

                  drTrack::IsPVContributor, // D 	isPVContributor 	                bool 	Run 3: Has this track contributed to the collision vertex fit

                  drTrack::DetectorMask,
                  // drTrack::HasITS, 	                        //D 	hasITS 	                            bool 	Flag to check if track has a ITS match
                  // drTrack::HasTPC, 	                        //D 	hasTPC 	                            bool 	Flag to check if track has a TPC match
                  // drTrack::HasTRD, 	                        //D 	hasTRD 	                            bool 	Flag to check if track has a TRD match
                  // drTrack::HasTOF, 	                        //D 	hasTOF 	                            bool 	Flag to check if track has a TOF measurement

                  drTrack::TPCNClsFound,       // D 	tpcNClsFound 	                    int16_t 	Number of found TPC clusters
                  drTrack::TPCNClsCrossedRows, // D 	tpcNClsCrossedRows 	                int16_t 	Number of crossed TPC Rows
                  drTrack::ITSClusterMap,      // D 	itsClusterMap 	                    uint8_t 	ITS cluster mapone bit per a layerstarting from the innermost
                  drTrack::ITSNCls,            // D 	itsNCls 	                        uint8_t 	Number of ITS clusters
                  drTrack::ITSNClsInnerBarrel, // D 	itsNClsInnerBarrel 	                uint8_t 	Number of ITS clusters in the Inner Barrel
                  // drTrack::ITSClsSizeInLayer, 	            //D 	itsClsSizeInLayer 	                uint8_t 	Size of the ITS cluster in a given layer
                  // drTrack::IsITSAfterburner, 	                //D 	isITSAfterburner 	                bool 	If the track used the afterburner in the ITS
                  // drTrack::TOFExpTime, 	                    //D 	tofExpTime 	                        float 	Expected time for the track to reach the TOF
                  // drTrack::TOFExpTimePi, 	                    //D 	tofExpTimePi 	                    float 	Expected time for the track to reach the TOF under the pion hypothesis
                  // drTrack::TOFExpTimeKa, 	                    //D 	tofExpTimeKa 	                    float 	Expected time for the track to reach the TOF under the kaon hypothesis
                  // drTrack::TOFExpTimePr, 	                    //D 	tofExpTimePr 	                    float 	Expected time for the track to reach the TOF under the proton hypothesis
                  drTrack::TPCCrossedRowsOverFindableCls, // D 	tpcCrossedRowsOverFindableCls 	    float 	Ratio crossed rows over findable clusters
                  drTrack::TPCFoundOverFindableCls,       // D 	tpcFoundOverFindableCls             float 	Ratio of found over findable clusters
                  drTrack::TPCFractionSharedCls,          // D 	tpcFractionSharedCls 	            float 	Fraction of shared TPC clusters
                  drTrack::DetectorMap,                   // E 	detectorMap 	                    uint8_t 	Detector map version 1see enum DetectorMapEnum

                  // TracksQA
                  // o2::aod::trackqa::TrackId, 	            //I 	    trackId 	        //int32 	track to which this QA information belongs
                  o2::aod::trackqa::TPCTime0,           //        tpcTime0 	        //float 	tpc only time0 (mTime0 in TPC track)
                  o2::aod::trackqa::TPCDCAR,            //        tpcdcaR 	        //int16_t 	tpc only DCAr
                  o2::aod::trackqa::TPCDCAZ,            //        tpcdcaZ 	        //int16_t 	tpc only DCAz
                  o2::aod::trackqa::TPCClusterByteMask, //        tpcClusterByteMask 	//uint8_t 	tracklet bitmask - track defining 8 tracklets (152=8*19 rows) bit set if nCluster>thr (default 5)
                  o2::aod::trackqa::TPCdEdxMax0R,       //        tpcdEdxMax0R 	    //uint8_t 	TPC dEdxQMax -ROC0/dEdx
                  o2::aod::trackqa::TPCdEdxMax1R,       //        tpcdEdxMax1R 	    //uint8_t 	TPC dEdxQMax -ROC1/dEdx
                  o2::aod::trackqa::TPCdEdxMax2R,       //        tpcdEdxMax2R 	    //uint8_t 	TPC dEdxQMax -ROC2/dEdx
                  o2::aod::trackqa::TPCdEdxMax3R,       //        tpcdEdxMax3R 	    //uint8_t 	TPC dEdxQMax -ROC3/dEdx
                  o2::aod::trackqa::TPCdEdxTot0R,       //        tpcdEdxTot0R 	    //uint8_t 	TPC dEdxQtot -ROC0/dEdx
                  o2::aod::trackqa::TPCdEdxTot1R,       //        tpcdEdxTot1R 	    //uint8_t 	TPC dEdxQtot -ROC1/dEdx
                  o2::aod::trackqa::TPCdEdxTot2R,       //        tpcdEdxTot2R 	    //uint8_t 	TPC dEdxQtot -ROC2/dEdx
                  o2::aod::trackqa::TPCdEdxTot3R,       //        tpcdEdxTot3R 	    //uint8_t 	TPC dEdxQtot -ROC3/dEdx

                  drTrack::Beta, //      beta         float);
                  o2::aod::pidtofbeta::BetaError,
                  drTrack::Mass, //      mass         float);

                  // drTrack::TPCNSigmaPi,
                  // drTrack::TOFNSigmaPi,
                  // drTrack::TPCNSigmaKa,
                  // drTrack::TOFNSigmaKa,
                  // drTrack::TPCNSigmaPr,
                  // drTrack::TOFNSigmaPr,
                  // drTrack::TPCNSigmaEl,
                  // drTrack::TOFNSigmaEl,
                  // drTrack::TPCNSigmaDe,
                  // drTrack::TOFNSigmaDe,

                  drTrack::CTPCNSigmaPi,
                  drTrack::CTOFNSigmaPi,
                  drTrack::CTPCNSigmaKa,
                  drTrack::CTOFNSigmaKa,
                  drTrack::CTPCNSigmaPr,
                  drTrack::CTOFNSigmaPr,
                  drTrack::CTPCNSigmaEl,
                  drTrack::CTOFNSigmaEl,
                  drTrack::CTPCNSigmaDe,
                  drTrack::CTOFNSigmaDe,

                  aod::pidtofsignal::TOFSignal, // tofSignal	float	TOF signal from track time
                  // drTrack::EventCollisionTime,
                  o2::aod::pidflags::GoodTOFMatch,

                  o2::aod::pidflags::TOFFlags, //  tofFlags	uint8_t	Flag for the complementary TOF PID information for the event time
                  drTrack::IsEvTimeDefined,    //  D	isEvTimeDefined	bool	True if the Event Time was computed with any method i.e. there is a usable event time
                  drTrack::IsEvTimeTOF,        //  D	isEvTimeTOF	bool	True if the Event Time was computed with the TOF
                  drTrack::IsEvTimeT0AC,       //  D	isEvTimeT0AC	bool	True if the Event Time was computed with the T0AC
                  drTrack::IsEvTimeTOFT0AC,    //

                  o2::aod::timestamp::Timestamp,
                  o2::aod::drcolmn::TfId,
                  o2::aod::drcolmn::BCinTF,

                  o2::aod::trackmeanocc::TrackId,

                  // o2::aod::trackmeanocc::TmoPrimUnfm80,
                  // o2::aod::trackmeanocc::TmoFV0AUnfm80,
                  // o2::aod::trackmeanocc::TmoFV0CUnfm80,
                  // o2::aod::trackmeanocc::TmoFT0AUnfm80,
                  // o2::aod::trackmeanocc::TmoFT0CUnfm80,
                  // o2::aod::trackmeanocc::TmoFDDAUnfm80,
                  // o2::aod::trackmeanocc::TmoFDDCUnfm80,

                  // o2::aod::trackmeanocc::TmoNTrackITSUnfm80,
                  // o2::aod::trackmeanocc::TmoNTrackTPCUnfm80,
                  // o2::aod::trackmeanocc::TmoNTrackTRDUnfm80,
                  // o2::aod::trackmeanocc::TmoNTrackTOFUnfm80,
                  // o2::aod::trackmeanocc::TmoNTrackSizeUnfm80,
                  // o2::aod::trackmeanocc::TmoNTrackTPCAUnfm80,
                  // o2::aod::trackmeanocc::TmoNTrackTPCCUnfm80,
                  // o2::aod::trackmeanocc::TmoNTrackITSTPCUnfm80,
                  // o2::aod::trackmeanocc::TmoNTrackITSTPCAUnfm80,
                  // o2::aod::trackmeanocc::TmoNTrackITSTPCCUnfm80,

                  // o2::aod::trackmeanocc::TmoMultNTracksHasITSUnfm80,
                  // o2::aod::trackmeanocc::TmoMultNTracksHasTPCUnfm80,
                  // o2::aod::trackmeanocc::TmoMultNTracksHasTOFUnfm80,
                  // o2::aod::trackmeanocc::TmoMultNTracksHasTRDUnfm80,
                  // o2::aod::trackmeanocc::TmoMultNTracksITSOnlyUnfm80,
                  // o2::aod::trackmeanocc::TmoMultNTracksTPCOnlyUnfm80,
                  // o2::aod::trackmeanocc::TmoMultNTracksITSTPCUnfm80,
                  // o2::aod::trackmeanocc::TmoMultAllTracksTPCOnlyUnfm80,

                  // o2::aod::trackmeanocc::TmoRobustT0V0PrimUnfm80,
                  // o2::aod::trackmeanocc::TmoRobustFDDT0V0PrimUnfm80,
                  // o2::aod::trackmeanocc::TmoRobustNtrackDetUnfm80,
                  // o2::aod::trackmeanocc::TmoRobustMultExtraTableUnfm80,

                  // o2::aod::trackmeanocc::TwmoPrimUnfm80,
                  // o2::aod::trackmeanocc::TwmoFV0AUnfm80,
                  // o2::aod::trackmeanocc::TwmoFV0CUnfm80,
                  // o2::aod::trackmeanocc::TwmoFT0AUnfm80,
                  // o2::aod::trackmeanocc::TwmoFT0CUnfm80,
                  // o2::aod::trackmeanocc::TwmoFDDAUnfm80,
                  // o2::aod::trackmeanocc::TwmoFDDCUnfm80,

                  // o2::aod::trackmeanocc::TwmoNTrackITSUnfm80,
                  // o2::aod::trackmeanocc::TwmoNTrackTPCUnfm80,
                  // o2::aod::trackmeanocc::TwmoNTrackTRDUnfm80,
                  // o2::aod::trackmeanocc::TwmoNTrackTOFUnfm80,
                  // o2::aod::trackmeanocc::TwmoNTrackSizeUnfm80,
                  // o2::aod::trackmeanocc::TwmoNTrackTPCAUnfm80,
                  // o2::aod::trackmeanocc::TwmoNTrackTPCCUnfm80,
                  // o2::aod::trackmeanocc::TwmoNTrackITSTPCUnfm80,
                  // o2::aod::trackmeanocc::TwmoNTrackITSTPCAUnfm80,
                  // o2::aod::trackmeanocc::TwmoNTrackITSTPCCUnfm80,

                  // o2::aod::trackmeanocc::TwmoMultNTracksHasITSUnfm80,
                  // o2::aod::trackmeanocc::TwmoMultNTracksHasTPCUnfm80,
                  // o2::aod::trackmeanocc::TwmoMultNTracksHasTOFUnfm80,
                  // o2::aod::trackmeanocc::TwmoMultNTracksHasTRDUnfm80,
                  // o2::aod::trackmeanocc::TwmoMultNTracksITSOnlyUnfm80,
                  // o2::aod::trackmeanocc::TwmoMultNTracksTPCOnlyUnfm80,
                  // o2::aod::trackmeanocc::TwmoMultNTracksITSTPCUnfm80,
                  // o2::aod::trackmeanocc::TwmoMultAllTracksTPCOnlyUnfm80,

                  // o2::aod::trackmeanocc::TwmoRobustT0V0PrimUnfm80,
                  // o2::aod::trackmeanocc::TwmoRobustFDDT0V0PrimUnfm80,
                  // o2::aod::trackmeanocc::TwmoRobustNtrackDetUnfm80,
                  // o2::aod::trackmeanocc::TwmoRobustMultExtraTableUnfm80,

                  // Lossy Compression [START]
                  o2::aod::trackmeanocc::CmoPrimUnfm80,
                  o2::aod::trackmeanocc::CmoFV0AUnfm80,
                  o2::aod::trackmeanocc::CmoFV0CUnfm80,
                  o2::aod::trackmeanocc::CmoFT0AUnfm80,
                  o2::aod::trackmeanocc::CmoFT0CUnfm80,
                  o2::aod::trackmeanocc::CmoFDDAUnfm80,
                  o2::aod::trackmeanocc::CmoFDDCUnfm80,

                  o2::aod::trackmeanocc::CmoNTrackITSUnfm80,
                  o2::aod::trackmeanocc::CmoNTrackTPCUnfm80,
                  o2::aod::trackmeanocc::CmoNTrackTRDUnfm80,
                  o2::aod::trackmeanocc::CmoNTrackTOFUnfm80,
                  o2::aod::trackmeanocc::CmoNTrackSizeUnfm80,
                  o2::aod::trackmeanocc::CmoNTrackTPCAUnfm80,
                  o2::aod::trackmeanocc::CmoNTrackTPCCUnfm80,
                  o2::aod::trackmeanocc::CmoNTrackITSTPCUnfm80,
                  o2::aod::trackmeanocc::CmoNTrackITSTPCAUnfm80,
                  o2::aod::trackmeanocc::CmoNTrackITSTPCCUnfm80,

                  o2::aod::trackmeanocc::CmoMultNTracksHasITSUnfm80,
                  o2::aod::trackmeanocc::CmoMultNTracksHasTPCUnfm80,
                  o2::aod::trackmeanocc::CmoMultNTracksHasTOFUnfm80,
                  o2::aod::trackmeanocc::CmoMultNTracksHasTRDUnfm80,
                  o2::aod::trackmeanocc::CmoMultNTracksITSOnlyUnfm80,
                  o2::aod::trackmeanocc::CmoMultNTracksTPCOnlyUnfm80,
                  o2::aod::trackmeanocc::CmoMultNTracksITSTPCUnfm80,
                  o2::aod::trackmeanocc::CmoMultAllTracksTPCOnlyUnfm80,

                  o2::aod::trackmeanocc::CmoRobustT0V0PrimUnfm80,
                  o2::aod::trackmeanocc::CmoRobustFDDT0V0PrimUnfm80,
                  o2::aod::trackmeanocc::CmoRobustNtrackDetUnfm80,
                  o2::aod::trackmeanocc::CmoRobustMultExtraTableUnfm80,

                  o2::aod::trackmeanocc::CwmoPrimUnfm80,
                  o2::aod::trackmeanocc::CwmoFV0AUnfm80,
                  o2::aod::trackmeanocc::CwmoFV0CUnfm80,
                  o2::aod::trackmeanocc::CwmoFT0AUnfm80,
                  o2::aod::trackmeanocc::CwmoFT0CUnfm80,
                  o2::aod::trackmeanocc::CwmoFDDAUnfm80,
                  o2::aod::trackmeanocc::CwmoFDDCUnfm80,

                  o2::aod::trackmeanocc::CwmoNTrackITSUnfm80,
                  o2::aod::trackmeanocc::CwmoNTrackTPCUnfm80,
                  o2::aod::trackmeanocc::CwmoNTrackTRDUnfm80,
                  o2::aod::trackmeanocc::CwmoNTrackTOFUnfm80,
                  o2::aod::trackmeanocc::CwmoNTrackSizeUnfm80,
                  o2::aod::trackmeanocc::CwmoNTrackTPCAUnfm80,
                  o2::aod::trackmeanocc::CwmoNTrackTPCCUnfm80,
                  o2::aod::trackmeanocc::CwmoNTrackITSTPCUnfm80,
                  o2::aod::trackmeanocc::CwmoNTrackITSTPCAUnfm80,
                  o2::aod::trackmeanocc::CwmoNTrackITSTPCCUnfm80,

                  o2::aod::trackmeanocc::CwmoMultNTracksHasITSUnfm80,
                  o2::aod::trackmeanocc::CwmoMultNTracksHasTPCUnfm80,
                  o2::aod::trackmeanocc::CwmoMultNTracksHasTOFUnfm80,
                  o2::aod::trackmeanocc::CwmoMultNTracksHasTRDUnfm80,
                  o2::aod::trackmeanocc::CwmoMultNTracksITSOnlyUnfm80,
                  o2::aod::trackmeanocc::CwmoMultNTracksTPCOnlyUnfm80,
                  o2::aod::trackmeanocc::CwmoMultNTracksITSTPCUnfm80,
                  o2::aod::trackmeanocc::CwmoMultAllTracksTPCOnlyUnfm80,

                  o2::aod::trackmeanocc::CwmoRobustT0V0PrimUnfm80,
                  o2::aod::trackmeanocc::CwmoRobustFDDT0V0PrimUnfm80,
                  o2::aod::trackmeanocc::CwmoRobustNtrackDetUnfm80,
                  o2::aod::trackmeanocc::CwmoRobustMultExtraTableUnfm80,
                  // Lossy Compression [END]

                  o2::aod::track::DcaXY, //	dcaXY	float	Impact parameter in XY of the track to the primary vertex
                  o2::aod::track::DcaZ,  //	dcaZ	float	Impact parameter in Z of the track to the primary vertex

                  o2::aod::track::SigmaDcaXY2, //		sigmaDcaXY2	float	Impact parameter sigma^2 in XY of the track to the primary vertex
                  o2::aod::track::SigmaDcaZ2,  //		sigmaDcaZ2	float	Impact parameter sigma^2 in Z of the track to the primary verte

                  o2::aod::pidtofevtime::UsedForTOFEvTime,
                  o2::aod::pidtofevtime::EvTimeTOF,
                  o2::aod::pidtofevtime::EvTimeTOFErr,
                  o2::aod::pidtofevtime::EvTimeTOFMult,

                  // From Track Extra
                  drTrack::TOFExpTimeEl, // D	tofExpTimeEl	float	Expected time for the track to reach the TOF under the electron hypothesis
                  drTrack::TOFExpTimeMu, // D	tofExpTimeMu	float	Expected time for the track to reach the TOF under the muon hypothesis
                  drTrack::TOFExpTimePi, // D	tofExpTimePi	float	Expected time for the track to reach the TOF under the pion hypothesis
                  drTrack::TOFExpTimeKa, // D	tofExpTimeKa	float	Expected time for the track to reach the TOF under the kaon hypothesis
                  drTrack::TOFExpTimePr, // D	tofExpTimePr	float	Expected time for the track to reach the TOF under the proton hypothesis
                  drTrack::TOFExpTimeDe, // D	tofExpTimeDe	float	Expected time for the track to reach the TOF under the deuteron hypothesis
                  drTrack::TOFExpTimeTr, // D	tofExpTimeTr	float	Expected time for the track to reach the TOF under the triton hypothesis
                  drTrack::TOFExpTimeHe, // D	tofExpTimeHe	float	Expected time for the track to reach the TOF under the helium3 hypothesis
                  drTrack::TOFExpTimeAl, // D	tofExpTimeAl

                  // drTrack::ITSNSigmaEl,
                  // drTrack::ITSNSigmaMu,
                  // drTrack::ITSNSigmaPi,
                  // drTrack::ITSNSigmaKa,
                  // drTrack::ITSNSigmaPr,
                  // drTrack::ITSNSigmaDe,
                  // drTrack::ITSNSigmaTr,
                  // drTrack::ITSNSigmaHe,
                  // drTrack::ITSNSigmaAl,

                  drTrack::CITSNSigmaEl,
                  drTrack::CITSNSigmaMu,
                  drTrack::CITSNSigmaPi,
                  drTrack::CITSNSigmaKa,
                  drTrack::CITSNSigmaPr,
                  drTrack::CITSNSigmaDe,
                  drTrack::CITSNSigmaTr,
                  drTrack::CITSNSigmaHe,
                  drTrack::CITSNSigmaAl,
                  // o2::soa::Index<>, trackqa::TrackId, trackqa::TPCTime0, trackqa::TPCDCAR, trackqa::TPCDCAZ, trackqa::TPCClusterByteMask,
                  // trackqa::TPCdEdxMax0R, trackqa::TPCdEdxMax1R, trackqa::TPCdEdxMax2R, trackqa::TPCdEdxMax3R,
                  // trackqa::TPCdEdxTot0R, trackqa::TPCdEdxTot1R, trackqa::TPCdEdxTot2R, trackqa::TPCdEdxTot3R,
                  // trackqa::DeltaRefContParamY, trackqa::DeltaRefContParamZ, trackqa::DeltaRefContParamSnp, trackqa::DeltaRefContParamTgl, trackqa::DeltaRefContParamQ2Pt,
                  // trackqa::DeltaRefGloParamY, trackqa::DeltaRefGloParamZ, trackqa::DeltaRefGloParamSnp, trackqa::DeltaRefGloParamTgl, trackqa::DeltaRefGloParamQ2Pt,
                  // trackqa::DeltaTOFdX, trackqa::DeltaTOFdZ,
                  // trackqa::IsDummy<trackqa::DeltaRefContParamY, trackqa::DeltaRefContParamZ, trackqa::DeltaRefContParamSnp, trackqa::DeltaRefContParamTgl, trackqa::DeltaRefContParamQ2Pt,
                  //                  trackqa::DeltaRefGloParamY, trackqa::DeltaRefGloParamZ, trackqa::DeltaRefGloParamSnp, trackqa::DeltaRefGloParamTgl, trackqa::DeltaRefGloParamQ2Pt>);
                  trackqa::DeltaRefContParamY,
                  trackqa::DeltaRefContParamZ,
                  trackqa::DeltaRefContParamSnp,
                  trackqa::DeltaRefContParamTgl,
                  trackqa::DeltaRefContParamQ2Pt,
                  trackqa::DeltaRefGloParamY,
                  trackqa::DeltaRefGloParamZ,
                  trackqa::DeltaRefGloParamSnp,
                  trackqa::DeltaRefGloParamTgl,
                  trackqa::DeltaRefGloParamQ2Pt,
                  trackqa::DeltaTOFdX,
                  trackqa::DeltaTOFdZ,

                  o2::aod::cent::CentFV0A,
                  o2::aod::cent::CentFT0M,
                  o2::aod::cent::CentFT0A,
                  o2::aod::cent::CentFT0C,
                  o2::aod::cent::CentFDDM,
                  o2::aod::cent::CentNTPV,
                  aod::evsel::NumTracksInTimeRange,  // trackOccupancyInTimeRange	int	Occupancy in specified time interval by a number of tracks from nearby collisions
                  aod::evsel::SumAmpFT0CInTimeRange, // ft0cOccupancyInTimeRange	float	Occupancy in specified time interval by a sum of FT0C amplitudes from nearby collisions
                  drTrack::WeightDSPt,
                  drTrack::WeightDSQPt,
                  drTrack::TriggerMaskDS,
                  drcolmn::HEntriesTrkPt,
                  drcolmn::HEntriesTrkQPt,
                  drcolmn::Pol1FitTrackPt,
                  drcolmn::Pol1FitTrackQPt

                  // ,drTrack::IsDummy
                  //  <trackqa::DeltaRefContParamY, trackqa::DeltaRefContParamZ, trackqa::DeltaRefContParamSnp, trackqa::DeltaRefContParamTgl, trackqa::DeltaRefContParamQ2Pt,
                  //                                              trackqa::DeltaRefGloParamY, trackqa::DeltaRefGloParamZ, trackqa::DeltaRefGloParamSnp, trackqa::DeltaRefGloParamTgl, trackqa::DeltaRefGloParamQ2Pt>);
);

using DrTrack = DrTracks::iterator;

// using DrPosTrack = DrTracks::iterator;
// using DrNegTrack = DrTracks::iterator;

namespace drcolmn
{
DECLARE_SOA_INDEX_COLUMN_FULL(DrPosTrack, drPosTrack, int, DrTracks, "_Pos"); //! Positive track
DECLARE_SOA_INDEX_COLUMN_FULL(DrNegTrack, drNegTrack, int, DrTracks, "_Neg"); //! Negative track

DECLARE_SOA_COLUMN(GPosTrackId, gPosTrackId, int64_t);
DECLARE_SOA_COLUMN(GNegTrackId, gNegTrackId, int64_t);
DECLARE_SOA_COLUMN(GV0Id, gV0Id, int64_t);

DECLARE_SOA_INDEX_COLUMN_FULL(DrProng0, drProng0, int, DrTracks, ""); //! Positive track
DECLARE_SOA_INDEX_COLUMN_FULL(DrProng1, drProng1, int, DrTracks, ""); //! Negative track

DECLARE_SOA_COLUMN(GProng0Id, gProng0Id, int64_t);
DECLARE_SOA_COLUMN(GProng1Id, gProng1Id, int64_t);
DECLARE_SOA_COLUMN(GD0Id, gD0Id, int64_t);
// DECLARE_SOA_COLUMN(DrProng0Id, drProng0Id, int64_t);
// DECLARE_SOA_COLUMN(DrProng1Id, drProng1Id, int64_t);

DECLARE_SOA_INDEX_COLUMN_FULL(DrPosKaon, drPosKaon, int, DrTracks, "_PosKaon"); //! Positive track
DECLARE_SOA_INDEX_COLUMN_FULL(DrNegKaon, drNegKaon, int, DrTracks, "_NegKaon"); //! Negative track

DECLARE_SOA_COLUMN(PosKaonId, posKaon, int64_t); //! Positive track
DECLARE_SOA_COLUMN(NegKaonId, negKaon, int64_t); //! Negative track

DECLARE_SOA_COLUMN(GPosKaonId, gPosKaonId, int64_t);
DECLARE_SOA_COLUMN(GNegKaonId, gNegKaonId, int64_t);
DECLARE_SOA_COLUMN(GPhiId, gPhiId, int64_t);

DECLARE_SOA_COLUMN(CollisionId, collisionId, int64_t);

} // namespace drcolmn

// o2::aod::V0Datas = soa::Join<o2::aod::V0Indices, o2::aod::V0TrackXs, o2::aod::V0Cores>

DECLARE_SOA_TABLE(DrV0s, "AOD", "DRV0s", o2::soa::Index<>,
                  o2::aod::drcolmn::GlobalIndex,
                  o2::aod::drcolmn::OriginalIndex,
                  o2::aod::drcolmn::GBCid,
                  o2::aod::drcolmn::GDFId,
                  o2::aod::drcolmn::GCollId,
                  o2::aod::drcolmn::GV0Id,
                  o2::aod::drcolmn::GPosTrackId,
                  o2::aod::drcolmn::GNegTrackId,

                  o2::aod::drcolmn::DrPosTrackId,
                  o2::aod::drcolmn::DrNegTrackId,
                  o2::aod::v0data::PosTrackId,  // I 	posTrackId 	int 	Pointer into Tracks
                  o2::aod::v0data::NegTrackId,  // I 	negTrackId 	int 	Pointer into Tracks
                  o2::aod::v0data::CollisionId, // I 	collisionId 	int32 	Pointer into Collisions
                  o2::aod::v0data::V0Id,        // I 	v0Id 	int32 	Pointer into V0s

                  o2::aod::v0data::PosX, // posX 	float 	positive track X at min
                  o2::aod::v0data::NegX, // negX 	float 	negative track X at min

                  // //third Table
                  o2::aod::v0data::X, // x 	float 	decay position X
                  o2::aod::v0data::Y, // y 	float 	decay position Y
                  o2::aod::v0data::Z, // z 	float 	decay position Z
                  // o2::aod::v0data::PxPos, 		pxpos 	float 	positive track px at min
                  // o2::aod::v0data::PyPos, 		pypos 	float 	positive track py at min
                  // o2::aod::v0data::PzPos, 		pzpos 	float 	positive track pz at min
                  // o2::aod::v0data::PxNeg, 		pxneg 	float 	negative track px at min
                  // o2::aod::v0data::PyNeg, 		pyneg 	float 	negative track py at min
                  // o2::aod::v0data::PzNeg, 		pzneg 	float 	negative track pz at min
                  o2::aod::v0data::DCAV0Daughters, // dcaV0daughters 	float 	DCA between V0 daughters
                  o2::aod::v0data::DCAPosToPV,     // dcapostopv 	float 	DCA positive prong to PV
                  o2::aod::v0data::DCANegToPV,     // dcanegtopv 	float 	DCA negative prong to PV
                  o2::aod::v0data::V0CosPA,        // v0cosPA 	float 	V0 CosPA
                  o2::aod::v0data::DCAV0ToPV,      // dcav0topv 	float 	DCA V0 to PV

                  o2::aod::drcolmn::V0Radius,         // D 	v0radius 	float 	V0 decay radius (2D, centered at zero)
                  o2::aod::drcolmn::MLambda,          // D 	mLambda 	float 	mass under lambda hypothesis
                  o2::aod::drcolmn::MAntiLambda,      // D 	mAntiLambda 	float 	mass under antilambda hypothesis
                  o2::aod::drcolmn::MK0Short,         // D 	mK0Short 	float 	mass under K0short hypothesis
                  o2::aod::drcolmn::MGamma,           // D 	mGamma 	float 	mass under gamma hypothesis
                  o2::aod::drcolmn::MHypertriton,     // D 	mHypertriton 	float 	mass under hypertriton hypothesis
                  o2::aod::drcolmn::MAntiHypertriton, // D 	mAntiHypertriton 	float 	mass under antihypertriton hypothesis
                  o2::aod::drcolmn::NegativePt,       // D 	negativept 	float 	negative daughter pT
                  o2::aod::drcolmn::PositivePt,       // D 	positivept 	float 	positive daughter pT
                  o2::aod::drcolmn::NegativeEta,      // D 	negativeeta 	float 	negative daughter eta
                  o2::aod::drcolmn::NegativePhi,      // D 	negativephi 	float 	negative daughter phi
                  o2::aod::drcolmn::PositiveEta,      // D 	positiveeta 	float 	positive daughter eta
                  o2::aod::drcolmn::PositivePhi,      // D 	positivephi 	float 	positive daughter phi
                  o2::aod::drcolmn::IsStandardV0,     // D 	isStandardV0 	bool 	is standard V0
                  o2::aod::drcolmn::IsPhotonTPConly,  // D 	isPhotonTPConly 	bool 	is tpc-only photon V0
                  o2::aod::timestamp::Timestamp,
                  o2::aod::drcolmn::TfId,
                  o2::aod::drcolmn::BCinTF,

                  o2::aod::drcolmn::Px,
                  o2::aod::drcolmn::Py,
                  o2::aod::drcolmn::Pz,

                  o2::aod::drcolmn::Pt,
                  o2::aod::drcolmn::Eta,
                  o2::aod::drcolmn::Phi,

                  o2::aod::drcolmn::Alpha,
                  o2::aod::drcolmn::QtArm,
                  // o2::aod::drcolmn::DistOverTotMom,
                  o2::aod::drcolmn::PsiPair,
                  o2::aod::drcolmn::PFracPos,
                  o2::aod::drcolmn::PFracNeg,
                  //
                  o2::aod::drcolmn::WeightDSPtV0,
                  o2::aod::drcolmn::WeightDSPtLambda,
                  o2::aod::drcolmn::WeightDSPtAntiLambda,
                  o2::aod::drcolmn::WeightDSPtK0Short,
                  o2::aod::drcolmn::WeightDSPtGamma,
                  o2::aod::drcolmn::WeightDSPtHypertriton,
                  o2::aod::drcolmn::WeightDSPtAntiHypertriton,

                  o2::aod::drcolmn::TriggerMaskV0,
                  // o2::aod::drcolmn::TriggerMaskV0PT,

                  o2::aod::drcolmn::HEntriesV0Pt,
                  o2::aod::drcolmn::HEntriesLambdaPt,
                  o2::aod::drcolmn::HEntriesAntiLambdaPt,
                  o2::aod::drcolmn::HEntriesK0ShortPt,

                  o2::aod::drcolmn::Pol1FitV0Pt,
                  o2::aod::drcolmn::Pol1FitLambdaPt,
                  o2::aod::drcolmn::Pol1FitAntiLambdaPt,
                  o2::aod::drcolmn::Pol1FitK0ShortPt

                  // o2::aod::drcolmn::HEntriesGammaPt
                  // o2::aod::drcolmn::HEntriesHypertritonPt
                  // o2::aod::drcolmn::HEntriesAntiHypertritonPt

                  // o2::aod::v0data::V0Radius
                  // // o2::aod::v0data::V0Type 		v0Type 	uint8_t 	type of V0. 0: built solely for cascades (does not pass standard V0 cuts), 1: standard 2, 3: photon-like with TPC-only use. Regular analysis should always use type 1.
                  // // o2::aod::v0data::PtHypertriton 	D 	ptHypertriton 	float 	V0 pT
                  // // o2::aod::v0data::PtAntiHypertriton 	D 	ptAntiHypertriton 	float 	V0 pT
                  //----> o2::aod::v0data::V0Radius 	//D 	v0radius 	float 	V0 decay radius (2D, centered at zero)
                  // // o2::aod::v0data::DistOverTotMom, 	D 	distovertotmom 	? 	PV to V0decay distance over total momentum
                  // // o2::aod::v0data::Alpha, 	D 	alpha 	? 	Armenteros Alpha
                  // // o2::aod::v0data::QtArm, 	D 	qtarm 	? 	Armenteros Qt
                  // // o2::aod::v0data::PsiPair, 	D 	psipair 	? 	psi pair angle
                  // // o2::aod::v0data::PFracPos, 	D 	pfracpos 	?
                  // // o2::aod::v0data::PFracNeg, 	D 	pfracneg 	?
                  // o2::aod::v0data::MLambda, 	          //D 	mLambda 	float 	mass under lambda hypothesis
                  // o2::aod::v0data::MAntiLambda, 	      //D 	mAntiLambda 	float 	mass under antilambda hypothesis
                  // o2::aod::v0data::MK0Short, 	        //D 	mK0Short 	float 	mass under K0short hypothesis
                  // o2::aod::v0data::MGamma, 	          //D 	mGamma 	float 	mass under gamma hypothesis
                  // o2::aod::v0data::MHypertriton, 	    //D 	mHypertriton 	float 	mass under hypertriton hypothesis
                  // o2::aod::v0data::MAntiHypertriton, 	//D 	mAntiHypertriton 	float 	mass under antihypertriton hypothesis
                  // // o2::aod::v0data::M 	D 	m 	float 	mass under a certain hypothesis (0:K0, 1:L, 2:Lbar, 3:gamma, 4:hyp, 5:ahyp)
                  // // o2::aod::v0data::YK0Short, 	D 	yK0Short 	float 	V0 y with K0short hypothesis
                  // // o2::aod::v0data::YLambda, 	D 	yLambda 	float 	V0 y with lambda or antilambda hypothesis
                  // // o2::aod::v0data::YHypertriton, 	D 	yHypertriton 	float 	V0 y with hypertriton hypothesis
                  // // o2::aod::v0data::YAntiHypertriton, 	D 	yAntiHypertriton 	float 	V0 y with antihypertriton hypothesis
                  // // o2::aod::v0data::Rapidity, 	D 	rapidity 	float 	rapidity (0:K0, 1:L, 2:Lbar)
                  // o2::aod::v0data::NegativePt, 	    //D 	negativept 	float 	negative daughter pT
                  // o2::aod::v0data::PositivePt, 	    //D 	positivept 	float 	positive daughter pT
                  // o2::aod::v0data::NegativeEta, 	    //D 	negativeeta 	float 	negative daughter eta
                  // o2::aod::v0data::NegativePhi, 	    //D 	negativephi 	float 	negative daughter phi
                  // o2::aod::v0data::PositiveEta, 	    //D 	positiveeta 	float 	positive daughter eta
                  // o2::aod::v0data::PositivePhi, 	    //D 	positivephi 	float 	positive daughter phi
                  // o2::aod::v0data::IsStandardV0, 	  //D 	isStandardV0 	bool 	is standard V0
                  // o2::aod::v0data::IsPhotonTPConly, 	//D 	isPhotonTPConly 	bool 	is tpc-only photon V0
);

DECLARE_SOA_TABLE(DrPhis, "AOD", "DRPhis", o2::soa::Index<>,
                  o2::aod::drcolmn::GlobalIndex,
                  o2::aod::drcolmn::OriginalIndex,
                  o2::aod::drcolmn::GBCid,
                  o2::aod::drcolmn::GDFId,
                  o2::aod::drcolmn::GCollId,
                  o2::aod::drcolmn::GPhiId,

                  o2::aod::drcolmn::GPosKaonId, // global Id in merged DF in derived Tracks
                  o2::aod::drcolmn::GNegKaonId, // global Id in merged DF in derived Tracks

                  o2::aod::drcolmn::DrPosKaonId, // global Id in single DF in derived Tracks
                  o2::aod::drcolmn::DrNegKaonId, // global Id in single DF in derived Tracks

                  o2::aod::drcolmn::PosKaonId, // original Id in single DF in original Tracks
                  o2::aod::drcolmn::NegKaonId, // original Id in single DF in original Tracks

                  o2::aod::drcolmn::CollisionId, // I 	collisionId 	int32 	Pointer into Collisions

                  o2::aod::drcolmn::Px,
                  o2::aod::drcolmn::Py,
                  o2::aod::drcolmn::Pz,
                  o2::aod::drcolmn::Pt,
                  o2::aod::drcolmn::Eta,
                  o2::aod::drcolmn::Phi,
                  o2::aod::drcolmn::Y,
                  o2::aod::drcolmn::MPhi,
                  o2::aod::drcolmn::E,

                  o2::aod::drcolmn::PosKaonPx,
                  o2::aod::drcolmn::PosKaonPy,
                  o2::aod::drcolmn::PosKaonPz,
                  o2::aod::drcolmn::PosKaonPt,
                  o2::aod::drcolmn::PosKaonEta,
                  // o2::aod::drcolmn::PosKaonY,
                  // o2::aod::drcolmn::PosKaonIsPVContributor,
                  // o2::aod::drcolmn::PosKaonIsGlobalTrack,

                  o2::aod::drcolmn::NegKaonPx,
                  o2::aod::drcolmn::NegKaonPy,
                  o2::aod::drcolmn::NegKaonPz,
                  o2::aod::drcolmn::NegKaonPt,
                  o2::aod::drcolmn::NegKaonEta,
                  // o2::aod::drcolmn::NegKaonY,
                  // o2::aod::drcolmn::NegKaonIsPVContributor,
                  // o2::aod::drcolmn::NegKaonIsGlobalTrack,

                  o2::aod::timestamp::Timestamp,
                  o2::aod::drcolmn::TfId,
                  o2::aod::drcolmn::BCinTF,

                  o2::aod::drcolmn::WeightDSPtPhi,
                  o2::aod::drcolmn::TriggerMaskPhi,
                  o2::aod::drcolmn::HEntriesPhiPt,

                  drcolmn::Pol1FitPhiPt

);

DECLARE_SOA_TABLE(DrD0s, "AOD", "DRD0s", o2::soa::Index<>,
                  o2::aod::drcolmn::GlobalIndex,
                  o2::aod::drcolmn::OriginalIndex,
                  o2::aod::drcolmn::GBCid,
                  o2::aod::drcolmn::GDFId,
                  o2::aod::drcolmn::GCollId,
                  o2::aod::drcolmn::GD0Id,
                  o2::aod::drcolmn::GProng0Id,
                  o2::aod::drcolmn::GProng1Id,
                  o2::aod::drcolmn::DrProng0Id,
                  o2::aod::drcolmn::DrProng1Id,

                  o2::aod::hf_track_index::Prong0Id,            // I 	prong0Id 	int 	Index to first prong
                  o2::aod::hf_track_index::Prong1Id,            // I 	prong1Id 	int 	Index to second prong
                  o2::aod::hf_cand::CollisionId,                // I 	collisionId 	int32 	Pointer into Collisions
                  o2::aod::collision::PosX,                     // posX 	float 	X Vertex position in cm
                  o2::aod::collision::PosY,                     // posY 	float 	Y Vertex position in cm
                  o2::aod::collision::PosZ,                     // posZ 	float 	Z Vertex position in cm
                  o2::aod::hf_cand::XSecondaryVertex,           // xSecondaryVertex 	float
                  o2::aod::hf_cand::YSecondaryVertex,           // ySecondaryVertex 	float
                  o2::aod::hf_cand::ZSecondaryVertex,           // zSecondaryVertex 	float
                  o2::aod::hf_cand::ErrorDecayLength,           // errorDecayLength 	float
                  o2::aod::hf_cand::ErrorDecayLengthXY,         // errorDecayLengthXY 	float
                  o2::aod::hf_cand::Chi2PCA,                    // chi2PCA 	float 	sum of (non-weighted) distances of the secondary vertex to its prongs
                  o2::aod::drcolmn::RSecondaryVertex,           // GI 		?
                  o2::aod::drcolmn::DecayLength,                // D 	decayLength 	float
                  o2::aod::drcolmn::DecayLengthXY,              // D 	decayLengthXY 	float
                  o2::aod::drcolmn::DecayLengthNormalised,      // D 	decayLengthNormalised 	float
                  o2::aod::drcolmn::DecayLengthXYNormalised,    // D 	decayLengthXYNormalised 	float
                  o2::aod::drcolmn::ImpactParameterNormalised0, // GI 		?
                  o2::aod::drcolmn::PtProng0,                   // D 	ptProng0 	float
                  o2::aod::drcolmn::Pt2Prong0,                  // D 	pt2Prong0 	float
                  // o2::aod::drcolmn::PVectorProng0, 	              //D 	pVectorProng0 	std::array<float,3>
                  o2::aod::drcolmn::ImpactParameterNormalised1, // GI 		?
                  o2::aod::drcolmn::PtProng1,                   // D 	ptProng1 	float
                  o2::aod::drcolmn::Pt2Prong1,                  // D 	pt2Prong1 	float
                  // o2::aod::hf_cand::PVectorProng1, 	              //D 	pVectorProng1 	std::array<float,3>
                  o2::aod::hf_cand::PxProng0,               // pxProng0 	float
                  o2::aod::hf_cand::PyProng0,               // pyProng0 	float
                  o2::aod::hf_cand::PzProng0,               // pzProng0 	float
                  o2::aod::hf_cand::PxProng1,               // pxProng1 	float
                  o2::aod::hf_cand::PyProng1,               // pyProng1 	float
                  o2::aod::hf_cand::PzProng1,               // pzProng1 	float
                  o2::aod::hf_cand::ImpactParameter0,       // impactParameter0 	float
                  o2::aod::hf_cand::ImpactParameter1,       // impactParameter1 	float
                  o2::aod::hf_cand::ErrorImpactParameter0,  // errorImpactParameter0 	float
                  o2::aod::hf_cand::ErrorImpactParameter1,  // errorImpactParameter1 	float
                  o2::aod::hf_cand::ImpactParameterZ0,      // impactParameterZ0 	float
                  o2::aod::hf_cand::ImpactParameterZ1,      // impactParameterZ1 	float
                  o2::aod::hf_cand::ErrorImpactParameterZ0, // errorImpactParameterZ0 	float
                  o2::aod::hf_cand::ErrorImpactParameterZ1, // errorImpactParameterZ1 	float
                  o2::aod::hf_cand::NProngsContributorsPV,  // nProngsContributorsPV 	uint8_t 	number of prongs contributing to the primary-vertex reconstruction
                  o2::aod::hf_track_index::HFflag,          // hfflag 	uint8_t
                  // o2::aod::hf_cand_2prong::M, 	                          //    GI 		?
                  // o2::aod::hf_cand_2prong::M2, 	                          //D 	m2 	float
                  o2::aod::drcolmn::ImpactParameterProduct, // D 	impactParameterProduct 	float
                  // // o2::aod::hf_cand_2prong::CosThetaStar, 	                //D 	cosThetaStar 	float
                  o2::aod::drcolmn::ImpactParameterProngSqSum, // D 	impactParameterProngSqSum 	float
                  o2::aod::drcolmn::Pt,                        // GI 		?
                  o2::aod::drcolmn::Pt2,                       // D 	pt2 	float
                  o2::aod::drcolmn::P,                         // D 	p 	float
                  o2::aod::drcolmn::P2,                        // D 	p2 	float
                  // // // o2::aod::hf_cand::PVector, 	                            //D 	pVector 	std::array<float,3>
                  o2::aod::drcolmn::CPA,   // D 	cpa 	float
                  o2::aod::drcolmn::CPAXY, // D 	cpaXY 	float
                  // // o2::aod::hf_cand::Ct, 	                                //D 	ct 	float
                  o2::aod::drcolmn::ImpactParameterXY,    // D 	impactParameterXY 	float
                  o2::aod::drcolmn::MaxNormalisedDeltaIP, // D 	maxNormalisedDeltaIP 	float
                  o2::aod::drcolmn::Px,
                  o2::aod::drcolmn::Py,
                  o2::aod::drcolmn::Pz,

                  o2::aod::drcolmn::Eta, // D 	eta 	float
                  o2::aod::drcolmn::Phi, // D 	phi 	float
                  o2::aod::drcolmn::Y,   // D 	y 	float
                  o2::aod::drcolmn::E,   // D 	e 	float

                  o2::aod::timestamp::Timestamp,
                  o2::aod::drcolmn::TfId,
                  o2::aod::drcolmn::BCinTF,

                  o2::aod::drcolmn::WeightDSPtD0,
                  o2::aod::drcolmn::WeightDSPtD0Bar,
                  o2::aod::drcolmn::TriggerMaskD0,

                  o2::aod::drcolmn::HEntriesD0Pt,
                  o2::aod::drcolmn::HEntriesD0BarPt,

                  o2::aod::drcolmn::Pol1FitD0Pt,
                  o2::aod::drcolmn::Pol1FitD0BarPt

);

DECLARE_SOA_TABLE(OCCS, "AOD", "OCCS", o2::soa::Index<>,
                  o2::aod::occp::TfId,
                  o2::aod::drcolmn::GDFId,
                  o2::aod::occp::BcsInTFList,

                  o2::aod::occp::OccPrimUnfm80,
                  o2::aod::occp::OccFV0AUnfm80,
                  o2::aod::occp::OccFV0CUnfm80,
                  o2::aod::occp::OccFT0AUnfm80,
                  o2::aod::occp::OccFT0CUnfm80,
                  o2::aod::occp::OccFDDAUnfm80,
                  o2::aod::occp::OccFDDCUnfm80,

                  o2::aod::occp::OccNTrackITSUnfm80,
                  o2::aod::occp::OccNTrackTPCUnfm80,
                  o2::aod::occp::OccNTrackTRDUnfm80,
                  o2::aod::occp::OccNTrackTOFUnfm80,
                  o2::aod::occp::OccNTrackSizeUnfm80,
                  o2::aod::occp::OccNTrackTPCAUnfm80,
                  o2::aod::occp::OccNTrackTPCCUnfm80,
                  o2::aod::occp::OccNTrackITSTPCUnfm80,
                  o2::aod::occp::OccNTrackITSTPCAUnfm80,
                  o2::aod::occp::OccNTrackITSTPCCUnfm80,

                  o2::aod::occp::OccMultNTracksHasITSUnfm80,
                  o2::aod::occp::OccMultNTracksHasTPCUnfm80,
                  o2::aod::occp::OccMultNTracksHasTOFUnfm80,
                  o2::aod::occp::OccMultNTracksHasTRDUnfm80,
                  o2::aod::occp::OccMultNTracksITSOnlyUnfm80,
                  o2::aod::occp::OccMultNTracksTPCOnlyUnfm80,
                  o2::aod::occp::OccMultNTracksITSTPCUnfm80,
                  o2::aod::occp::OccMultAllTracksTPCOnlyUnfm80,

                  o2::aod::occp::OccRobustT0V0PrimUnfm80,
                  o2::aod::occp::OccRobustFDDT0V0PrimUnfm80,
                  o2::aod::occp::OccRobustNtrackDetUnfm80,
                  o2::aod::occp::OccRobustMultExtraTableUnfm80,

                  o2::aod::occp::MeanOccPrimUnfm80,
                  o2::aod::occp::MeanOccFV0AUnfm80,
                  o2::aod::occp::MeanOccFV0CUnfm80,
                  o2::aod::occp::MeanOccFT0AUnfm80,
                  o2::aod::occp::MeanOccFT0CUnfm80,
                  o2::aod::occp::MeanOccFDDAUnfm80,
                  o2::aod::occp::MeanOccFDDCUnfm80,

                  o2::aod::occp::MeanOccNTrackITSUnfm80,
                  o2::aod::occp::MeanOccNTrackTPCUnfm80,
                  o2::aod::occp::MeanOccNTrackTRDUnfm80,
                  o2::aod::occp::MeanOccNTrackTOFUnfm80,
                  o2::aod::occp::MeanOccNTrackSizeUnfm80,
                  o2::aod::occp::MeanOccNTrackTPCAUnfm80,
                  o2::aod::occp::MeanOccNTrackTPCCUnfm80,
                  o2::aod::occp::MeanOccNTrackITSTPCUnfm80,
                  o2::aod::occp::MeanOccNTrackITSTPCAUnfm80,
                  o2::aod::occp::MeanOccNTrackITSTPCCUnfm80,

                  o2::aod::occp::MeanOccMultNTracksHasITSUnfm80,
                  o2::aod::occp::MeanOccMultNTracksHasTPCUnfm80,
                  o2::aod::occp::MeanOccMultNTracksHasTOFUnfm80,
                  o2::aod::occp::MeanOccMultNTracksHasTRDUnfm80,
                  o2::aod::occp::MeanOccMultNTracksITSOnlyUnfm80,
                  o2::aod::occp::MeanOccMultNTracksTPCOnlyUnfm80,
                  o2::aod::occp::MeanOccMultNTracksITSTPCUnfm80,
                  o2::aod::occp::MeanOccMultAllTracksTPCOnlyUnfm80,

                  o2::aod::occp::MeanOccRobustT0V0PrimUnfm80,
                  o2::aod::occp::MeanOccRobustFDDT0V0PrimUnfm80,
                  o2::aod::occp::MeanOccRobustNtrackDetUnfm80,
                  o2::aod::occp::MeanOccRobustMultExtraTableUnfm80);

// Z0 reconstruction table
DECLARE_SOA_TABLE(DrZ0s, "AOD", "DRZ0S", o2::soa::Index<>,
                  o2::aod::drcolmn::GlobalIndex,
                  o2::aod::drcolmn::OriginalIndex,
                  o2::aod::drcolmn::GBCid,
                  o2::aod::drcolmn::GDFId,
                  o2::aod::drcolmn::GCollId,

                  o2::aod::drcolmn::GZ0Id,
                  o2::aod::drcolmn::GPosDauId,
                  o2::aod::drcolmn::GNegDauId,

                  o2::aod::drcolmn::DrPosDauId,
                  o2::aod::drcolmn::DrNegDauId,

                  o2::aod::drcolmn::PosTrackId, // original Id in single DF in original Tracks
                  o2::aod::drcolmn::NegTrackId, // original Id in single DF in original Tracks

                  o2::aod::drcolmn::DecayFlag, // muon decay or electron decay

                  o2::aod::drcolmn::Px,
                  o2::aod::drcolmn::Py,
                  o2::aod::drcolmn::Pz,
                  o2::aod::drcolmn::Pt,
                  o2::aod::drcolmn::Eta,
                  o2::aod::drcolmn::Phi,
                  o2::aod::drcolmn::Y,
                  o2::aod::drcolmn::MZ0,
                  o2::aod::drcolmn::E,

                  o2::aod::drcolmn::PosDauPx,
                  o2::aod::drcolmn::PosDauPy,
                  o2::aod::drcolmn::PosDauPz,
                  o2::aod::drcolmn::PosDauPt,
                  o2::aod::drcolmn::PosDauEta,
                  // o2::aod::drcolmn::PosKaonY,
                  // o2::aod::drcolmn::PosKaonIsPVContributor,
                  // o2::aod::drcolmn::PosKaonIsGlobalTrack,

                  o2::aod::drcolmn::NegDauPx,
                  o2::aod::drcolmn::NegDauPy,
                  o2::aod::drcolmn::NegDauPz,
                  o2::aod::drcolmn::NegDauPt,
                  o2::aod::drcolmn::NegDauEta,
                  // o2::aod::drcolmn::NegKaonY,
                  // o2::aod::drcolmn::NegKaonIsPVContributor,
                  // o2::aod::drcolmn::NegKaonIsGlobalTrack,

                  o2::aod::timestamp::Timestamp,
                  o2::aod::drcolmn::TfId,
                  o2::aod::drcolmn::BCinTF);

namespace drcolmn
{
DECLARE_SOA_COLUMN(UTrkX, uTrkX, float);
DECLARE_SOA_COLUMN(UTrkY, uTrkY, float);
DECLARE_SOA_COLUMN(UTrkZ, uTrkZ, float);
DECLARE_SOA_COLUMN(UTrkSnp, uTrkSnp, float);
DECLARE_SOA_COLUMN(UTrkSigned1Pt, uTrkSigned1Pt, float);
DECLARE_SOA_COLUMN(UTrkP, uTrkP, float);
DECLARE_SOA_COLUMN(UTrkIsWithinBeamPipe, uTrkIsWithinBeamPipe, uint8_t);

DECLARE_SOA_COLUMN(UTrkSigmaY, uTrkSigmaY, float);
DECLARE_SOA_COLUMN(UTrkSigmaZ, uTrkSigmaZ, float);
DECLARE_SOA_COLUMN(UTrkSigmaSnp, uTrkSigmaSnp, float);
DECLARE_SOA_COLUMN(UTrkSigmaTgl, uTrkSigmaTgl, float);
DECLARE_SOA_COLUMN(UTrkSigma1Pt, uTrkSigma1Pt, float);

DECLARE_SOA_COLUMN(UTrkRhoZY, uTrkRhoZY, float);
DECLARE_SOA_COLUMN(UTrkRhoSnpY, uTrkRhoSnpY, float);
DECLARE_SOA_COLUMN(UTrkRhoSnpZ, uTrkRhoSnpZ, float);
DECLARE_SOA_COLUMN(UTrkRhoTglY, uTrkRhoTglY, float);
DECLARE_SOA_COLUMN(UTrkRhoTglZ, uTrkRhoTglZ, float);
DECLARE_SOA_COLUMN(UTrkRhoTglSnp, uTrkRhoTglSnp, float);
DECLARE_SOA_COLUMN(UTrkRho1PtY, uTrkRho1PtY, float);
DECLARE_SOA_COLUMN(UTrkRho1PtZ, uTrkRho1PtZ, float);
DECLARE_SOA_COLUMN(UTrkRho1PtSnp, uTrkRho1PtSnp, float);
DECLARE_SOA_COLUMN(UTrkRho1PtTgl, uTrkRho1PtTgl, float);

DECLARE_SOA_COLUMN(UTrkCYY, uTrkCYY, float);
DECLARE_SOA_COLUMN(UTrkCZY, uTrkCZY, float);
DECLARE_SOA_COLUMN(UTrkCZZ, uTrkCZZ, float);
DECLARE_SOA_COLUMN(UTrkCSnpY, uTrkCSnpY, float);
DECLARE_SOA_COLUMN(UTrkCSnpZ, uTrkCSnpZ, float);
DECLARE_SOA_COLUMN(UTrkCSnpSnp, uTrkCSnpSnp, float);
DECLARE_SOA_COLUMN(UTrkCTglY, uTrkCTglY, float);
DECLARE_SOA_COLUMN(UTrkCTglZ, uTrkCTglZ, float);
DECLARE_SOA_COLUMN(UTrkCTglSnp, uTrkCTglSnp, float);
DECLARE_SOA_COLUMN(UTrkCTglTgl, uTrkCTglTgl, float);
DECLARE_SOA_COLUMN(UTrkC1PtY, uTrkC1PtY, float);
DECLARE_SOA_COLUMN(UTrkC1PtZ, uTrkC1PtZ, float);
DECLARE_SOA_COLUMN(UTrkC1PtSnp, uTrkC1PtSnp, float);
DECLARE_SOA_COLUMN(UTrkC1PtTgl, uTrkC1PtTgl, float);
DECLARE_SOA_COLUMN(UTrkC1Pt21Pt2, uTrkC1Pt21Pt2, float);

DECLARE_SOA_COLUMN(UTrkTpcInnerParam, uTrkTpcInnerParam, float); // if structured, adjust type
DECLARE_SOA_COLUMN(UTrkFlags, uTrkFlags, uint32_t);
DECLARE_SOA_COLUMN(UTrkItsClusterSizes, uTrkItsClusterSizes, uint32_t);
DECLARE_SOA_COLUMN(UTrkTpcNClsFindable, uTrkTpcNClsFindable, uint8_t);
DECLARE_SOA_COLUMN(UTrkTpcNClsFindableMinusFound, uTrkTpcNClsFindableMinusFound, int8_t);
DECLARE_SOA_COLUMN(UTrkTpcNClsFindableMinusCrossedRows, uTrkTpcNClsFindableMinusCrossedRows, int8_t);
DECLARE_SOA_COLUMN(UTrkTpcNClsShared, uTrkTpcNClsShared, uint8_t);
DECLARE_SOA_COLUMN(UTrkTrdPattern, uTrkTrdPattern, uint8_t);

DECLARE_SOA_COLUMN(UTrkDetectorMask, uTrkDetectorMask, uint8_t);

DECLARE_SOA_COLUMN(UTrkItsChi2NCl, uTrkItsChi2NCl, float);
DECLARE_SOA_COLUMN(UTrkTpcChi2NCl, uTrkTpcChi2NCl, float);
DECLARE_SOA_COLUMN(UTrkTrdChi2, uTrkTrdChi2, float);
DECLARE_SOA_COLUMN(UTrkTofChi2, uTrkTofChi2, float);
DECLARE_SOA_COLUMN(UTrkTpcSignal, uTrkTpcSignal, float);
DECLARE_SOA_COLUMN(UTrkTrdSignal, uTrkTrdSignal, float);
DECLARE_SOA_COLUMN(UTrkLength, uTrkLength, float);
DECLARE_SOA_COLUMN(UTrkTofExpMom, uTrkTofExpMom, float);
DECLARE_SOA_COLUMN(UTrkEtaEmcal, uTrkEtaEmcal, float);
DECLARE_SOA_COLUMN(UTrkPhiEmcal, uTrkPhiEmcal, float);
DECLARE_SOA_COLUMN(UTrkTime, uTrkTime, float);
DECLARE_SOA_COLUMN(UTrkTimeRes, uTrkTimeRes, float);
DECLARE_SOA_COLUMN(UTrkPidForTracking, uTrkPidForTracking, int8_t);
DECLARE_SOA_COLUMN(UTrkIsPVContributor, uTrkIsPVContributor, uint8_t);
DECLARE_SOA_COLUMN(UTrkTpcNClsFound, uTrkTPCNClsFound, uint8_t);
DECLARE_SOA_COLUMN(UTrkTpcNClsCrossedRows, uTrkTPCNClsCrossedRows, uint8_t);
DECLARE_SOA_COLUMN(UTrkItsClusterMap, uTrkITSClusterMap, uint8_t);
DECLARE_SOA_COLUMN(UTrkItsNCls, uTrkITSNCls, uint8_t);
DECLARE_SOA_COLUMN(UTrkItsNClsInnerBarrel, uTrkITSNClsInnerBarrel, uint8_t);
DECLARE_SOA_COLUMN(UTrkTpcCrossedRowsOverFindableCls, uTrkTPCCrossedRowsOverFindableCls, float);
DECLARE_SOA_COLUMN(UTrkTpcFoundOverFindableCls, uTrkTPCFoundOverFindableCls, float);
DECLARE_SOA_COLUMN(UTrkTpcFractionSharedCls, uTrkTPCFractionSharedCls, float);
DECLARE_SOA_COLUMN(UTrkDetectorMap, uTrkDetectorMap, uint32_t);

DECLARE_SOA_COLUMN(UTrkTpcTime0, uTrkTpcTime0, float);
DECLARE_SOA_COLUMN(UTrkTpcDcaR, uTrkTpcDcaR, float);
DECLARE_SOA_COLUMN(UTrkTpcDcaZ, uTrkTpcDcaZ, float);
DECLARE_SOA_COLUMN(UTrkTpcClusterByteMask, uTrkTpcClusterByteMask, uint8_t);

DECLARE_SOA_COLUMN(UTrkTpcdEdxMax0R, uTrkTpcdEdxMax0R, float);
DECLARE_SOA_COLUMN(UTrkTpcdEdxMax1R, uTrkTpcdEdxMax1R, float);
DECLARE_SOA_COLUMN(UTrkTpcdEdxMax2R, uTrkTpcdEdxMax2R, float);
DECLARE_SOA_COLUMN(UTrkTpcdEdxMax3R, uTrkTpcdEdxMax3R, float);
DECLARE_SOA_COLUMN(UTrkTpcdEdxTot0R, uTrkTpcdEdxTot0R, float);
DECLARE_SOA_COLUMN(UTrkTpcdEdxTot1R, uTrkTpcdEdxTot1R, float);
DECLARE_SOA_COLUMN(UTrkTpcdEdxTot2R, uTrkTpcdEdxTot2R, float);
DECLARE_SOA_COLUMN(UTrkTpcdEdxTot3R, uTrkTpcdEdxTot3R, float);

DECLARE_SOA_COLUMN(UTrkDeltaRefContParamY, uTrkDeltaRefContParamY, float);
DECLARE_SOA_COLUMN(UTrkDeltaRefContParamZ, uTrkDeltaRefContParamZ, float);
DECLARE_SOA_COLUMN(UTrkDeltaRefContParamSnp, uTrkDeltaRefContParamSnp, float);
DECLARE_SOA_COLUMN(UTrkDeltaRefContParamTgl, uTrkDeltaRefContParamTgl, float);
DECLARE_SOA_COLUMN(UTrkDeltaRefContParamQ2Pt, uTrkDeltaRefContParamQ2Pt, float);
DECLARE_SOA_COLUMN(UTrkDeltaRefGloParamY, uTrkDeltaRefGloParamY, float);
DECLARE_SOA_COLUMN(UTrkDeltaRefGloParamZ, uTrkDeltaRefGloParamZ, float);
DECLARE_SOA_COLUMN(UTrkDeltaRefGloParamSnp, uTrkDeltaRefGloParamSnp, float);
DECLARE_SOA_COLUMN(UTrkDeltaRefGloParamTgl, uTrkDeltaRefGloParamTgl, float);
DECLARE_SOA_COLUMN(UTrkDeltaRefGloParamQ2Pt, uTrkDeltaRefGloParamQ2Pt, float);
DECLARE_SOA_COLUMN(UTrkDeltaTOFdX, uTrkDeltaTOFdX, float);
DECLARE_SOA_COLUMN(UTrkDeltaTOFdZ, uTrkDeltaTOFdZ, float);

DECLARE_SOA_COLUMN(UTrkBeta, uTrkBeta, float);
DECLARE_SOA_COLUMN(UTrkBetaerror, uTrkBetaerror, float);
DECLARE_SOA_COLUMN(UTrkMass, uTrkMass, float);

DECLARE_SOA_COLUMN(UTrkCTpcNSigmaPi, uTrkCTpcNSigmaPi, uint8_t);
DECLARE_SOA_COLUMN(UTrkCTofNSigmaPi, uTrkCTofNSigmaPi, uint8_t);
DECLARE_SOA_COLUMN(UTrkCTpcNSigmaKa, uTrkCTpcNSigmaKa, uint8_t);
DECLARE_SOA_COLUMN(UTrkCTofNSigmaKa, uTrkCTofNSigmaKa, uint8_t);
DECLARE_SOA_COLUMN(UTrkCTpcNSigmaPr, uTrkCTpcNSigmaPr, uint8_t);
DECLARE_SOA_COLUMN(UTrkCTofNSigmaPr, uTrkCTofNSigmaPr, uint8_t);
DECLARE_SOA_COLUMN(UTrkCTpcNSigmaEl, uTrkCTpcNSigmaEl, uint8_t);
DECLARE_SOA_COLUMN(UTrkCTofNSigmaEl, uTrkCTofNSigmaEl, uint8_t);
DECLARE_SOA_COLUMN(UTrkCTpcNSigmaDe, uTrkCTpcNSigmaDe, uint8_t);
DECLARE_SOA_COLUMN(UTrkCTofNSigmaDe, uTrkCTofNSigmaDe, uint8_t);

DECLARE_SOA_COLUMN(UTrkTofSignal, uTrkTofSignal, float);

DECLARE_SOA_COLUMN(UTrkPx, uTrkPx, float);
DECLARE_SOA_COLUMN(UTrkPy, uTrkPy, float);
DECLARE_SOA_COLUMN(UTrkPz, uTrkPz, float);
DECLARE_SOA_COLUMN(UTrkPt, uTrkPt, float);
DECLARE_SOA_COLUMN(UTrkEta, uTrkEta, float);
DECLARE_SOA_COLUMN(UTrkPhi, uTrkPhi, float);
DECLARE_SOA_COLUMN(UTrkQPt, uTrkQPt, float);
DECLARE_SOA_COLUMN(UTrkTgl, uTrkTgl, float);
DECLARE_SOA_COLUMN(UTrkDcaCalc, uTrkDcaCalc, float);
DECLARE_SOA_COLUMN(UTrkDcaXY, uTrkDcaXY, float);
DECLARE_SOA_COLUMN(UTrkDcaZ, uTrkDcaZ, float);
DECLARE_SOA_COLUMN(UTrkSigmaDcaXY2, uTrkSigmaDcaXY2, float);
DECLARE_SOA_COLUMN(UTrkSigmaDcaZ2, uTrkSigmaDcaZ2, float);
DECLARE_SOA_COLUMN(UTrkAlpha, uTrkAlpha, float);

DECLARE_SOA_COLUMN(LTrkX, lTrkX, float);
DECLARE_SOA_COLUMN(LTrkY, lTrkY, float);
DECLARE_SOA_COLUMN(LTrkZ, lTrkZ, float);
DECLARE_SOA_COLUMN(LTrkSnp, lTrkSnp, float);
DECLARE_SOA_COLUMN(LTrkSigned1Pt, lTrkSigned1Pt, float);
DECLARE_SOA_COLUMN(LTrkP, lTrkP, float);
DECLARE_SOA_COLUMN(LTrkIsWithinBeamPipe, lTrkIsWithinBeamPipe, uint8_t);

DECLARE_SOA_COLUMN(LTrkSigmaY, lTrkSigmaY, float);
DECLARE_SOA_COLUMN(LTrkSigmaZ, lTrkSigmaZ, float);
DECLARE_SOA_COLUMN(LTrkSigmaSnp, lTrkSigmaSnp, float);
DECLARE_SOA_COLUMN(LTrkSigmaTgl, lTrkSigmaTgl, float);
DECLARE_SOA_COLUMN(LTrkSigma1Pt, lTrkSigma1Pt, float);

DECLARE_SOA_COLUMN(LTrkRhoZY, lTrkRhoZY, float);
DECLARE_SOA_COLUMN(LTrkRhoSnpY, lTrkRhoSnpY, float);
DECLARE_SOA_COLUMN(LTrkRhoSnpZ, lTrkRhoSnpZ, float);
DECLARE_SOA_COLUMN(LTrkRhoTglY, lTrkRhoTglY, float);
DECLARE_SOA_COLUMN(LTrkRhoTglZ, lTrkRhoTglZ, float);
DECLARE_SOA_COLUMN(LTrkRhoTglSnp, lTrkRhoTglSnp, float);
DECLARE_SOA_COLUMN(LTrkRho1PtY, lTrkRho1PtY, float);
DECLARE_SOA_COLUMN(LTrkRho1PtZ, lTrkRho1PtZ, float);
DECLARE_SOA_COLUMN(LTrkRho1PtSnp, lTrkRho1PtSnp, float);
DECLARE_SOA_COLUMN(LTrkRho1PtTgl, lTrkRho1PtTgl, float);

DECLARE_SOA_COLUMN(LTrkCYY, lTrkCYY, float);
DECLARE_SOA_COLUMN(LTrkCZY, lTrkCZY, float);
DECLARE_SOA_COLUMN(LTrkCZZ, lTrkCZZ, float);
DECLARE_SOA_COLUMN(LTrkCSnpY, lTrkCSnpY, float);
DECLARE_SOA_COLUMN(LTrkCSnpZ, lTrkCSnpZ, float);
DECLARE_SOA_COLUMN(LTrkCSnpSnp, lTrkCSnpSnp, float);
DECLARE_SOA_COLUMN(LTrkCTglY, lTrkCTglY, float);
DECLARE_SOA_COLUMN(LTrkCTglZ, lTrkCTglZ, float);
DECLARE_SOA_COLUMN(LTrkCTglSnp, lTrkCTglSnp, float);
DECLARE_SOA_COLUMN(LTrkCTglTgl, lTrkCTglTgl, float);
DECLARE_SOA_COLUMN(LTrkC1PtY, lTrkC1PtY, float);
DECLARE_SOA_COLUMN(LTrkC1PtZ, lTrkC1PtZ, float);
DECLARE_SOA_COLUMN(LTrkC1PtSnp, lTrkC1PtSnp, float);
DECLARE_SOA_COLUMN(LTrkC1PtTgl, lTrkC1PtTgl, float);
DECLARE_SOA_COLUMN(LTrkC1Pt21Pt2, lTrkC1Pt21Pt2, float);

DECLARE_SOA_COLUMN(LTrkTpcInnerParam, lTrkTpcInnerParam, float); // if structured, adjust type
DECLARE_SOA_COLUMN(LTrkFlags, lTrkFlags, uint32_t);
DECLARE_SOA_COLUMN(LTrkItsClusterSizes, lTrkItsClusterSizes, uint32_t);
DECLARE_SOA_COLUMN(LTrkTpcNClsFindable, lTrkTpcNClsFindable, uint8_t);
DECLARE_SOA_COLUMN(LTrkTpcNClsFindableMinusFound, lTrkTpcNClsFindableMinusFound, int8_t);
DECLARE_SOA_COLUMN(LTrkTpcNClsFindableMinusCrossedRows, lTrkTpcNClsFindableMinusCrossedRows, int8_t);
DECLARE_SOA_COLUMN(LTrkTpcNClsShared, lTrkTpcNClsShared, uint8_t);
DECLARE_SOA_COLUMN(LTrkTrdPattern, lTrkTrdPattern, uint8_t);

DECLARE_SOA_COLUMN(LTrkDetectorMask, lTrkDetectorMask, uint8_t);

DECLARE_SOA_COLUMN(LTrkItsChi2NCl, lTrkItsChi2NCl, float);
DECLARE_SOA_COLUMN(LTrkTpcChi2NCl, lTrkTpcChi2NCl, float);
DECLARE_SOA_COLUMN(LTrkTrdChi2, lTrkTrdChi2, float);
DECLARE_SOA_COLUMN(LTrkTofChi2, lTrkTofChi2, float);
DECLARE_SOA_COLUMN(LTrkTpcSignal, lTrkTpcSignal, float);
DECLARE_SOA_COLUMN(LTrkTrdSignal, lTrkTrdSignal, float);
DECLARE_SOA_COLUMN(LTrkLength, lTrkLength, float);
DECLARE_SOA_COLUMN(LTrkTofExpMom, lTrkTofExpMom, float);
DECLARE_SOA_COLUMN(LTrkEtaEmcal, lTrkEtaEmcal, float);
DECLARE_SOA_COLUMN(LTrkPhiEmcal, lTrkPhiEmcal, float);
DECLARE_SOA_COLUMN(LTrkTime, lTrkTime, float);
DECLARE_SOA_COLUMN(LTrkTimeRes, lTrkTimeRes, float);
DECLARE_SOA_COLUMN(LTrkPidForTracking, lTrkPidForTracking, int8_t);
DECLARE_SOA_COLUMN(LTrkIsPVContributor, lTrkIsPVContributor, uint8_t);
DECLARE_SOA_COLUMN(LTrkTpcNClsFound, lTrkTPCNClsFound, uint8_t);
DECLARE_SOA_COLUMN(LTrkTpcNClsCrossedRows, lTrkTPCNClsCrossedRows, uint8_t);
DECLARE_SOA_COLUMN(LTrkItsClusterMap, lTrkITSClusterMap, uint8_t);
DECLARE_SOA_COLUMN(LTrkItsNCls, lTrkITSNCls, uint8_t);
DECLARE_SOA_COLUMN(LTrkItsNClsInnerBarrel, lTrkITSNClsInnerBarrel, uint8_t);
DECLARE_SOA_COLUMN(LTrkTpcCrossedRowsOverFindableCls, lTrkTPCCrossedRowsOverFindableCls, float);
DECLARE_SOA_COLUMN(LTrkTpcFoundOverFindableCls, lTrkTPCFoundOverFindableCls, float);
DECLARE_SOA_COLUMN(LTrkTpcFractionSharedCls, lTrkTPCFractionSharedCls, float);
DECLARE_SOA_COLUMN(LTrkDetectorMap, lTrkDetectorMap, uint32_t);

DECLARE_SOA_COLUMN(LTrkTpcTime0, lTrkTpcTime0, float);
DECLARE_SOA_COLUMN(LTrkTpcDcaR, lTrkTpcDcaR, float);
DECLARE_SOA_COLUMN(LTrkTpcDcaZ, lTrkTpcDcaZ, float);
DECLARE_SOA_COLUMN(LTrkTpcClusterByteMask, lTrkTpcClusterByteMask, uint8_t);

DECLARE_SOA_COLUMN(LTrkTpcdEdxMax0R, lTrkTpcdEdxMax0R, float);
DECLARE_SOA_COLUMN(LTrkTpcdEdxMax1R, lTrkTpcdEdxMax1R, float);
DECLARE_SOA_COLUMN(LTrkTpcdEdxMax2R, lTrkTpcdEdxMax2R, float);
DECLARE_SOA_COLUMN(LTrkTpcdEdxMax3R, lTrkTpcdEdxMax3R, float);
DECLARE_SOA_COLUMN(LTrkTpcdEdxTot0R, lTrkTpcdEdxTot0R, float);
DECLARE_SOA_COLUMN(LTrkTpcdEdxTot1R, lTrkTpcdEdxTot1R, float);
DECLARE_SOA_COLUMN(LTrkTpcdEdxTot2R, lTrkTpcdEdxTot2R, float);
DECLARE_SOA_COLUMN(LTrkTpcdEdxTot3R, lTrkTpcdEdxTot3R, float);

DECLARE_SOA_COLUMN(LTrkDeltaRefContParamY, lTrkDeltaRefContParamY, float);
DECLARE_SOA_COLUMN(LTrkDeltaRefContParamZ, lTrkDeltaRefContParamZ, float);
DECLARE_SOA_COLUMN(LTrkDeltaRefContParamSnp, lTrkDeltaRefContParamSnp, float);
DECLARE_SOA_COLUMN(LTrkDeltaRefContParamTgl, lTrkDeltaRefContParamTgl, float);
DECLARE_SOA_COLUMN(LTrkDeltaRefContParamQ2Pt, lTrkDeltaRefContParamQ2Pt, float);
DECLARE_SOA_COLUMN(LTrkDeltaRefGloParamY, lTrkDeltaRefGloParamY, float);
DECLARE_SOA_COLUMN(LTrkDeltaRefGloParamZ, lTrkDeltaRefGloParamZ, float);
DECLARE_SOA_COLUMN(LTrkDeltaRefGloParamSnp, lTrkDeltaRefGloParamSnp, float);
DECLARE_SOA_COLUMN(LTrkDeltaRefGloParamTgl, lTrkDeltaRefGloParamTgl, float);
DECLARE_SOA_COLUMN(LTrkDeltaRefGloParamQ2Pt, lTrkDeltaRefGloParamQ2Pt, float);
DECLARE_SOA_COLUMN(LTrkDeltaTOFdX, lTrkDeltaTOFdX, float);
DECLARE_SOA_COLUMN(LTrkDeltaTOFdZ, lTrkDeltaTOFdZ, float);

DECLARE_SOA_COLUMN(LTrkBeta, lTrkBeta, float);
DECLARE_SOA_COLUMN(LTrkBetaerror, lTrkBetaerror, float);
DECLARE_SOA_COLUMN(LTrkMass, lTrkMass, float);

DECLARE_SOA_COLUMN(LTrkCTpcNSigmaPi, lTrkCTpcNSigmaPi, uint8_t);
DECLARE_SOA_COLUMN(LTrkCTofNSigmaPi, lTrkCTofNSigmaPi, uint8_t);
DECLARE_SOA_COLUMN(LTrkCTpcNSigmaKa, lTrkCTpcNSigmaKa, uint8_t);
DECLARE_SOA_COLUMN(LTrkCTofNSigmaKa, lTrkCTofNSigmaKa, uint8_t);
DECLARE_SOA_COLUMN(LTrkCTpcNSigmaPr, lTrkCTpcNSigmaPr, uint8_t);
DECLARE_SOA_COLUMN(LTrkCTofNSigmaPr, lTrkCTofNSigmaPr, uint8_t);
DECLARE_SOA_COLUMN(LTrkCTpcNSigmaEl, lTrkCTpcNSigmaEl, uint8_t);
DECLARE_SOA_COLUMN(LTrkCTofNSigmaEl, lTrkCTofNSigmaEl, uint8_t);
DECLARE_SOA_COLUMN(LTrkCTpcNSigmaDe, lTrkCTpcNSigmaDe, uint8_t);
DECLARE_SOA_COLUMN(LTrkCTofNSigmaDe, lTrkCTofNSigmaDe, uint8_t);

DECLARE_SOA_COLUMN(LTrkTofSignal, lTrkTofSignal, float);

DECLARE_SOA_COLUMN(LTrkPx, lTrkPx, float);
DECLARE_SOA_COLUMN(LTrkPy, lTrkPy, float);
DECLARE_SOA_COLUMN(LTrkPz, lTrkPz, float);
DECLARE_SOA_COLUMN(LTrkPt, lTrkPt, float);
DECLARE_SOA_COLUMN(LTrkEta, lTrkEta, float);
DECLARE_SOA_COLUMN(LTrkPhi, lTrkPhi, float);
DECLARE_SOA_COLUMN(LTrkQPt, lTrkQPt, float);
DECLARE_SOA_COLUMN(LTrkTgl, lTrkTgl, float);

DECLARE_SOA_COLUMN(LTrkDcaCalc, lTrkDcaCalc, float);
DECLARE_SOA_COLUMN(LTrkDcaXY, lTrkDcaXY, float);
DECLARE_SOA_COLUMN(LTrkDcaZ, lTrkDcaZ, float);
DECLARE_SOA_COLUMN(LTrkSigmaDcaXY2, lTrkSigmaDcaXY2, float);
DECLARE_SOA_COLUMN(LTrkSigmaDcaZ2, lTrkSigmaDcaZ2, float);
DECLARE_SOA_COLUMN(LTrkAlpha, lTrkAlpha, float);

DECLARE_SOA_COLUMN(SumPtPair, sumPtPair, float);
DECLARE_SOA_COLUMN(SumQPtPair, sumQPtPair, float);
DECLARE_SOA_COLUMN(SumTglPair, sumTglPair, float);
DECLARE_SOA_COLUMN(SumDcaCalcPair, sumDcaCalcPair, float);
DECLARE_SOA_COLUMN(SumDcaXYPair, sumDcaXYPair, float);
DECLARE_SOA_COLUMN(SumAlphaPair, sumAlphaPair, float);
DECLARE_SOA_COLUMN(DiffPtPair, diffPtPair, float);
DECLARE_SOA_COLUMN(DiffQPtPair, diffQPtPair, float);
DECLARE_SOA_COLUMN(DiffTglPair, diffTglPair, float);
DECLARE_SOA_COLUMN(DiffDcaXYPair, diffDcaXYPair, float);
DECLARE_SOA_COLUMN(DiffAlphaPair, diffAlphaPair, float);
DECLARE_SOA_COLUMN(B_kG, b_kG, int8_t);

}; // namespace drcolmn

DECLARE_SOA_TABLE(DrCosmicPairs, "AOD", "DRCOSMICPAIRS", o2::soa::Index<>,
                  o2::aod::drcolmn::GlobalIndex,
                  o2::aod::drcolmn::OriginalIndex,
                  o2::aod::drcolmn::GDFId,

                  o2::aod::drcolmn::GUBCid,
                  o2::aod::drcolmn::GUCollId,

                  o2::aod::drcolmn::GUpperTrackId,
                  o2::aod::drcolmn::DrUpperTrackId,
                  o2::aod::drcolmn::UpperTrackId,

                  o2::aod::drcolmn::GLBCid,
                  o2::aod::drcolmn::GLCollId,

                  o2::aod::drcolmn::GLowerTrackId,
                  o2::aod::drcolmn::DrLowerTrackId,
                  o2::aod::drcolmn::LowerTrackId,

                  o2::aod::drcolmn::UTrkX,         //
                  o2::aod::drcolmn::UTrkAlpha,     //
                  o2::aod::drcolmn::UTrkY,         //
                  o2::aod::drcolmn::UTrkZ,         //
                  o2::aod::drcolmn::UTrkSnp,       //
                  o2::aod::drcolmn::UTrkTgl,       //
                  o2::aod::drcolmn::UTrkSigned1Pt, // (sign of charge)/Pt in c/GeV. Use pt and sign instead
                  o2::aod::drcolmn::UTrkPt,        // Transverse momentum of the track in GeV/c
                  o2::aod::drcolmn::UTrkP,         // Momentum in Gev/c
                  o2::aod::drcolmn::UTrkEta,       // Pseudorapidity
                  o2::aod::drcolmn::UTrkPhi,       // Phi of the track, in radians within [0, 2pi)

                  o2::aod::drcolmn::UTrkIsWithinBeamPipe, // Is the track within the beam pipe (= successfully propagated to a collision vertex)
                  o2::aod::drcolmn::UTrkPx,               // Momentum in x-direction in GeV/c
                  o2::aod::drcolmn::UTrkPy,               // Momentum in y-direction in GeV/c
                  o2::aod::drcolmn::UTrkPz,               // Momentum in z-direction in GeV/c

                  // // TracksCov
                  o2::aod::drcolmn::UTrkSigmaY,    // Covariance matrix
                  o2::aod::drcolmn::UTrkSigmaZ,    // Covariance matrix
                  o2::aod::drcolmn::UTrkSigmaSnp,  // Covariance matrix
                  o2::aod::drcolmn::UTrkSigmaTgl,  // Covariance matrix
                  o2::aod::drcolmn::UTrkSigma1Pt,  // Covariance matrix
                  o2::aod::drcolmn::UTrkRhoZY,     // Covariance matrix in compressed form
                  o2::aod::drcolmn::UTrkRhoSnpY,   // Covariance matrix in compressed form
                  o2::aod::drcolmn::UTrkRhoSnpZ,   // Covariance matrix in compressed form
                  o2::aod::drcolmn::UTrkRhoTglY,   // Covariance matrix in compressed form
                  o2::aod::drcolmn::UTrkRhoTglZ,   // Covariance matrix in compressed form
                  o2::aod::drcolmn::UTrkRhoTglSnp, // Covariance matrix in compressed form
                  o2::aod::drcolmn::UTrkRho1PtY,   // Covariance matrix in compressed form
                  o2::aod::drcolmn::UTrkRho1PtZ,   // Covariance matrix in compressed form
                  o2::aod::drcolmn::UTrkRho1PtSnp, // Covariance matrix in compressed form
                  o2::aod::drcolmn::UTrkRho1PtTgl, // Covariance matrix in compressed form
                  o2::aod::drcolmn::UTrkCYY,       //
                  o2::aod::drcolmn::UTrkCZY,
                  o2::aod::drcolmn::UTrkCZZ,
                  o2::aod::drcolmn::UTrkCSnpY,
                  o2::aod::drcolmn::UTrkCSnpZ,
                  o2::aod::drcolmn::UTrkCSnpSnp,
                  o2::aod::drcolmn::UTrkCTglY,
                  o2::aod::drcolmn::UTrkCTglZ,
                  o2::aod::drcolmn::UTrkCTglSnp,
                  o2::aod::drcolmn::UTrkCTglTgl,
                  o2::aod::drcolmn::UTrkC1PtY,
                  o2::aod::drcolmn::UTrkC1PtZ,
                  o2::aod::drcolmn::UTrkC1PtSnp,
                  o2::aod::drcolmn::UTrkC1PtTgl,
                  o2::aod::drcolmn::UTrkC1Pt21Pt2,

                  // // TracksExtra
                  o2::aod::drcolmn::UTrkTpcInnerParam,                   // Momentum at inner wall of the TPC
                  o2::aod::drcolmn::UTrkFlags,                           // Track flags. Run 2: see TrackFlagsRun2Enum | Run 3: see TrackFlags
                  o2::aod::drcolmn::UTrkItsClusterSizes,                 // Clusters sizes, four bits per a layer, starting from the innermost
                  o2::aod::drcolmn::UTrkTpcNClsFindable,                 // Findable TPC clusters for this track geometry
                  o2::aod::drcolmn::UTrkTpcNClsFindableMinusFound,       // TPC Clusters: Findable - Found
                  o2::aod::drcolmn::UTrkTpcNClsFindableMinusCrossedRows, // TPC Clusters: Findable - crossed rows
                  o2::aod::drcolmn::UTrkTpcNClsShared,                   // Number of shared TPC clusters
                  o2::aod::drcolmn::UTrkTrdPattern,                      // Contributor to the track on TRD layer in bits 0-5, starting from the innermost, bit 6 indicates a potentially split tracklet, bit 7 if the track crossed a padrow
                  o2::aod::drcolmn::UTrkItsChi2NCl,                      // Chi2 / cluster for the ITS track segment
                  o2::aod::drcolmn::UTrkTpcChi2NCl,                      // Chi2 / cluster for the TPC track segment
                  o2::aod::drcolmn::UTrkTrdChi2,                         // Chi2 for the TRD track segment
                  o2::aod::drcolmn::UTrkTofChi2,                         // Chi2 for the TOF track segment
                  o2::aod::drcolmn::UTrkTpcSignal,                       // dE/dx signal in the TPC
                  o2::aod::drcolmn::UTrkTrdSignal,                       // PID signal in the TRD
                  o2::aod::drcolmn::UTrkLength,                          // Track length
                  o2::aod::drcolmn::UTrkTofExpMom,                       // TOF expected momentum obtained in tracking, used to compute the expected times
                  o2::aod::drcolmn::UTrkEtaEmcal,                        //
                  o2::aod::drcolmn::UTrkPhiEmcal,                        //
                  o2::aod::drcolmn::UTrkTime,                            // Estimated time of the track in ns wrt collision.bc or ambiguousupperTrk.bcSlice[0]
                  o2::aod::drcolmn::UTrkTimeRes,                         // Resolution of the track time in ns (see TrackFlags::TrackTimeResIsRange)

                  o2::aod::drcolmn::UTrkPidForTracking,  // PID hypothesis used during tracking. See the constants in the class PID in PID.h
                  o2::aod::drcolmn::UTrkIsPVContributor, // Run 3: Has this track contributed to the collision vertex fit
                  o2::aod::drcolmn::UTrkDetectorMask,
                  o2::aod::drcolmn::UTrkTpcNClsFound,                  // Number of found TPC clusters
                  o2::aod::drcolmn::UTrkTpcNClsCrossedRows,            // Number of crossed TPC Rows
                  o2::aod::drcolmn::UTrkItsClusterMap,                 // ITS cluster map, one bit per a layer, starting from the innermost
                  o2::aod::drcolmn::UTrkItsNCls,                       // Number of ITS clusters
                  o2::aod::drcolmn::UTrkItsNClsInnerBarrel,            // Number of ITS clusters in the Inner Barrel
                  o2::aod::drcolmn::UTrkTpcCrossedRowsOverFindableCls, // Ratio crossed rows over findable clusters
                  o2::aod::drcolmn::UTrkTpcFoundOverFindableCls,       // Ratio of found over findable clusters
                  o2::aod::drcolmn::UTrkTpcFractionSharedCls,          // Fraction of shared TPC clusters
                  o2::aod::drcolmn::UTrkDetectorMap,                   // Detector map version 1, see enum DetectorMapEnum

                  // TracksQA
                  o2::aod::drcolmn::UTrkTpcTime0,           // tpc only time0 (mTime0 in TPC track)
                  o2::aod::drcolmn::UTrkTpcDcaR,            // tpc only DCAr
                  o2::aod::drcolmn::UTrkTpcDcaZ,            // tpc only DCAz
                  o2::aod::drcolmn::UTrkTpcClusterByteMask, // tracklet bitmask - track defining 8 tracklets (152=8*19 rows) bit set if nCluster>thr (default 5)
                  o2::aod::drcolmn::UTrkTpcdEdxMax0R,       // TPC dEdxQMax -ROC0/dEdx
                  o2::aod::drcolmn::UTrkTpcdEdxMax1R,       // TPC dEdxQMax -ROC1/dEdx
                  o2::aod::drcolmn::UTrkTpcdEdxMax2R,       // TPC dEdxQMax -ROC2/dEdx
                  o2::aod::drcolmn::UTrkTpcdEdxMax3R,       // TPC dEdxQMax -ROC3/dEdx
                  o2::aod::drcolmn::UTrkTpcdEdxTot0R,       // TPC dEdxQtot -ROC0/dEdx
                  o2::aod::drcolmn::UTrkTpcdEdxTot1R,       // TPC dEdxQtot -ROC1/dEdx
                  o2::aod::drcolmn::UTrkTpcdEdxTot2R,       // TPC dEdxQtot -ROC2/dEdx
                  o2::aod::drcolmn::UTrkTpcdEdxTot3R,       // TPC dEdxQtot -ROC3/dEdx
                  //
                //   o2::aod::drcolmn::UTrkBeta,
                //   o2::aod::drcolmn::UTrkBetaerror,
                //   o2::aod::drcolmn::UTrkMass,

                //   o2::aod::drcolmn::UTrkCTpcNSigmaPi, //, cfgPIDSigma0, cfgPIDSigma1, cfgPIDClampMin, cfgPIDClampMax)
                //   o2::aod::drcolmn::UTrkCTofNSigmaPi, //, cfgPIDSigma0, cfgPIDSigma1, cfgPIDClampMin, cfgPIDClampMax)
                //   o2::aod::drcolmn::UTrkCTpcNSigmaKa, //, cfgPIDSigma0, cfgPIDSigma1, cfgPIDClampMin, cfgPIDClampMax)
                //   o2::aod::drcolmn::UTrkCTofNSigmaKa, //, cfgPIDSigma0, cfgPIDSigma1, cfgPIDClampMin, cfgPIDClampMax)
                //   o2::aod::drcolmn::UTrkCTpcNSigmaPr, //, cfgPIDSigma0, cfgPIDSigma1, cfgPIDClampMin, cfgPIDClampMax)
                //   o2::aod::drcolmn::UTrkCTofNSigmaPr, //, cfgPIDSigma0, cfgPIDSigma1, cfgPIDClampMin, cfgPIDClampMax)
                //   o2::aod::drcolmn::UTrkCTpcNSigmaEl, //, cfgPIDSigma0, cfgPIDSigma1, cfgPIDClampMin, cfgPIDClampMax)
                //   o2::aod::drcolmn::UTrkCTofNSigmaEl, //, cfgPIDSigma0, cfgPIDSigma1, cfgPIDClampMin, cfgPIDClampMax)
                //   o2::aod::drcolmn::UTrkCTpcNSigmaDe, //, cfgPIDSigma0, cfgPIDSigma1, cfgPIDClampMin, cfgPIDClampMax)
                //   o2::aod::drcolmn::UTrkCTofNSigmaDe, //, cfgPIDSigma0, cfgPIDSigma1, cfgPIDClampMin, cfgPIDClampMax)

                  o2::aod::drcolmn::UTrkTofSignal,

                  o2::aod::drcolmn::UTrkDcaCalc ,       //	dcaXY	float	Impact parameter in XY of the track to the primary vertex                  
                  o2::aod::drcolmn::UTrkDcaXY,       //	dcaXY	float	Impact parameter in XY of the track to the primary vertex
                  o2::aod::drcolmn::UTrkDcaZ,        //	dcaZ	float	Impact parameter in Z of the track to the primary vertex
                  o2::aod::drcolmn::UTrkSigmaDcaXY2, //		sigmaDcaXY2	float	Impact parameter sigma^2 in XY of the track to the primary vertex
                  o2::aod::drcolmn::UTrkSigmaDcaZ2,  //		sigmaDcaZ2	float	Impact parameter sigma^2 in Z of the track to the primary verte

                  o2::aod::drcolmn::UTrkDeltaRefContParamY,
                  o2::aod::drcolmn::UTrkDeltaRefContParamZ,
                  o2::aod::drcolmn::UTrkDeltaRefContParamSnp,
                  o2::aod::drcolmn::UTrkDeltaRefContParamTgl,
                  o2::aod::drcolmn::UTrkDeltaRefContParamQ2Pt,
                  o2::aod::drcolmn::UTrkDeltaRefGloParamY,
                  o2::aod::drcolmn::UTrkDeltaRefGloParamZ,
                  o2::aod::drcolmn::UTrkDeltaRefGloParamSnp,
                  o2::aod::drcolmn::UTrkDeltaRefGloParamTgl,
                  o2::aod::drcolmn::UTrkDeltaRefGloParamQ2Pt,
                  o2::aod::drcolmn::UTrkDeltaTOFdX,
                  o2::aod::drcolmn::UTrkDeltaTOFdZ,

                  o2::aod::drcolmn::LTrkX,         //
                  o2::aod::drcolmn::LTrkAlpha,     //
                  o2::aod::drcolmn::LTrkY,         //
                  o2::aod::drcolmn::LTrkZ,         //
                  o2::aod::drcolmn::LTrkSnp,       //
                  o2::aod::drcolmn::LTrkTgl,       //
                  o2::aod::drcolmn::LTrkSigned1Pt, // (sign of charge)/Pt in c/GeV. Use pt and sign instead
                  o2::aod::drcolmn::LTrkPt,        // Transverse momentum of the track in GeV/c
                  o2::aod::drcolmn::LTrkP,         // Momentum in Gev/c
                  o2::aod::drcolmn::LTrkEta,       // Pseudorapidity
                  o2::aod::drcolmn::LTrkPhi,       // Phi of the track, in radians within [0, 2pi)

                  o2::aod::drcolmn::LTrkIsWithinBeamPipe, // Is the track within the beam pipe (= successfully propagated to a collision vertex)
                  o2::aod::drcolmn::LTrkPx,               // Momentum in x-direction in GeV/c
                  o2::aod::drcolmn::LTrkPy,               // Momentum in y-direction in GeV/c
                  o2::aod::drcolmn::LTrkPz,               // Momentum in z-direction in GeV/c

                  // // TracksCov
                  o2::aod::drcolmn::LTrkSigmaY,    // Covariance matrix
                  o2::aod::drcolmn::LTrkSigmaZ,    // Covariance matrix
                  o2::aod::drcolmn::LTrkSigmaSnp,  // Covariance matrix
                  o2::aod::drcolmn::LTrkSigmaTgl,  // Covariance matrix
                  o2::aod::drcolmn::LTrkSigma1Pt,  // Covariance matrix
                  o2::aod::drcolmn::LTrkRhoZY,     // Covariance matrix in compressed form
                  o2::aod::drcolmn::LTrkRhoSnpY,   // Covariance matrix in compressed form
                  o2::aod::drcolmn::LTrkRhoSnpZ,   // Covariance matrix in compressed form
                  o2::aod::drcolmn::LTrkRhoTglY,   // Covariance matrix in compressed form
                  o2::aod::drcolmn::LTrkRhoTglZ,   // Covariance matrix in compressed form
                  o2::aod::drcolmn::LTrkRhoTglSnp, // Covariance matrix in compressed form
                  o2::aod::drcolmn::LTrkRho1PtY,   // Covariance matrix in compressed form
                  o2::aod::drcolmn::LTrkRho1PtZ,   // Covariance matrix in compressed form
                  o2::aod::drcolmn::LTrkRho1PtSnp, // Covariance matrix in compressed form
                  o2::aod::drcolmn::LTrkRho1PtTgl, // Covariance matrix in compressed form
                  o2::aod::drcolmn::LTrkCYY,       //
                  o2::aod::drcolmn::LTrkCZY,
                  o2::aod::drcolmn::LTrkCZZ,
                  o2::aod::drcolmn::LTrkCSnpY,
                  o2::aod::drcolmn::LTrkCSnpZ,
                  o2::aod::drcolmn::LTrkCSnpSnp,
                  o2::aod::drcolmn::LTrkCTglY,
                  o2::aod::drcolmn::LTrkCTglZ,
                  o2::aod::drcolmn::LTrkCTglSnp,
                  o2::aod::drcolmn::LTrkCTglTgl,
                  o2::aod::drcolmn::LTrkC1PtY,
                  o2::aod::drcolmn::LTrkC1PtZ,
                  o2::aod::drcolmn::LTrkC1PtSnp,
                  o2::aod::drcolmn::LTrkC1PtTgl,
                  o2::aod::drcolmn::LTrkC1Pt21Pt2,

                  // // TracksExtra
                  o2::aod::drcolmn::LTrkTpcInnerParam,                   // Momentum at inner wall of the TPC
                  o2::aod::drcolmn::LTrkFlags,                           // Track flags. Run 2: see TrackFlagsRun2Enum | Run 3: see TrackFlags
                  o2::aod::drcolmn::LTrkItsClusterSizes,                 // Clusters sizes, four bits per a layer, starting from the innermost
                  o2::aod::drcolmn::LTrkTpcNClsFindable,                 // Findable TPC clusters for this track geometry
                  o2::aod::drcolmn::LTrkTpcNClsFindableMinusFound,       // TPC Clusters: Findable - Found
                  o2::aod::drcolmn::LTrkTpcNClsFindableMinusCrossedRows, // TPC Clusters: Findable - crossed rows
                  o2::aod::drcolmn::LTrkTpcNClsShared,                   // Number of shared TPC clusters
                  o2::aod::drcolmn::LTrkTrdPattern,                      // Contributor to the track on TRD layer in bits 0-5, starting from the innermost, bit 6 indicates a potentially split tracklet, bit 7 if the track crossed a padrow
                  o2::aod::drcolmn::LTrkItsChi2NCl,                      // Chi2 / cluster for the ITS track segment
                  o2::aod::drcolmn::LTrkTpcChi2NCl,                      // Chi2 / cluster for the TPC track segment
                  o2::aod::drcolmn::LTrkTrdChi2,                         // Chi2 for the TRD track segment
                  o2::aod::drcolmn::LTrkTofChi2,                         // Chi2 for the TOF track segment
                  o2::aod::drcolmn::LTrkTpcSignal,                       // dE/dx signal in the TPC
                  o2::aod::drcolmn::LTrkTrdSignal,                       // PID signal in the TRD
                  o2::aod::drcolmn::LTrkLength,                          // Track length
                  o2::aod::drcolmn::LTrkTofExpMom,                       // TOF expected momentum obtained in tracking, used to compute the expected times
                  o2::aod::drcolmn::LTrkEtaEmcal,                        //
                  o2::aod::drcolmn::LTrkPhiEmcal,                        //
                  o2::aod::drcolmn::LTrkTime,                            // Estimated time of the track in ns wrt collision.bc or ambiguousupperTrk.bcSlice[0]
                  o2::aod::drcolmn::LTrkTimeRes,                         // Resolution of the track time in ns (see TrackFlags::TrackTimeResIsRange)

                  o2::aod::drcolmn::LTrkPidForTracking,  // PID hypothesis used during tracking. See the constants in the class PID in PID.h
                  o2::aod::drcolmn::LTrkIsPVContributor, // Run 3: Has this track contributed to the collision vertex fit
                  o2::aod::drcolmn::LTrkDetectorMask,
                  o2::aod::drcolmn::LTrkTpcNClsFound,                  // Number of found TPC clusters
                  o2::aod::drcolmn::LTrkTpcNClsCrossedRows,            // Number of crossed TPC Rows
                  o2::aod::drcolmn::LTrkItsClusterMap,                 // ITS cluster map, one bit per a layer, starting from the innermost
                  o2::aod::drcolmn::LTrkItsNCls,                       // Number of ITS clusters
                  o2::aod::drcolmn::LTrkItsNClsInnerBarrel,            // Number of ITS clusters in the Inner Barrel
                  o2::aod::drcolmn::LTrkTpcCrossedRowsOverFindableCls, // Ratio crossed rows over findable clusters
                  o2::aod::drcolmn::LTrkTpcFoundOverFindableCls,       // Ratio of found over findable clusters
                  o2::aod::drcolmn::LTrkTpcFractionSharedCls,          // Fraction of shared TPC clusters
                  o2::aod::drcolmn::LTrkDetectorMap,                   // Detector map version 1, see enum DetectorMapEnum

                  // TracksQA
                  o2::aod::drcolmn::LTrkTpcTime0,           // tpc only time0 (mTime0 in TPC track)
                  o2::aod::drcolmn::LTrkTpcDcaR,            // tpc only DCAr
                  o2::aod::drcolmn::LTrkTpcDcaZ,            // tpc only DCAz
                  o2::aod::drcolmn::LTrkTpcClusterByteMask, // tracklet bitmask - track defining 8 tracklets (152=8*19 rows) bit set if nCluster>thr (default 5)
                  o2::aod::drcolmn::LTrkTpcdEdxMax0R,       // TPC dEdxQMax -ROC0/dEdx
                  o2::aod::drcolmn::LTrkTpcdEdxMax1R,       // TPC dEdxQMax -ROC1/dEdx
                  o2::aod::drcolmn::LTrkTpcdEdxMax2R,       // TPC dEdxQMax -ROC2/dEdx
                  o2::aod::drcolmn::LTrkTpcdEdxMax3R,       // TPC dEdxQMax -ROC3/dEdx
                  o2::aod::drcolmn::LTrkTpcdEdxTot0R,       // TPC dEdxQtot -ROC0/dEdx
                  o2::aod::drcolmn::LTrkTpcdEdxTot1R,       // TPC dEdxQtot -ROC1/dEdx
                  o2::aod::drcolmn::LTrkTpcdEdxTot2R,       // TPC dEdxQtot -ROC2/dEdx
                  o2::aod::drcolmn::LTrkTpcdEdxTot3R,       // TPC dEdxQtot -ROC3/dEdx
                  //
                //   o2::aod::drcolmn::LTrkBeta,
                //   o2::aod::drcolmn::LTrkBetaerror,
                //   o2::aod::drcolmn::LTrkMass,

                //   o2::aod::drcolmn::LTrkCTpcNSigmaPi, //, cfgPIDSigma0, cfgPIDSigma1, cfgPIDClampMin, cfgPIDClampMax)
                //   o2::aod::drcolmn::LTrkCTofNSigmaPi, //, cfgPIDSigma0, cfgPIDSigma1, cfgPIDClampMin, cfgPIDClampMax)
                //   o2::aod::drcolmn::LTrkCTpcNSigmaKa, //, cfgPIDSigma0, cfgPIDSigma1, cfgPIDClampMin, cfgPIDClampMax)
                //   o2::aod::drcolmn::LTrkCTofNSigmaKa, //, cfgPIDSigma0, cfgPIDSigma1, cfgPIDClampMin, cfgPIDClampMax)
                //   o2::aod::drcolmn::LTrkCTpcNSigmaPr, //, cfgPIDSigma0, cfgPIDSigma1, cfgPIDClampMin, cfgPIDClampMax)
                //   o2::aod::drcolmn::LTrkCTofNSigmaPr, //, cfgPIDSigma0, cfgPIDSigma1, cfgPIDClampMin, cfgPIDClampMax)
                //   o2::aod::drcolmn::LTrkCTpcNSigmaEl, //, cfgPIDSigma0, cfgPIDSigma1, cfgPIDClampMin, cfgPIDClampMax)
                //   o2::aod::drcolmn::LTrkCTofNSigmaEl, //, cfgPIDSigma0, cfgPIDSigma1, cfgPIDClampMin, cfgPIDClampMax)
                //   o2::aod::drcolmn::LTrkCTpcNSigmaDe, //, cfgPIDSigma0, cfgPIDSigma1, cfgPIDClampMin, cfgPIDClampMax)
                //   o2::aod::drcolmn::LTrkCTofNSigmaDe, //, cfgPIDSigma0, cfgPIDSigma1, cfgPIDClampMin, cfgPIDClampMax)

                  o2::aod::drcolmn::LTrkTofSignal,

                  o2::aod::drcolmn::LTrkDcaCalc ,       //                
                  o2::aod::drcolmn::LTrkDcaXY,       //	dcaXY	float	Impact parameter in XY of the track to the primary vertex
                  o2::aod::drcolmn::LTrkDcaZ,        //	dcaZ	float	Impact parameter in Z of the track to the primary vertex
                  o2::aod::drcolmn::LTrkSigmaDcaXY2, //		sigmaDcaXY2	float	Impact parameter sigma^2 in XY of the track to the primary vertex
                  o2::aod::drcolmn::LTrkSigmaDcaZ2,  //		sigmaDcaZ2	float	Impact parameter sigma^2 in Z of the track to the primary verte

                  o2::aod::drcolmn::LTrkDeltaRefContParamY,
                  o2::aod::drcolmn::LTrkDeltaRefContParamZ,
                  o2::aod::drcolmn::LTrkDeltaRefContParamSnp,
                  o2::aod::drcolmn::LTrkDeltaRefContParamTgl,
                  o2::aod::drcolmn::LTrkDeltaRefContParamQ2Pt,
                  o2::aod::drcolmn::LTrkDeltaRefGloParamY,
                  o2::aod::drcolmn::LTrkDeltaRefGloParamZ,
                  o2::aod::drcolmn::LTrkDeltaRefGloParamSnp,
                  o2::aod::drcolmn::LTrkDeltaRefGloParamTgl,
                  o2::aod::drcolmn::LTrkDeltaRefGloParamQ2Pt,
                  o2::aod::drcolmn::LTrkDeltaTOFdX,
                  o2::aod::drcolmn::LTrkDeltaTOFdZ,

                  o2::aod::drcolmn::SumPtPair,
                  o2::aod::drcolmn::SumQPtPair,
                  o2::aod::drcolmn::SumTglPair,
                  o2::aod::drcolmn::SumDcaCalcPair,
                  o2::aod::drcolmn::SumDcaXYPair,
                  o2::aod::drcolmn::SumAlphaPair,

                  o2::aod::drcolmn::DiffPtPair,
                  o2::aod::drcolmn::DiffQPtPair,
                  o2::aod::drcolmn::DiffTglPair,
                  o2::aod::drcolmn::DiffDcaXYPair,
                  o2::aod::drcolmn::DiffAlphaPair,

                  o2::aod::drcolmn::UTimestamp,
                  o2::aod::drcolmn::UTfId,
                  o2::aod::drcolmn::UBCinTF,

                  o2::aod::drcolmn::LTimestamp,
                  o2::aod::drcolmn::LTfId,
                  o2::aod::drcolmn::LBCinTF,
                  o2::aod::drcolmn::B_kG
                );
} // namespace o2::aod
