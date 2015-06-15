// framework include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

// L1IT include files
#include "L1Trigger/L1IntegratedMuonTrigger/interface/TriggerPrimitive.h"
#include "L1Trigger/L1IntegratedMuonTrigger/interface/TriggerPrimitiveFwd.h"

#include "L1Trigger/L1IntegratedMuonTrigger/interface/MBLTCollection.h"
#include "L1Trigger/L1IntegratedMuonTrigger/interface/MBLTCollectionFwd.h"

#include "DataFormats/MuonDetId/interface/RPCDetId.h"
#include "DataFormats/MuonDetId/interface/DTChamberId.h"

#include "DataFormats/L1DTTrackFinder/interface/L1MuDTChambPhDigi.h"
#include "DataFormats/L1DTTrackFinder/interface/L1MuDTChambPhContainer.h"

#include "DataFormats/L1DTTrackFinder/interface/L1MuDTChambThDigi.h"
#include "DataFormats/L1DTTrackFinder/interface/L1MuDTChambThContainer.h"

#include "Geometry/HcalTowerAlgo/interface/HcalGeometry.h"
#include "Geometry/Records/interface/HcalGeometryRecord.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "L1Trigger/L1IntegratedMuonTrigger/interface/PrimitiveCombiner.h"

// user include files
#include "DataFormats/Math/interface/deltaPhi.h"
bool doDebug = false;

class L1ITMuonBarrelPrimitiveProducer : public edm::EDProducer {

public:
  L1ITMuonBarrelPrimitiveProducer(const edm::ParameterSet&);
  ~L1ITMuonBarrelPrimitiveProducer();
  virtual void produce(edm::Event&, const edm::EventSetup&);

private:
  edm::InputTag _mbltCollectionInput;
  const L1ITMu::PrimitiveCombiner::resolutions _resol;
  const int _qualityRemappingMode;
  const int _qualityHORemappingMode;
  const int _useRpcBxForDtBelowQuality;
  const bool _is7QualityCodes;
  const bool _isHOQualityCodes;
  const bool _addHitsFromHO;
  double _minDistForRpcHOMatch;
  const double _minDistForHORpcClusterMatch;
  edm::ESHandle<DTGeometry> _muonGeom;
  edm::ESHandle<CaloGeometry> _hoGeom;
};

std::ostream & operator<< (std::ostream & out, const L1ITMu::TriggerPrimitiveList & rpc )
{
  std::vector<L1ITMu::TriggerPrimitiveRef>::const_iterator it = rpc.begin();
  std::vector<L1ITMu::TriggerPrimitiveRef>::const_iterator end = rpc.end();
  for ( ; it != end ; ++it ) out << (*it)->getCMSGlobalPhi() << '\t';
  out << std::endl;
  return out;
}

L1ITMuonBarrelPrimitiveProducer::~L1ITMuonBarrelPrimitiveProducer()
{
}

L1ITMuonBarrelPrimitiveProducer::L1ITMuonBarrelPrimitiveProducer( const edm::ParameterSet& iConfig )
  : _mbltCollectionInput( iConfig.getParameter<edm::InputTag>("MBLTCollection") ),
    _resol( iConfig.getParameter<double>("xDtResol"), 
	    iConfig.getParameter<double>("xRpcResol"),
	    iConfig.getParameter<double>("phibDtCorrResol"),
	    iConfig.getParameter<double>("phibDtUnCorrResol"),
	    iConfig.getParameter<double>("xHOResol")
	    ),
    _qualityRemappingMode( iConfig.getParameter<int>("qualityRemappingMode") ),
    _qualityHORemappingMode( iConfig.getParameter<int>("qualityHORemappingMode") ),
    _useRpcBxForDtBelowQuality( iConfig.getParameter<int>("useRpcBxForDtBelowQuality") ),
    _is7QualityCodes( iConfig.getParameter<bool>("is7QualityCodes") ),
    _isHOQualityCodes( iConfig.getParameter<bool>("isHOQualityCodes") ),
    _addHitsFromHO( iConfig.getUntrackedParameter<bool>("addHitsFromHO", false) ),
    _minDistForRpcHOMatch( iConfig.getParameter<double>("minDistForRpcHOMatch") ),
    _minDistForHORpcClusterMatch(iConfig.getParameter<double>("minDistForHORpcClusterMatch") )

{
  produces<L1MuDTChambPhContainer>();
  //produces<L1MuDTChambThContainer>();
  // produces<std::vector<L1MuDTChambPhDigi> >();
}


void L1ITMuonBarrelPrimitiveProducer::produce( edm::Event& iEvent, 
					      const edm::EventSetup& iSetup )
{
  //  std::cout<<"Event Id---------------------------------->"<< iEvent.id().event() << std::endl;
  iSetup.get<MuonGeometryRecord>().get(_muonGeom);
  iSetup.get<CaloGeometryRecord>().get(_hoGeom);

  std::auto_ptr<L1MuDTChambPhContainer> out(new L1MuDTChambPhContainer);
  std::vector<L1MuDTChambPhDigi> phiChambVector;
  
  edm::Handle<L1ITMu::MBLTContainer> mbltContainer;
  iEvent.getByLabel( _mbltCollectionInput, mbltContainer );

  L1ITMu::MBLTContainer::const_iterator st = mbltContainer->begin();
  L1ITMu::MBLTContainer::const_iterator stend = mbltContainer->end();

  L1MuDTChambPhContainer phiContainer;
  std::vector<L1MuDTChambPhDigi> phiVector;
  
  for ( ; st != stend; ++st ) {

    const L1ITMu::MBLTCollection & mbltStation = st->second;

    /// useful index
    int station = mbltStation.station();
    int wheel   = mbltStation.wheel();
    int sector  = mbltStation.sector();
    ///

    /// get dt to rpc associations
    size_t dtListSize = mbltStation.getDtSegments().size();
    std::vector<size_t> uncorrelated;
    std::vector<size_t> correlated;

    if(doDebug)  std::cout<<" ################# DtSegments size --> "<< dtListSize <<" for (Wh,Se,St) -------> ("<< wheel<< ", "<< sector<<", "<< station <<")"<< std::endl;
    
    
    for ( size_t iDt = 0; iDt < dtListSize; ++iDt ) {
      const L1ITMu::TriggerPrimitiveRef & dt = mbltStation.getDtSegments().at(iDt);
      int dtquality = dt->getDTData().qualityCode;

      switch ( dtquality ) {
      case -1 : continue;/// -1 are theta
      case 0 : /* qualityCode = -2;*/ break;
      case 1 : /* qualityCode = -2;*/ break;
      case 2 : uncorrelated.push_back( iDt ); continue;  // HI
      case 3 : uncorrelated.push_back( iDt ); continue;  // HO
      case 4 : correlated.push_back( iDt ); continue;    // LL
      case 5 : correlated.push_back( iDt ); continue;    // HL
      case 6 : correlated.push_back( iDt ); continue;    // HH
      default : /* qualityCode = dtquality; */ break;
      }
    }
    
    if(doDebug)  std::cout<<"SIZE OF THE (cor, uncor) --> ("<<correlated.size() <<", "<< uncorrelated.size()<<")"<< std::endl; 
    if(doDebug)  std::cout<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ Correlated with size :: "<<correlated.size() << std::endl;
    
    // START OF BX ANALYSIS FOR CORRELATED TRIGGER
    size_t cSize = correlated.size();
    for ( size_t idxDt = 0; idxDt < cSize; ++idxDt ) {
      int bx=-999;
      int iDt = correlated.at(idxDt);
      if ( iDt < 0 ) continue;
      const L1ITMu::TriggerPrimitive & dt = *mbltStation.getDtSegments().at(iDt);
 
      L1ITMu::TriggerPrimitiveList rpcInMatch  = mbltStation.getRpcInAssociatedStubs( iDt );
      L1ITMu::TriggerPrimitiveList rpcOutMatch = mbltStation.getRpcOutAssociatedStubs( iDt );
      L1ITMu::TriggerPrimitiveList hoMatch;
      if(_addHitsFromHO) hoMatch    = mbltStation.getHOAssociatedStubs( iDt);

      size_t rpcInMatchSize  = rpcInMatch.size();
      size_t rpcOutMatchSize = rpcOutMatch.size();
      size_t hoMatchSize     = 0;
      if(_addHitsFromHO) hoMatchSize =  hoMatch.size();
      
      if(_addHitsFromHO && doDebug)  
	std::cout<<"MatchSize (rpcIn, rpcOut, ho) : (" << rpcInMatchSize <<", "<< rpcOutMatchSize << ", " << hoMatchSize <<")"<< ", DtSt -> "<< dt.getDTData().station <<  std::endl;
      else if( doDebug)  
	std::cout<<"MatchSize (rpcIn, rpcOut) : (" << rpcInMatchSize <<", "<< rpcOutMatchSize << ")"<<  std::endl;

      int qualityCode = 6;      

      if ( rpcInMatchSize && rpcOutMatchSize ) {
	const L1ITMu::TriggerPrimitive & rpcIn = *rpcInMatch.front();
	const L1ITMu::TriggerPrimitive & rpcOut = *rpcOutMatch.front();

	/// only the first is real...
	// LG try also to reassign BX to single H using RPC BX, e.g. do not ask for DT and RPC to have the same BX
	if ( ( dt.getBX() == rpcIn.getBX() && dt.getBX() == rpcOut.getBX() )
	     || (_qualityRemappingMode > 1 && rpcIn.getBX()==rpcOut.getBX() && abs(dt.getBX()-rpcIn.getBX())<=1) ) {
	  bx = rpcIn.getBX();
	  if(_isHOQualityCodes)  qualityCode = 14;
	  
	  //confirming MIP
	  if(_addHitsFromHO && hoMatchSize) {
	    const L1ITMu::TriggerPrimitive & hoIn   = *hoMatch.front();
	    if (  ( hoIn.getBX() == rpcIn.getBX() && hoIn.getBX() == rpcOut.getBX() ) ||
		  (_qualityHORemappingMode >  1 && rpcIn.getBX()==rpcOut.getBX()  && abs(hoIn.getBX()-rpcIn.getBX()) <= 1 ) ) 
	      {
		bool  isMip =  hoIn.isMIP(hoIn.getHOData().Ehit, hoIn.getHOData().Emax, hoIn.getHOData().Emin);
		if (isMip) qualityCode = 15;
	      }
	  }
	}
	
      } else if (rpcInMatchSize){
	const L1ITMu::TriggerPrimitive & rpcIn = *rpcInMatch.front();
	
	if ( dt.getBX() == rpcIn.getBX() || (_qualityRemappingMode > 1 && abs(dt.getBX()-rpcIn.getBX())<=1)) {
	  bx = rpcIn.getBX();
	  if(_isHOQualityCodes) qualityCode = 12;
	  
	  //confirming MIP                                                                                                                                                                               
          if(_addHitsFromHO && hoMatchSize) {
	    const L1ITMu::TriggerPrimitive & hoIn   = *hoMatch.front();
            if (  ( hoIn.getBX() == rpcIn.getBX() ) ||
                  (_qualityHORemappingMode >  1 && abs(hoIn.getBX()-rpcIn.getBX()) <= 1 ) )
	      {
		bool isMip =  hoIn.isMIP(hoIn.getHOData().Ehit, hoIn.getHOData().Emax, hoIn.getHOData().Emin);
		if (isMip) qualityCode = 13;
	      }
	  }
	}
      }
      else if (rpcOutMatchSize){
	const L1ITMu::TriggerPrimitive & rpcOut = *rpcOutMatch.front();

	if ( dt.getBX() == rpcOut.getBX() || (_qualityRemappingMode > 1 && abs(dt.getBX()-rpcOut.getBX())<=1)) {
	  bx = rpcOut.getBX();
	  if(_isHOQualityCodes)  qualityCode = 12;
	  
	  //confirming MIP                                                                                                                                                                               
          if(_addHitsFromHO && hoMatchSize) {
	    const L1ITMu::TriggerPrimitive & hoIn   = *hoMatch.front();
            if (  ( hoIn.getBX() == rpcOut.getBX() ) ||
                  (_qualityHORemappingMode >  1 && abs(hoIn.getBX()-rpcOut.getBX()) <= 1 ) ) {
	      bool isMip =  hoIn.isMIP(hoIn.getHOData().Ehit, hoIn.getHOData().Emax, hoIn.getHOData().Emin);
	      if (isMip) qualityCode = 13;
	    }
          }
	  
	}
      }
      
      
      // add primitive here
      int newBx = dt.getBX();
      if (bx > -999 && dt.getDTData().qualityCode < _useRpcBxForDtBelowQuality){
	newBx=bx;
      }

      if ( ! _is7QualityCodes && !_isHOQualityCodes) { 
	std::cout<<"Does it enter Here"<< std::endl;
	qualityCode = 13;
	switch ( dt.getDTData().qualityCode ) {
	case 4 : qualityCode = 14; break;    // LL // TODO: LL+rpc=13
	case 5 : qualityCode = 15; break;    // HL
	case 6 : qualityCode = 15; break;    // HH 
	default : break;
	}
      }
      
      if(doDebug) std::cout<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ Correlated"<< std::endl;
      if(doDebug) std::cout << "DataFormat:: (bx, Wh, Se, St, RadAng, BenAng, QC, Ts2Code, BxCntCode) --> (" <<newBx <<", "<< wheel <<", " <<sector<<", "<<station <<
		    ","<< dt.getDTData().radialAngle<<", " <<dt.getDTData().bendingAngle<<", "<<  qualityCode <<", "<< dt.getDTData().Ts2TagCode <<", "<< dt.getDTData().BxCntCode <<")"<< std::endl; 
      
      L1MuDTChambPhDigi chamb( newBx, wheel, sector-1, station, dt.getDTData().radialAngle,
			       dt.getDTData().bendingAngle, qualityCode,
			       dt.getDTData().Ts2TagCode, dt.getDTData().BxCntCode );
      phiChambVector.push_back( chamb );
    }
    // END OF BX ANALYSIS FOR CORRELATED TRIGGER

    if(doDebug)     std::cout<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ UnCorrelated with size :: "<< uncorrelated.size() << std::endl;



    

    // BEGIN OF BX ANALYSIS FOR UNCORRELATED TRIGGER
    size_t uncSize = uncorrelated.size();
    for ( size_t idxDt = 0; idxDt < uncSize; ++idxDt ) {
      
      int iDt = uncorrelated.at(idxDt);
      if ( iDt < 0 ) continue;
      const L1ITMu::TriggerPrimitive & dt = *mbltStation.getDtSegments().at(iDt);
      
      /// check if there is a pair of HI+HO at different bx
      int closest = -1;
      int closestIdx = -1;
      double minDeltaPhiDt = 9999999999;
      for ( size_t jdxDt = idxDt+1; jdxDt < uncSize; ++jdxDt ) {
	
	int jDt = uncorrelated.at(jdxDt);
	if ( jDt < 0 ) continue;

	const L1ITMu::TriggerPrimitiveRef & dtM = mbltStation.getDtSegments().at(jDt);
	if ( dt.getBX() == dtM->getBX() || dt.getDTData().qualityCode == dtM->getDTData().qualityCode )
	  continue;

	double deltaPhiDt = fabs( reco::deltaPhi( dt.getCMSGlobalPhi(), dtM->getCMSGlobalPhi() ) );
	if ( deltaPhiDt < minDeltaPhiDt ) {
	  closest = jDt;
	  closestIdx = jdxDt;
	  minDeltaPhiDt=deltaPhiDt;
	}
      }

      /// check if the pair shares the closest rpc/HO hit
      L1ITMu::MBLTCollection::bxMatch match = L1ITMu::MBLTCollection::NOMATCH;
      L1ITMu::MBLTCollection::bxMatch HOmatch = L1ITMu::MBLTCollection::NOMATCH;
      
      if ( closest > 0 && minDeltaPhiDt < 0.05 ) {
	// if ( closest > 0 ) {
	match = mbltStation.haveCommonRpc( iDt, closest );
	
	if( _addHitsFromHO && dt.getDTData().station ==1 )   
	  HOmatch = mbltStation.haveCommonHO( iDt, closest );
      }
      
      // this is just a set of output variables for building L1ITMuDTChambPhDigi
      // int qualityCode = dt.getDTData().qualityCode;
      int bx = -2;
      int radialAngle = 0;
      int bendingAngle = 0;
      
      L1ITMu::PrimitiveCombiner combiner( _resol, _muonGeom);
      if(_addHitsFromHO) 
	L1ITMu::PrimitiveCombiner combiner(_resol, _muonGeom, _hoGeom);
      
      /// association HI/HO provided by the tool
      combiner.addDt( dt );
      
      /// there is a pair HI+HO with a shared inner RPC hit
      if ( match != L1ITMu::MBLTCollection::NOMATCH ) {
	if(doDebug) 	std::cout<<"HI+HO with shared innner RPC hit is found"<< std::endl;
	
	uncorrelated[closestIdx] = -1;
	
	/// association HI/HO provided by the tool
	combiner.addDt( *mbltStation.getDtSegments().at(closest) );
	

	/// redefine quality
	/// qualityCode = 4;
	L1ITMu::TriggerPrimitiveList rpcInMatch  = mbltStation.getRpcInAssociatedStubs( iDt );
	L1ITMu::TriggerPrimitiveList rpcOutMatch = mbltStation.getRpcOutAssociatedStubs( iDt );
	L1ITMu::TriggerPrimitiveList hoMatch     = mbltStation.getHOAssociatedStubs( iDt );

	// add hadron-outer hits
	if(_addHitsFromHO)  {
	  
	  if(HOmatch == L1ITMu::MBLTCollection::INMATCH && dt.getDTData().station == 1)   {
	    const L1ITMu::TriggerPrimitive & hoIn  = *hoMatch.front();	     
	    // std::cout<<"hit, max, min " << ( hoIn.getHOData().Ehit ) <<", "<<  ( hoIn.getHOData().Emax ) <<", " << ( hoIn.getHOData().Emin )<< std::endl;
	    if(  hoIn.isMIP(hoIn.getHOData().Ehit, hoIn.getHOData().Emax, hoIn.getHOData().Emin) )  
	      combiner.addHO( hoIn );
	    // for ( size_t i = 0; i < hoMatch.size(); ++i ) {
	    // const L1ITMu::TriggerPrimitiveRef & dt = hoMatch.at(i); 
	    // std::cout<<"ho_hit : "<< (dt->getHOData().Ehit) << std::endl;
	  }
	}

	
	/// there is a pair HI+HO with a shared inner RPC hit
	if ( match == L1ITMu::MBLTCollection::INMATCH ) {
	  
	  const L1ITMu::TriggerPrimitive & rpcIn = *rpcInMatch.front();
	  combiner.addRpcIn( rpcIn );
	  bx = rpcIn.getBX();
	  
	  /// there is a pair HI+HO with a shared outer RPC hit
	} else if ( match == L1ITMu::MBLTCollection::OUTMATCH ) {
	  
	  const L1ITMu::TriggerPrimitive & rpcOut = *rpcOutMatch.front();
	  combiner.addRpcOut( rpcOut );
	  bx = rpcOut.getBX();
	  
	  /// there is a pair HI+HO with both shared inner and outer RPC hit
	} else if ( match == L1ITMu::MBLTCollection::FULLMATCH ) {
	  
	  const L1ITMu::TriggerPrimitive & rpcIn = *rpcInMatch.front();
	  const L1ITMu::TriggerPrimitive & rpcOut = *rpcOutMatch.front();
	  combiner.addRpcIn( rpcIn );
	  combiner.addRpcOut( rpcOut );
	  bx = rpcIn.getBX();
	}
	
	
      } else { /// there is no match
	
	L1ITMu::TriggerPrimitiveList rpcInMatch    = mbltStation.getRpcInAssociatedStubs( iDt );
	L1ITMu::TriggerPrimitiveList rpcOutMatch   = mbltStation.getRpcOutAssociatedStubs( iDt );
	L1ITMu::TriggerPrimitiveList 	  hoMatch  = mbltStation.getHOAssociatedStubs( iDt);     
	
	size_t rpcInMatchSize  = rpcInMatch.size();
	size_t rpcOutMatchSize = rpcOutMatch.size();
	size_t hoMatchSize     = hoMatch.size();
	
	// add HO, depending if Dt_matched or unmatched
	if ( _addHitsFromHO ) { 
	  if(HOmatch == L1ITMu::MBLTCollection::INMATCH && dt.getDTData().station == 1)   {
	    const L1ITMu::TriggerPrimitive & hoIn  = *hoMatch.front();	     
	    if(  hoIn.isMIP(hoIn.getHOData().Ehit, hoIn.getHOData().Emax, hoIn.getHOData().Emin) )  
	      combiner.addHO( hoIn );
	  } else if (hoMatchSize) {
	    const L1ITMu::TriggerPrimitive &  hoIn = *hoMatch.front();
	    if ( (dt.getBX() == hoIn.getBX() ) || ( _qualityRemappingMode>1 && abs(dt.getBX() - hoIn.getBX()) <=1))  {
	      if(  hoIn.isMIP(hoIn.getHOData().Ehit, hoIn.getHOData().Emax, hoIn.getHOData().Emin) )  
		combiner.addHO( hoIn );
	    }
	  }
	}
	
	/// the uncorrelated has possibly inner and outer confirmation
	if ( rpcInMatchSize && rpcOutMatchSize ) {
	  const L1ITMu::TriggerPrimitive & rpcIn = *rpcInMatch.front();
	  const L1ITMu::TriggerPrimitive & rpcOut = *rpcOutMatch.front();
	  
	  /// only the first is real...
	  // LG try also to reassign BX to single H using RPC BX, e.g. do not ask for DT and RPC to have the same BX
	  if (( dt.getBX() == rpcIn.getBX() && dt.getBX() == rpcOut.getBX() )
	      || (_qualityRemappingMode > 1 && rpcIn.getBX()==rpcOut.getBX() && abs(dt.getBX()-rpcIn.getBX()) <= 1)) {
	    bx = rpcIn.getBX();
	    combiner.addRpcIn( rpcIn );
	    combiner.addRpcOut( rpcOut );
	  } else if ( dt.getBX() == rpcIn.getBX() ) {
	    bx = rpcIn.getBX();
	    combiner.addRpcIn( rpcIn );
	  } else if ( dt.getBX() == rpcOut.getBX() ) {
	    bx = rpcOut.getBX();
	    combiner.addRpcOut( rpcOut );
	  }
	  
	  /// the uncorrelated has a possible inner confirmation
	} else if ( rpcInMatchSize ) {
	  const L1ITMu::TriggerPrimitive & rpcIn = *rpcInMatch.front();
	  if ( dt.getBX() == rpcIn.getBX() || (_qualityRemappingMode>1 && abs(dt.getBX()-rpcIn.getBX()) <= 1 )) {
	    bx = rpcIn.getBX();
	    combiner.addRpcIn( rpcIn );
	  }
	  
	  /// the uncorrelated has a possible outer confirmation
	} else if ( rpcOutMatchSize ) {
	  const L1ITMu::TriggerPrimitive & rpcOut = *rpcOutMatch.front();
	  if ( dt.getBX() == rpcOut.getBX()|| (_qualityRemappingMode>1  && abs(dt.getBX()-rpcOut.getBX()) <= 1)) {
	    bx = rpcOut.getBX();
	    combiner.addRpcOut( rpcOut );
	  }
	  
	}
      }
      
      
      // match found, PrimitiveCombiner has the needed variables already calculated
      if ( combiner.isValid() ) {
	if(doDebug)       std::cout<<"=== I am making a combination ==="<<std::endl;
	combiner.combine();
	radialAngle = combiner.radialAngle();
	bendingAngle = (combiner.bendingAngle() < -511 || combiner.bendingAngle() > 511) ? dt.getDTData().bendingAngle : combiner.bendingAngle( );
      } else {
	// no match found, keep the primitive as it is
	bx = dt.getBX();
	radialAngle = dt.getDTData().radialAngle;
	bendingAngle = dt.getDTData().bendingAngle;
	//if (_qualityRemappingMode==0) 
	// qualityCode = ( qualityCode == 2 ) ? 0 : 1;
      }
      
      int qualityCode = -1;
      if( _isHOQualityCodes ) qualityCode =  combiner.getUncorrelatedQuality16_withHO();
      else  qualityCode = ( _is7QualityCodes ?
			    combiner.getUncorrelatedQuality7() :
			    combiner.getUncorrelatedQuality16() );
      
      if(doDebug)     std::cout<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ Uncorrelated"<< std::endl;       
      if(doDebug)     std::cout << "DataFormat:: (bx, Wh, Se, St, RadAng, BenAng, QC, Ts2Code, BxCntCode) --> (" <<bx <<", "<< wheel <<", " <<sector<<", "<<station <<","<< radialAngle<<
			", " <<bendingAngle<<", "<<  qualityCode <<", "<< dt.getDTData().Ts2TagCode <<", "<< dt.getDTData().BxCntCode <<")"<< std::endl; 
      //     std::cout << "[n]" << qualityCode << std::endl; /// GC
      L1MuDTChambPhDigi chamb( bx, wheel, sector-1, station, radialAngle,
			       bendingAngle, qualityCode,
			       dt.getDTData().Ts2TagCode, dt.getDTData().BxCntCode );
      phiChambVector.push_back( chamb );
      if (abs(bendingAngle)>511||1==1){
	//	std::cout<<"@@@@@@@@@@@@@@@@@@@@@@@@ Uncorrelated Data Format"<< std::endl;
	if(doDebug)       std::cout<<"Got bending angle: "<<bendingAngle<<std::endl;
	if(doDebug)       std::cout<<"Original DT primitive had bending angle: "<<dt.getDTData().bendingAngle<<std::endl;
	if(doDebug)       std::cout<<"Original radial angle: "<<radialAngle<<std::endl;
	if(doDebug)       std::cout<<"Quality: "<<qualityCode<<std::endl;
	if(doDebug)       std::cout<<"Station: "<<station<<std::endl;
      }
      
    } //// end of the Uncorrelated loop

    
    if(_addHitsFromHO ) {
      
      const std::map< std::pair< L1ITMu::TriggerPrimitiveList, L1ITMu::TriggerPrimitiveList >, L1ITMu::TriggerPrimitiveList > rpcHOPairList =  
	mbltStation.getUnassociatedHORpcClusterss( _minDistForHORpcClusterMatch );
      
      std::map< std::pair< L1ITMu::TriggerPrimitiveList, L1ITMu::TriggerPrimitiveList >, L1ITMu::TriggerPrimitiveList >::const_iterator rpcHOPair ;
      
      for( rpcHOPair = rpcHOPairList.begin(); rpcHOPair != rpcHOPairList.end(); ++ rpcHOPair) {
	
	const std::pair <L1ITMu::TriggerPrimitiveList, L1ITMu::TriggerPrimitiveList > & rpcPair = rpcHOPair->first;
	const L1ITMu::TriggerPrimitiveList & inHO = rpcHOPair->second;
	const L1ITMu::TriggerPrimitiveList & inRpc = rpcPair.first;
	const L1ITMu::TriggerPrimitiveList & outRpc = rpcPair.second;
	
      	if( inRpc.empty() && outRpc.empty() && inHO.empty() ) continue;
	
	L1ITMu::PrimitiveCombiner combiner( _resol, _muonGeom, _hoGeom );		
	size_t inSize = inRpc.size();
      	size_t outSize = outRpc.size(); 
      	size_t hoSize = inHO.size(); 	    
      	int station = -1;
      	int sector = -1;
      	int wheel = -5;
	
	if(doDebug) std::cout<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ In&Out&HO Algorithm with (In, Out, ho) :: ("<<inSize <<", "<< outSize <<", "<< hoSize <<")"<< std::endl;
	
	if(inSize && outSize && hoSize) {
	  
	  size_t inPos = 99;
	  size_t outPos = 99;
	  size_t hoPos = 99;
	  
	  for ( size_t i = 0; i < inSize; ++i ) { 
	    const GlobalPoint GPrpcIn =  inRpc.at( i )->getCMSGlobalPoint();
	    
	    for ( size_t j = 0; j < outSize; ++j ) {             
	      const GlobalPoint GPrpcOut   =  outRpc.at( j )->getCMSGlobalPoint();  
	      double FromRpcInToRpcOutdist = sqrt( pow(GPrpcIn.eta()-GPrpcOut.eta(),2)+pow(GPrpcIn.phi()-GPrpcOut.phi(),2) );
	      
	      for( size_t k = 0; k < hoSize; k++ ) { 
		const GlobalPoint GPhoIn    = inHO.at( k )->getCMSGlobalPoint();
		double FromRpcInToHoIndist  =  sqrt( pow(GPrpcOut.eta()-GPhoIn.eta(),2)+pow(GPrpcOut.phi()-GPhoIn.phi(),2) );
		double RpcToHOdist          = FromRpcInToRpcOutdist + FromRpcInToHoIndist ;        
		
		if( RpcToHOdist < _minDistForRpcHOMatch) {
		  inPos = i;
		  outPos = j;
		  hoPos = k;
		  // std::cout<<" distances : "<<  RpcToHOdist <<", "<< _minDistForRpcHOMatch << std::endl;
		  _minDistForRpcHOMatch = RpcToHOdist;
		}
	      }
	    }
	  }
	  
	  if(doDebug)  std::cout<<"HO+RPC case, (i, j, k) : \t"<<inPos<<"\t"<<outPos<<"\t"<<hoPos<< std::endl;

	  if(inPos > 50 && outPos >  50 && hoPos > 50) continue; //quick-fix
	  
	  L1ITMu::TriggerPrimitive rpc1 = (*inRpc.at(inPos));
	  L1ITMu::TriggerPrimitive rpc2 = (*outRpc.at(outPos));
	  L1ITMu::TriggerPrimitive rpc3 = (*inHO.at(hoPos));

	  //shift it above
	  bool IsMipMuon = false;
	  IsMipMuon =  ( (rpc3.getHOData().Ehit > rpc3.getHOData().Emin) && (rpc3.getHOData().Ehit < rpc3.getHOData().Emax) );
	  
	  if(doDebug) 	  std::cout<<"(Emax, Emin, Ehit) :\t"<<(rpc3.getHOData().Emax) <<", \t"<< (rpc3.getHOData().Emin)<<", \t"<<(rpc3.getHOData().Ehit) << std::endl;
	  if (!IsMipMuon) continue;
	  if(doDebug) 	  std::cout<<"MIP found with values (Emax, Emin, Ehit) :\t"<<(rpc3.getHOData().Emax) <<", \t"<< (rpc3.getHOData().Emin)<<", \t"<<(rpc3.getHOData().Ehit) << std::endl;
	  
	  station = rpc1.detId<RPCDetId>().station();
	  sector  = rpc1.detId<RPCDetId>().sector();
	  wheel = rpc1.detId<RPCDetId>().ring();
	  combiner.addRpcIn( rpc1 );
	  combiner.addRpcOut( rpc2 );
	  combiner.addHO( rpc3 );
	}
      
	int bx = -2;
	int radialAngle = 0;
	int bendingAngle = 0;
	double BxCntCode = -1;
	double Ts2TagCode = -1;
	
	if ( combiner.isValid() ) {
	  //	  std::cout<<"=== I am making a combination ==="<<std::endl;
	  combiner.combine();
	  radialAngle = combiner.radialAngle();
	  bx  = combiner.bx(); 
	  Ts2TagCode = 0;                                       
	  BxCntCode = 0;                                   
	  bendingAngle = combiner.bendingAngle();
	}
	

	int qualityCode = -1;
	if( _isHOQualityCodes ) qualityCode =  combiner.getUncorrelatedQuality16_withHO();
	else  qualityCode = ( _is7QualityCodes ?
			      combiner.getUncorrelatedQuality7() :
			      combiner.getUncorrelatedQuality16() );
	if(doDebug)	std::cout<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ HO+RPC Algorithm"<< std::endl;       
	if(doDebug)	std::cout << "DataFormat:: (bx, Wh, Se, St, RadAng, BenAng, QC, Ts2Code, BxCntCode) --> (" <<bx <<", "<< wheel <<", " <<sector<<", "<<station <<","<< radialAngle<<", " <<
			  bendingAngle<<", "<<  qualityCode <<", "<< Ts2TagCode <<", "<< BxCntCode <<")"<< std::endl; 

	L1MuDTChambPhDigi chamb( bx, wheel, sector-1, station, radialAngle,
				 bendingAngle, qualityCode,
				 Ts2TagCode, BxCntCode );
	phiChambVector.push_back( chamb );
      }
      
    } else {
      const std::vector< std::pair< L1ITMu::TriggerPrimitiveList, L1ITMu::TriggerPrimitiveList > >
      rpcPairList = mbltStation.getUnassociatedRpcClusters( 0.05 );

      auto rpcPair = rpcPairList.cbegin();
      auto rpcPairEnd = rpcPairList.cend();
      
      for ( ; rpcPair != rpcPairEnd; ++ rpcPair ) {
	const L1ITMu::TriggerPrimitiveList & inRpc = rpcPair->first;
	const L1ITMu::TriggerPrimitiveList & outRpc = rpcPair->second;
	
	if ( inRpc.empty() && outRpc.empty() ) continue;
	
	L1ITMu::PrimitiveCombiner combiner( _resol, _muonGeom );
	size_t inSize = inRpc.size();
	size_t outSize = outRpc.size();
	int station = -1;
	int sector  = -1;
	int wheel = -5;
	// double qualityCode = 0;
	
	if(doDebug)       std::cout<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ In&Out Algorithm with (In, Out) :: ("<<inSize <<", "<< outSize <<")"<< std::endl;      
	if ( inSize ) {
	if(doDebug) 	std::cout<<"Producer working on IN&&!OUT"<<std::endl;
        size_t inPos = 0;
        double avPhiSin = 0;
        double avPhiCos = 0;
        for ( size_t i = 0; i < inSize; ++i ) {
          double locPhi = inRpc.at(i)->getCMSGlobalPhi();
          // avPhiIn += ( locPhi > 0 ? locPhi : 2*M_PI + locPhi );
	  avPhiSin += sin( locPhi );
	  avPhiCos += cos( locPhi ); 
        }
	avPhiSin /= inSize;
	avPhiCos /= inSize;
	double avPhiIn = atan2( avPhiSin, avPhiCos );
	
        double minDist = fabs( inRpc.at(0)->getCMSGlobalPhi() - avPhiIn );
        for ( size_t i = 1; i < inSize; ++i ) {
          double dist = fabs( inRpc.at(i)->getCMSGlobalPhi() - avPhiIn );
          if ( dist < minDist ) {
            inPos = i;
            minDist = dist;
          }
        }
	
	// const L1ITMu::TriggerPrimitive & rpc = (*inRpc.at(inPos));
	L1ITMu::TriggerPrimitive rpc = (*inRpc.at(inPos));
	rpc.setCMSGlobalPhi( avPhiIn );
        station = rpc.detId<RPCDetId>().station();
        sector  = rpc.detId<RPCDetId>().sector();
        wheel = rpc.detId<RPCDetId>().ring();
        combiner.addRpcIn( rpc );
      }

      if ( outSize ) {
	if(doDebug) 	std::cout<<"Producer working on OUT&&!IN"<<std::endl;
        size_t outPos = 0;
        double avPhiSin = 0;
        double avPhiCos = 0;
        for ( size_t i = 0; i < outSize; ++i ) {
          double locPhi = outRpc.at(i)->getCMSGlobalPhi();
          // avPhiOut += ( locPhi > 0 ? locPhi : 2*M_PI + locPhi );
	  avPhiSin += sin( locPhi );
	  avPhiCos += cos( locPhi );
        }

        //avPhiOut /= outSize;
	avPhiSin /= outSize;
	avPhiCos /= outSize;
	double avPhiOut = atan2( avPhiSin, avPhiCos );
        double minDist = fabs( outRpc.at(0)->getCMSGlobalPhi() - avPhiOut );
        for ( size_t i = 1; i < outSize; ++i ) {
          double dist = fabs( outRpc.at(i)->getCMSGlobalPhi() - avPhiOut );
          if ( dist < minDist ) {
            outPos = i;
            minDist = dist;
          }
	}
        // const L1ITMu::TriggerPrimitive & rpc = (*outRpc.at(outPos));
	L1ITMu::TriggerPrimitive rpc = (*outRpc.at(outPos));
	rpc.setCMSGlobalPhi( avPhiOut );
        station = rpc.detId<RPCDetId>().station();
        sector  = rpc.detId<RPCDetId>().sector();
        wheel = rpc.detId<RPCDetId>().ring();
        combiner.addRpcOut( rpc );
      }
      //if (inSize && outSize) qualityCode=1;
      combiner.combine();
      double radialAngle = combiner.radialAngle();
      double bendingAngle = combiner.bendingAngle();
      double bx = combiner.bx();
      double Ts2TagCode = 0;
      double BxCntCode = 0;
    
      
      int qualityCode = -1;
      if(_isHOQualityCodes) qualityCode =  combiner.getUncorrelatedQuality16_withHO();
      else  qualityCode = ( _is7QualityCodes ?
			    combiner.getUncorrelatedQuality7() :
			    combiner.getUncorrelatedQuality16() );
      
      if ( qualityCode >= 0 ) {
	//	std::cout << "[r]" << qualityCode << std::endl ; /// GC
	
	if(doDebug) 	std::cout << "DataFormat:: (bx, Wh, Se, St, RadAng, BenAng, QC, Ts2Code, BxCntCode) --> (" <<
			  bx <<", "<< wheel <<", " <<sector<<", "<< station <<","<< radialAngle<<", " <<bendingAngle<<", "<<  qualityCode <<", "<< Ts2TagCode <<", "<< BxCntCode <<")"<< std::endl; 

	L1MuDTChambPhDigi chamb( bx, wheel, sector-1, station, radialAngle,
				 bendingAngle, qualityCode,
				 Ts2TagCode, BxCntCode );
	phiChambVector.push_back( chamb );
      }
    }
    }
  }  

  out->setContainer( phiChambVector );
  /// fill event
  iEvent.put(out);

}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1ITMuonBarrelPrimitiveProducer);


// virtual bool nearElement(const GlobalPoint& point, const DetId& id, const double distance)
// {
//   GlobalPoint center = getPosition(id);
//   return sqrt(pow(point.eta()-center.eta(),2)+pow(point.phi()-center.phi(),2)) < distance;
// };
    
