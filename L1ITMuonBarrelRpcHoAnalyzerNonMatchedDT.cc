// -*- C++ -*-
//
// Package:    L1Trigger/L1IntegratedMuonTrigger
// Class:      L1ITMuonBarrelRpcHoAnalyzerNonMatchedDT
// 
/**\class L1ITMuonBarrelRpcHoAnalyzerNonMatchedDT L1ITMuonBarrelRpcHoAnalyzerNonMatchedDT.cc L1Trigger/L1IntegratedMuonTrigger/plugins/L1ITMuonBarrelRpcHoAnalyzerNonMatchedDT.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Pooja Saxena
//         Created:  Mon, 27 Apr 2015 07:12:16 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

// L1IT include files
#include "L1Trigger/L1IntegratedMuonTrigger/interface/MBLTCollection.h"
#include "L1Trigger/L1IntegratedMuonTrigger/interface/MBLTCollectionFwd.h"

// CMSSW include files
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/Records/interface/HcalGeometryRecord.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/MuonDetId/interface/RPCDetId.h"

#include "L1Trigger/L1IntegratedMuonTrigger/interface/PrimitiveCombiner.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TTree.h"

#define maxcount 1000 

//
// class declaration
//

class L1ITMuonBarrelRpcHoAnalyzerNonMatchedDT : public edm::EDAnalyzer {

public:

  //constructor
  explicit L1ITMuonBarrelRpcHoAnalyzerNonMatchedDT(const edm::ParameterSet&);

  //destructor
  ~L1ITMuonBarrelRpcHoAnalyzerNonMatchedDT();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  edm::InputTag _mbltCollectionInput;
  const L1ITMu::PrimitiveCombiner::resolutions _resol;
  const int _qualityRemappingMode;
  const int _qualityHORemappingMode;
  const int _useRpcBxForDtBelowQuality;
  double _minDistForRpcHOMatch;
  const double _minDistForHORpcClusterMatch;
  edm::ESHandle<DTGeometry> _muonGeom;
  edm::ESHandle<CaloGeometry> _hoGeom;

  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  // ----------member data ---------------------------
public:
  UInt_t event;
  UInt_t run;

  Int_t wodt_wheel;
  UInt_t wodt_sector;
  UInt_t wodt_station;
  Float_t wodt_radialAngle;
  Float_t wodt_bendingAngle;
  UInt_t wodt_MIP;
  Float_t wodt_Ehit;

  //  bool_t  dtBXmatched[count_dt];
  //  bool_t  hoBXmatched[count_dt];
  UInt_t dtBX;
  Float_t dtGEta;
  Float_t dtGPhi;
  Float_t rpcInGEta;
  Float_t rpcInGPhi;
  Float_t rpcOutGEta;
  Float_t rpcOutGPhi;
  Float_t HOGEta;
  Float_t HOGPhi;
  UInt_t rpcInBX;
  UInt_t rpcOutBX;
  UInt_t hoBX;
  UInt_t countHO;
  Int_t dtmatch_wheel;
  UInt_t dtmatch_sector;
  UInt_t dtmatch_station;
  Float_t dtmatch_radialAngle;
  Float_t dtmatch_bendingAngle;
  UInt_t dtmatch_qualityCode;
 
  
  // book histos
  void bookHistos();
  TTree* tree;
  TH1D*  nEvents;
  bool  doDebug = false;
  bool  tkEvtA = false;
  bool  tkEvtB = false;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
L1ITMuonBarrelRpcHoAnalyzerNonMatchedDT::L1ITMuonBarrelRpcHoAnalyzerNonMatchedDT(const edm::ParameterSet& iConfig)
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
    //    _minDistForRpcHOMatch( iConfig.getParameter<double>("minDistForRpcHOMatch") ),
    _minDistForHORpcClusterMatch(iConfig.getParameter<double>("minDistForHORpcClusterMatch") )
{
  //now do what ever initialization is needed
}


L1ITMuonBarrelRpcHoAnalyzerNonMatchedDT::~L1ITMuonBarrelRpcHoAnalyzerNonMatchedDT()
{
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called for each event  ------------
void L1ITMuonBarrelRpcHoAnalyzerNonMatchedDT::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  if(doDebug)  std::cout<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ event "<< iEvent.id().event() << std::endl;
  using namespace edm;
  
  iSetup.get<MuonGeometryRecord>().get(_muonGeom);
  iSetup.get<CaloGeometryRecord>().get(_hoGeom);
  event      = iEvent.id().event();
  run        = iEvent.id().run();
  
  edm::Handle<L1ITMu::MBLTContainer> mbltContainer;
  iEvent.getByLabel( _mbltCollectionInput, mbltContainer );
  
  L1ITMu::MBLTContainer::const_iterator st    = mbltContainer->begin();
  L1ITMu::MBLTContainer::const_iterator stend = mbltContainer->end();
  
  if(doDebug) std::cout<<"Entering into mblt loop"<< std::endl;  

  //initialization
  tkEvtA      = false;
  tkEvtB      = false;
  
  ////
  for ( ; st != stend; ++st ) {
    
    const L1ITMu::MBLTCollection & mbltStation = st->second;
    
    /// useful index
    int station = mbltStation.station();
    int wheel   = mbltStation.wheel();
    int sector  = mbltStation.sector();
    
    /// get dt to rpc associations
    size_t dtListSize = mbltStation.getDtSegments().size();
    std::vector<size_t> correlated;
    if(doDebug && dtListSize > 0 )    std::cout<<"mbltCollection  is :"<< dtListSize << std::endl;

    /// loop over DT collections
    for ( size_t iDt = 0; iDt < dtListSize; ++iDt ) {
      
      // only station == 1
      const L1ITMu::TriggerPrimitiveRef & dt = mbltStation.getDtSegments().at(iDt);
      int station = dt->getDTData().station;
      int dtQuality = dt->getDTData().qualityCode;
      
      if (station != 1) continue;
      if(dtQuality == 4 || dtQuality == 5 || dtQuality == 6)  correlated.push_back( iDt );
    }// for ( size_t iDt = 0; iDt < dtListSize; ++iDt )
    
    
    // do the analysis only for the correlated DT, so the results can be compare later wit RPC+HO case
    size_t cSize = correlated.size();
    
    
    for( size_t idxDt = 0; idxDt < cSize; ++idxDt )	{
      if(doDebug)    std::cout<<"Inside, idxDt: "<<idxDt <<" with correletaed Size: "<< cSize<<  std::endl; 
      
      int iDt = correlated.at(idxDt);
      if ( iDt < 0 ) continue;
      
      if(doDebug) std::cout<<"@@@ DT-algo::correlated DTs found at "<< idxDt << std::endl;
      
      const L1ITMu::TriggerPrimitive & dt = *mbltStation.getDtSegments().at(iDt);
      
      // get all collections
      L1ITMu::TriggerPrimitiveList rpcInMatch  = mbltStation.getRpcInAssociatedStubs( iDt );
      L1ITMu::TriggerPrimitiveList rpcOutMatch = mbltStation.getRpcOutAssociatedStubs( iDt );
      L1ITMu::TriggerPrimitiveList hoMatch     = mbltStation.getHOInner(); //mbltStation.getHOAssociatedStubs( iDt);
      
      size_t rpcInMatchSize  = rpcInMatch.size();
      size_t rpcOutMatchSize = rpcOutMatch.size();
      size_t hoMatchSize     = hoMatch.size();
      
      if(doDebug) std::cout<<"DT-algo::sizes are, "<<rpcInMatchSize<<", "<< rpcOutMatchSize<<", "<< hoMatchSize << std::endl;

      // collection should not be empty
      if (!( rpcInMatchSize && rpcOutMatchSize && hoMatchSize)) continue;
      //if (!( rpcInMatchSize && rpcOutMatchSize)) continue;
      
      const  L1ITMu::TriggerPrimitive & rpcIn  = *rpcInMatch.front();
      const  L1ITMu::TriggerPrimitive & rpcOut = *rpcOutMatch.front();
      
      bool dtBXmatched  =  ( ( dt.getBX() == rpcIn.getBX() && dt.getBX() == rpcOut.getBX() )
			     || (_qualityRemappingMode > 1 && rpcIn.getBX()==rpcOut.getBX() && abs(dt.getBX()-rpcIn.getBX())<=1) ) ;
      
      if(doDebug) std::cout<<"matching the BX"<< std::endl;
      if(!dtBXmatched) continue;
      
      tkEvtA = true;
      // lets store the info with DT matching
      dtBX                  = dt.getBX();
      dtGEta                = dt.getCMSGlobalEta();
      dtGPhi                = dt.getCMSGlobalPhi();
      dtmatch_sector        = sector-1;  //dt.getDTData().sector;
      dtmatch_wheel         = wheel;     //dt.getDTData().wheel;
      dtmatch_station       = station;   //dt.getDTData().station;
      dtmatch_radialAngle   = dt.getDTData().radialAngle;
      dtmatch_bendingAngle  = dt.getDTData().bendingAngle;
      dtmatch_qualityCode   = dt.getDTData().qualityCode;
      //      dtmatch_MIP           = hoIn.isMIP(hoIn.getHOData().Ehit, hoIn.getHOData().Emax, hoIn.getHOData().Emin);
      
      if(doDebug) std::cout <<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ DT Only"<< std::endl;
      if(doDebug) std::cout<<"(Wh, Se, St, benA, radA, mip): ("<<dtmatch_wheel <<", "<<dtmatch_sector<<", "
		    <<dtmatch_station <<", "<<dtmatch_bendingAngle <<", "<<dtmatch_radialAngle <<")"<< std::endl;
      if(doDebug) std::cout<<"Quality : "<<dt.getDTData().qualityCode<< std::endl;
      
      if(tkEvtA) break;
    }//for( size_t idxDt = 0; idxDt < cSize; ++idxDt )

    if(!tkEvtA) continue;

  
    //////// RPC+HO algorithem
    
    // get all collections
    const  L1ITMu::TriggerPrimitiveList & rpcInMatch  = mbltStation.getRpcInner();
    const  L1ITMu::TriggerPrimitiveList & rpcOutMatch = mbltStation.getRpcOuter();
    const  L1ITMu::TriggerPrimitiveList & hoMatch     = mbltStation.getHOInner();
    
    size_t rpcInMatchSize  = rpcInMatch.size();
    size_t rpcOutMatchSize = rpcOutMatch.size();
    size_t hoMatchSize     = hoMatch.size();
  
  
    // collection should not be empty
    if (!( rpcInMatchSize && rpcOutMatchSize && hoMatchSize)) continue;
    
    if(doDebug) std::cout <<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ RPC+HO"<< std::endl;
    if(doDebug)    std::cout<<"Size of the collection: "<<rpcInMatchSize<<", "<< rpcOutMatchSize<<", "<< hoMatchSize<< std::endl; 
    
    const  L1ITMu::TriggerPrimitive & rpcIn  = *rpcInMatch.front();
    const  L1ITMu::TriggerPrimitive & rpcOut = *rpcOutMatch.front();
    const  L1ITMu::TriggerPrimitive & hoIn   = *hoMatch.front();
    
    
    bool hoBXmatched  =  (  ( hoIn.getBX() == rpcIn.getBX() && hoIn.getBX() == rpcOut.getBX() ) ||
			    (_qualityHORemappingMode >  1 && rpcIn.getBX()==rpcOut.getBX()  && abs(hoIn.getBX()-rpcIn.getBX()) <= 1 ) );
    
    if(doDebug) std::cout<<"matching BX for Rpc&HO"<< std::endl;
    if(!hoBXmatched) continue;
    
    
    // collection should be provided here and then the checks should be done
    _minDistForRpcHOMatch = 0.5; 


    const std::map< std::pair< L1ITMu::TriggerPrimitiveList, L1ITMu::TriggerPrimitiveList >, L1ITMu::TriggerPrimitiveList > rpcHOPairList =  
      mbltStation.getAllHORpcClusterss(_minDistForHORpcClusterMatch );
    if(rpcHOPairList.empty()) continue;

    std::map< std::pair< L1ITMu::TriggerPrimitiveList, L1ITMu::TriggerPrimitiveList >, L1ITMu::TriggerPrimitiveList >::const_iterator rpcHOPair = rpcHOPairList.begin();

    const std::pair <L1ITMu::TriggerPrimitiveList, L1ITMu::TriggerPrimitiveList > & rpcPair = rpcHOPair->first;
    const L1ITMu::TriggerPrimitiveList & inHO = rpcHOPair->second;
    const L1ITMu::TriggerPrimitiveList & inRpc = rpcPair.first;
    const L1ITMu::TriggerPrimitiveList & outRpc = rpcPair.second;
    
    if( inRpc.empty() && outRpc.empty() && inHO.empty() ) continue;

    
    L1ITMu::PrimitiveCombiner combiner( _resol, _muonGeom, _hoGeom );		
    size_t inSize  = inRpc.size();
    size_t outSize = outRpc.size(); 
    size_t hoSize  = inHO.size(); 	    
    Float_t radialAngle = 0;
    Float_t bendingAngle = 0;
    
    if(doDebug) std::cout<<"Containor has hits, Num-hits (In, Out, ho):("<<inSize <<", "<< outSize <<", "<< hoSize <<
		  ") for threshold minDistRpcHODist: "<<_minDistForRpcHOMatch << std::endl;
    
    if(inSize && outSize && hoSize) {
      size_t inPos  = 999;
      size_t outPos = 999;
      size_t hoPos  = 999;

      for ( size_t i = 0; i < inSize; ++i ) { 
	double phi_RpcIn = inRpc.at( i )->getCMSGlobalPhi();
       	if(doDebug) std::cout<<"rpcIn id: \t" <<i <<" \t phi-I \t"<< phi_RpcIn << std::endl;
	
	for ( size_t j = 0; j < outSize; ++j ) {             
	  double phi_RpcOut = outRpc.at( j )->getCMSGlobalPhi();
	  double deltaPhi_RpcInRpcOut = fabs( reco::deltaPhi( phi_RpcIn, phi_RpcOut ) );
	  if(doDebug) std::cout<<"rpcOut id: \t" << j <<"\t  phi-II \t"<< phi_RpcOut << " dPhi \t"<< deltaPhi_RpcInRpcOut<< std::endl;
	  
	  for( size_t k = 0; k < hoSize; k++ ) { 
	    double phi_ho = inHO.at( k )->getCMSGlobalPhi(); 
	    double deltaPhi_RpcOutHO = fabs( reco::deltaPhi( phi_RpcIn, phi_ho));
	    
	    double  RpcToHOdist          = deltaPhi_RpcInRpcOut + deltaPhi_RpcOutHO;
	    if(doDebug) std::cout<<" ho id : \t"<<k << "\t phi-III \t "<<phi_ho <<"  phi_HORpcOut \t " << deltaPhi_RpcOutHO << " dPhi_rpcHO : \t"<< RpcToHOdist << std::endl;
	    
	    if( RpcToHOdist < _minDistForRpcHOMatch) {
	      inPos = i;
	      outPos = j;
	      hoPos = k;
	      if(doDebug) std::cout<<" RPC_to_HO distances : \t"<<  RpcToHOdist <<" threshold_dist, \t"<< _minDistForRpcHOMatch << std::endl;
	      _minDistForRpcHOMatch = RpcToHOdist;

	    }
	  }
	}
      } // for ( size_t i = 0; i < inSize; ++i )
      
      if(doDebug)  std::cout<<"Rpc-HO algo:: (i, j, k) : \t"<<inPos<<"\t"<<outPos<<"\t"<<hoPos<<" with minDistance :"<< _minDistForRpcHOMatch << std::endl;
      
      if(inPos > 50 && outPos >  50 && hoPos > 50) {
	if(doDebug) std::cout<<"Hits pointer is too Big, reversing!"<< std::endl;
	continue;
      }
      
      L1ITMu::TriggerPrimitive rpc1 = (*inRpc.at(inPos));
      L1ITMu::TriggerPrimitive rpc2 = (*outRpc.at(outPos));
      L1ITMu::TriggerPrimitive rpc3 = (*inHO.at(hoPos));
      
      
      if(doDebug)  std::cout<<"found the hits, adding to the Combiner"<< std::endl;
      combiner.addRpcIn( rpc1 );
      combiner.addRpcOut( rpc2 );
      combiner.addHO( rpc3 );

      if ( combiner.isValid() ) {
	tkEvtB = true;
	combiner.combine();
	radialAngle = combiner.radialAngle();
	bendingAngle = combiner.bendingAngle( ); 
      } else {
	// no match found, keep the primitive as it is
	radialAngle = -999.0;
	bendingAngle = -999.0;
      }

      countHO               = hoMatchSize;      
      rpcInBX               = rpc1.getBX();
      rpcOutBX              = rpc2.getBX();
      hoBX                  = rpc3.getBX();
      wodt_station       = rpc1.detId<RPCDetId>().station();
      wodt_sector        = rpc1.detId<RPCDetId>().sector()-1;
      wodt_wheel         = rpc1.detId<RPCDetId>().ring();
      wodt_radialAngle   = radialAngle;
      wodt_bendingAngle  = bendingAngle;
      wodt_MIP           = rpc3.isMIP(rpc3.getHOData().Ehit, rpc3.getHOData().Emax, rpc3.getHOData().Emin);
      wodt_Ehit          = rpc3.getHOData().Ehit;
      rpcInGEta          = rpc1.getCMSGlobalEta();
      rpcInGPhi          = rpc1.getCMSGlobalPhi();
      rpcOutGEta         = rpc2.getCMSGlobalEta();
      rpcOutGPhi         = rpc2.getCMSGlobalPhi();
      HOGEta             = rpc3.getCMSGlobalEta();    
      HOGPhi             = rpc3.getCMSGlobalPhi();    

      if(doDebug) std::cout<<"(Wh, Se, St, benA, radA, mip): ("<< wodt_wheel <<", "<< wodt_sector <<", "<< 
		    wodt_station  <<", "<< bendingAngle <<", "<< radialAngle <<", "<< wodt_MIP <<")"<< std::endl;
      if(doDebug) std::cout <<"MIP, Ehit, Emax, Emin :"<< rpc3.getHOData().Ehit <<", "<< rpc3.getHOData().Emax<<", "<< rpc3.getHOData().Emin << std::endl; 
    } // if(inSize && outSize && hoSize)
    
    if(tkEvtA && tkEvtB) break;
  }// for ( ; st != stend; ++st )
  
  if(doDebug) std::cout<<"takeEvent :"<<tkEvtA <<", "<< tkEvtB<< std::endl;
  if(tkEvtA && tkEvtB)    tree->Fill();
  if(tkEvtA && tkEvtB)     nEvents->Fill(1);
  else  nEvents->Fill(0);
}


// ------------ method called once each job just before starting event loop  ------------
void 
  L1ITMuonBarrelRpcHoAnalyzerNonMatchedDT::beginJob() {

  edm::Service<TFileService> FS;
  tree = FS->make<TTree>("l1HOmuonTrigger","l1HOmuonTrigger",1);
  nEvents = FS->make<TH1D>("nEvents","nEvents", 2, -0.5, 1.5);

  tree->Branch("event", &event, "event/i");
  tree->Branch("run", &run, "run/i");

  tree->Branch("dtmatch_wheel",&dtmatch_wheel,"dtmatch_wheel/I");
  tree->Branch("dtmatch_sector",&dtmatch_sector,"dtmatch_sector/i");
  tree->Branch("dtmatch_station", &dtmatch_station,"dtmatch_station/i");
  tree->Branch("dtmatch_radialAngle",&dtmatch_radialAngle,"dtmatch_radialAngle/f");
  tree->Branch("dtmatch_bendingAngle",&dtmatch_bendingAngle ,"dtmatch_bendingAngle/f");
  tree->Branch("dtmatch_qualityCode",&dtmatch_qualityCode,"dtmatch_qualityCode/i");
  //  tree->Branch("dtmatch_MIP",&dtmatch_MIP,"dtmatch_MIP/i");
  tree->Branch("dtBX",&dtBX,"dtBX/i");
  tree->Branch("dtGEta",&dtGEta,"dtGEta/f");
  tree->Branch("dtGPhi",&dtGPhi,"dtGPhi/f");
  tree->Branch("rpcInBX",&rpcInBX,"rpcInBX/i");
  tree->Branch("rpcOutBX",&rpcOutBX,"rpcOutBX/i");
  tree->Branch("hoBX",&hoBX,"hoBX/i");
  tree->Branch("countHO",&countHO,"countHO/I");

  tree->Branch("rpcInGEta",&rpcInGEta,"rpcInGEta/f");
  tree->Branch("rpcInGPhi",&rpcInGPhi,"rpcInGPhi/f");

  tree->Branch("rpcOutGEta",&rpcOutGEta,"rpcOutGEta/f");
  tree->Branch("rpcOutGPhi",&rpcOutGPhi,"rpcOutGPhi/f");

  tree->Branch("HOGEta",&HOGEta,"HOGEta/f");
  tree->Branch("HOGPhi",&HOGPhi,"HOGPhi/f");

  tree->Branch("wodt_wheel", &wodt_wheel, "wodt_wheel/I");
  tree->Branch("wodt_sector",&wodt_sector,"wodt_sector/i");
  tree->Branch("wodt_station",&wodt_station,"wodt_station/i");
  tree->Branch("wodt_radialAngle",&wodt_radialAngle,"wodt_radialAngle/f");
  tree->Branch("wodt_bendingAngle",&wodt_bendingAngle,"wodt_bendingAngle/f");
  tree->Branch("wodt_MIP",&wodt_MIP,"wodt_MIP/i");
  tree->Branch("wodt_Ehit",&wodt_Ehit,"wodt_Ehit/f");  
}
 
// ------------ method called once each job just after ending the event loop  ------------
void 
L1ITMuonBarrelRpcHoAnalyzerNonMatchedDT::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
void 
L1ITMuonBarrelRpcHoAnalyzerNonMatchedDT::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
L1ITMuonBarrelRpcHoAnalyzerNonMatchedDT::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
L1ITMuonBarrelRpcHoAnalyzerNonMatchedDT::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
L1ITMuonBarrelRpcHoAnalyzerNonMatchedDT::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
L1ITMuonBarrelRpcHoAnalyzerNonMatchedDT::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(L1ITMuonBarrelRpcHoAnalyzerNonMatchedDT);
