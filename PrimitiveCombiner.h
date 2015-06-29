#ifndef L1Trigger_L1IntegratedMuonTrigger_PrimitiveCombiner_h_
#define L1Trigger_L1IntegratedMuonTrigger_PrimitiveCombiner_h_
// 
// Class: L1ITMu:: PrimitiveCombiner
//
// Info: This class combine information from DT and/or RPC primitives
//       in order to calculate better phi/phiBending
//
// Author: 
//

#include "FWCore/Framework/interface/ESHandle.h"
#include <cmath>

class DTGeometry;
class CaloGeometry;

namespace L1ITMu {

  class TriggerPrimitive;
  class PrimitiveCombiner {

  public :
    
    /// a struct useful for resulution info sharing
    struct resolutions {
      double xDt;
      double xRpc;
      double phibDtCorr;
      double phibDtUnCorr;
      double xHO;
    resolutions( const double & xDtResol, const double & xRpcResol,
		 const double & phibDtCorrResol, const double & phibDtUnCorrResol ):
      xDt( xDtResol ), xRpc( xRpcResol ),
	phibDtCorr( phibDtCorrResol ), phibDtUnCorr( phibDtUnCorrResol ) {};
      
    resolutions( const double & xDtResol, const double & xRpcResol,
		 const double & phibDtCorrResol, const double & phibDtUnCorrResol, const double & xHOResol ):
      xDt( xDtResol ), xRpc( xRpcResol ),
	phibDtCorr( phibDtCorrResol ), phibDtUnCorr( phibDtUnCorrResol ), xHO(xHOResol) {};
    };

  public :
    explicit PrimitiveCombiner( const resolutions & res, edm::ESHandle<DTGeometry> & muonGeom );
    explicit PrimitiveCombiner( const resolutions & res, edm::ESHandle<DTGeometry> & muonGeom, edm::ESHandle<CaloGeometry> & hoGeom );

    /// feed the combiner with the available primitives
    void addDt( const L1ITMu::TriggerPrimitive & prim );
    void addDtHI( const L1ITMu::TriggerPrimitive & prim );
    void addDtHO( const L1ITMu::TriggerPrimitive & prim );
    void addRpcIn( const L1ITMu::TriggerPrimitive & prim );
    void addRpcOut( const L1ITMu::TriggerPrimitive & prim );
    void addHO( const L1ITMu::TriggerPrimitive & prim );

    /// do combine the primitives
    void combine();

    /// output result variables
    inline int bx() const { return _bx;};
    inline int radialAngle() const { return _radialAngle;};
    inline int bendingAngle() const { return _bendingAngle;};
    inline int bendingResol() const { return _bendingResol;};

    /// valid if we have at least: 1 rpc; 1 dt + 1 any
    bool isValid() const {
      int ret = _dtHI ? 1 : 0;
      ret += _dtHO ? 1 : 0;
      ret += _rpcIn ? 2 : 0;
      ret += _rpcOut ? 2 : 0;
      return ret > 1 ;
    };

    /// FIXME : Calculates new phiBending, check how to use 
    inline float phiBCombined( const float & xDt, const float & zDt,
			       const float & xRpc, const float & zRpc )
    {
      return (xRpc - xDt) / (zRpc - zDt);
    };
    /// FIXME END

    /// FIXME : Calculates new phiBending resolution
    inline float phiBCombinedResol( const float & resol_xDt,
				    const float & resol_xRpc,
				    const float & zDt,
				    const float & zRpc
				    )
    {
      return sqrt( resol_xRpc*resol_xRpc + resol_xDt*resol_xDt )/abs(zDt-zRpc);
    };
    /// FIXME END

    int getUncorrelatedQuality7() const {

      int qualityCode = 0;
      if ( _dtHI && _dtHO ) {
	if ( _rpcIn && _rpcOut ) qualityCode = 5;
	else qualityCode = 5;

      } else if ( _dtHO ) {// HO quality == 3
	if ( _rpcIn && _rpcOut ) qualityCode = 4;
	else if ( _rpcOut ) qualityCode = 4;
	else if ( _rpcIn ) qualityCode = 4;
	else  qualityCode = 2;

      } else if ( _dtHI ) {// HI, quality == 2
	if ( _rpcIn && _rpcOut ) qualityCode = 3;
	else if ( _rpcOut ) qualityCode = 3;
	else if ( _rpcIn ) qualityCode = 3;
	else  qualityCode = 1;
      } else {
	if ( _rpcIn && _rpcOut ) qualityCode = 0;
	else if ( _rpcOut ) qualityCode = -1;
	else if ( _rpcIn ) qualityCode = -1;
      }
      return qualityCode;
    }

    int getUncorrelatedQuality16() const {

      int qualityCode = 0;
      if ( _dtHI && _dtHO ) {
	if ( _rpcIn && _rpcOut ) qualityCode = 12;
	else qualityCode = 11;
      } else if ( _dtHO ) {// HO quality == 3
	if ( _rpcIn && _rpcOut ) qualityCode = 10;
	else if ( _rpcOut ) qualityCode = 8;
	else if ( _rpcIn ) qualityCode = 6;
	else  qualityCode = 4;
      } else if ( _dtHI ) {// HI, quality == 2
	if ( _rpcIn && _rpcOut ) qualityCode = 9;
	else if ( _rpcOut ) qualityCode = 7;
	else if ( _rpcIn ) qualityCode = 5;
	else  qualityCode = 3;
      } else {
	if ( _rpcIn && _rpcOut ) qualityCode = 2;
	else if ( _rpcOut ) qualityCode = 1;
	else if ( _rpcIn ) qualityCode = 0;
      }
      return qualityCode;
    }


    int getUncorrelatedQuality16_withHO() const {

      int qualityCode = 0;
      if ( _dtHI && _dtHO ) {
        if ( _rpcIn && _rpcOut && _HO ) qualityCode = 15;
	else if ( _rpcIn && _HO )            qualityCode = 13;
        else if ( _rpcOut && _HO )           qualityCode = 13;
        else if ( _rpcIn && _rpcOut )        qualityCode = 14;
	else                            qualityCode = 12;

      } else if ( _dtHO ) { // HO quality == 3                                                                                                                                                           
        if ( _rpcIn && _rpcOut && _HO)  qualityCode = 11;
        else if ( _rpcOut &&_rpcIn )    qualityCode = 10;
        else if ( _rpcIn && _HO )       qualityCode = 9;
        else if ( _rpcOut && _HO )      qualityCode = 9;
        else if ( _rpcOut )             qualityCode = 8;
        else if ( _rpcIn  )             qualityCode = 8;
        else                            qualityCode = 3;

      } else if ( _dtHI ) { // HO quality == 2
        if ( _rpcIn && _rpcOut && _HO)  qualityCode = 7;
        else if ( _rpcOut &&_rpcIn )    qualityCode = 6;
        else if ( _rpcIn && _HO )       qualityCode = 5;
        else if ( _rpcOut && _HO )      qualityCode = 5;
        else if ( _rpcOut )             qualityCode = 4;
        else if ( _rpcIn  )             qualityCode = 4;
        else                            qualityCode = 2;

      } else if ( _rpcOut ||_rpcIn ) {
	if (_rpcIn && _rpcOut && _HO )  qualityCode = 1;
	else if (_rpcIn && _rpcOut )    qualityCode = 0;
	else if ( _rpcOut )             qualityCode = -1;
	else if ( _rpcIn )              qualityCode = -1;
	else                            qualityCode = 0;
      }
      return qualityCode;
    }


    int getUncorrelatedQuality32_withHO() const {

      int qualityCode = 0;
      if ( _dtHI && _dtHO ) {
        if ( _rpcIn && _rpcOut && _HO ) qualityCode = 31;
	else if ( _rpcIn && _rpcOut )   qualityCode = 30;
	else if ( _rpcOut && _HO )      qualityCode = 29; // longer arm length, so higher quality
	else if ( _rpcIn && _HO )       qualityCode = 28;
	else if ( _rpcOut )             qualityCode = 27;
	else if ( _rpcIn )              qualityCode = 26;
	else                            qualityCode = 25;

      } else if ( _dtHO ) { // HO quality == 3                                                                                                                                                           
        if ( _rpcIn && _rpcOut && _HO)  qualityCode = 24;
        else if ( _rpcOut &&_rpcIn )    qualityCode = 23;
        else if ( _rpcOut && _HO )      qualityCode = 21;
        else if ( _rpcIn && _HO )       qualityCode = 20;
        else if ( _rpcOut )             qualityCode = 17;
        else if ( _rpcIn  )             qualityCode = 16;
        else                            qualityCode = 12;

      } else if ( _dtHI ) { // HO quality == 2
        if ( _rpcIn && _rpcOut && _HO)  qualityCode = 22;
        else if ( _rpcOut &&_rpcIn )    qualityCode = 19;
        else if ( _rpcIn && _HO )       qualityCode = 18;
        else if ( _rpcOut && _HO )      qualityCode = 15;
        else if ( _rpcOut )             qualityCode = 14;
        else if ( _rpcIn  )             qualityCode = 13;
        else                            qualityCode = 11;

      } else if ( _rpcOut ||_rpcIn ) {
	if (_rpcIn && _rpcOut && _HO )  qualityCode = 10;
	else if (_rpcIn && _rpcOut )    qualityCode = 9;
	else if ( _rpcOut )             qualityCode = 8;
	else if ( _rpcIn )              qualityCode = 8;
	else                            qualityCode = 0;
      }
      return qualityCode;
    }


  private :

    /// a struct for internal usage: store results
    struct results {
      double radialAngle;
      double bendingAngle;
      double bendingResol;
    results() : radialAngle(0), bendingAngle(0), bendingResol(0) {};
    };


    /// Calculates new phiBending, check how to use weights
    results combineDt( const L1ITMu::TriggerPrimitive * dt,
		       const L1ITMu::TriggerPrimitive * rpc );

    results dummyCombineDt( const L1ITMu::TriggerPrimitive * dt);

    /// Calculates new phiBending, check how to use weights
    results combineDtRpc( const L1ITMu::TriggerPrimitive * dt,
			  const L1ITMu::TriggerPrimitive * rpc );

    /// Calculates new phiBending, check how to use weights
    results combineRpcRpc( const L1ITMu::TriggerPrimitive * rpc1,
			   const L1ITMu::TriggerPrimitive * rpc2 );
 
 
    int radialAngleFromGlobalPhi( const L1ITMu::TriggerPrimitive * rpc );

  private :
    resolutions _resol;
    edm::ESHandle<DTGeometry> _muonGeom;
    edm::ESHandle<CaloGeometry> _hoGeom;

    int _bx;
    int _radialAngle;
    int _bendingAngle;
    int _bendingResol;

    const L1ITMu::TriggerPrimitive * _dtHI;
    const L1ITMu::TriggerPrimitive * _dtHO;
    const L1ITMu::TriggerPrimitive * _rpcIn;
    const L1ITMu::TriggerPrimitive * _rpcOut;
    const L1ITMu::TriggerPrimitive * _HO;

  };
}

#endif
