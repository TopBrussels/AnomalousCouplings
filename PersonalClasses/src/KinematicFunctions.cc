#include "../interface/KinematicFunctions.h"

float KinematicFunctions::CosTheta(TLorentzVector *leptTop, TLorentzVector *leptW, TLorentzVector *lepton){

    float CosTheta = 0;
    //-----    Initializing boost variables    -----//
    TLorentzVector leptonWRestFrame, WLeptTopRestFrame, topBooster;
    leptonWRestFrame = *lepton;
    topBooster = *leptTop;
    WLeptTopRestFrame = *leptW;

    //-----    Applying boost on muon and W    -----//
    leptonWRestFrame.Boost(-WLeptTopRestFrame.BoostVector());
    WLeptTopRestFrame.Boost(-topBooster.BoostVector());
    
    //-----   Calculating cos theta:   -----
    CosTheta = ((WLeptTopRestFrame.Vect()).Dot(leptonWRestFrame.Vect()))/(((WLeptTopRestFrame.Vect()).Mag())*((leptonWRestFrame.Vect()).Mag()));

    return CosTheta;
}
