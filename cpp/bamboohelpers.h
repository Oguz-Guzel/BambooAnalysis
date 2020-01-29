#pragma once

#include "Math/LorentzVector.h"
#include "Math/PtEtaPhiE4D.h"

#include "ROOT/RDF/InterfaceUtils.hxx"
#include "ROOT/RDF/HistoModels.hxx"

namespace rdfhelpers {

template<typename N1,typename N2,typename N3>
bool in_range( N1 lower, N2 value, N3 upper )
{
  return ( lower < value ) && ( value < upper );
}

template<typename LorentzVector>
ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float> > withMass( const LorentzVector& p4, float mass )
{
  return ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float> >{p4.Pt(), p4.Eta(), p4.Phi(), mass};
}

// helper methods for min/max (element) reduce expressions, see op.rng_max_element_by
template<typename T1,typename T2,typename T3>
typename std::pair<T1,T2> maxPairBySecond( typename std::pair<T1,T2> a, T1 bFirst, T3 bSecond )
{
  return ( bSecond > a.second ) ? std::pair<T1,T2>(bFirst, bSecond) : a;
}
template<typename T1,typename T2, typename T3>
typename std::pair<T1,T2> minPairBySecond( typename std::pair<T1,T2> a, T1 bFirst, T3 bSecond )
{
  return ( bSecond < a.second ) ? std::pair<T1,T2>(bFirst, bSecond) : a;
}

template<typename RDF,typename VAR1, typename WEIGHT>
ROOT::RDF::RResultPtr<TH1D> Histo1D(RDF& rdf, const ROOT::RDF::TH1DModel& model, std::string_view vName, std::string_view wName)
{
    return rdf.template Histo1D<VAR1,WEIGHT>(model, vName, wName);
}

template<typename RDF,typename VAR1>
ROOT::RDF::RResultPtr<TH1D> Histo1D(RDF& rdf, const ROOT::RDF::TH1DModel& model, std::string_view vName)
{
    return rdf.template Histo1D<VAR1>(model, vName);
}

};
