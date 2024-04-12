// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/AnalysisHandler.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "fastjet/contrib/SoftDrop.hh"
#include "HepMC/GenEvent.h"
#include "HepMC/GenParticle.h"

#include "iostream"
#include "fstream"
#include "random"
#include "vector"

using namespace std;

namespace fastjet
{
  class MyUserInfo : public PseudoJet::UserInfoBase
  {
  public:
    MyUserInfo(vector<double> q, const int &matching, const int &parton) : _q(q), _matching(matching), _parton(parton) {}
    vector<double> q() const { return _q; }
    int matching() const { return _matching; }
    int parton() const { return _parton; }

  protected:
    vector<double> _q;
    int _matching;
    int _parton;
  };

  class ParticleInfo : public PseudoJet::UserInfoBase
  {
  public:
    ParticleInfo(const double charge) : _charge(charge) {}
    double charge() const { return _charge; }

  protected:
    double _charge;
  };
}

namespace Rivet
{
  class FORWARD_JETS : public Analysis
  {
  public:
    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(FORWARD_JETS);

    void init()
    {
      // Initialise and register projections

      // final-state particles in mid rap
      const FinalState fs(Cuts::pt > 0.2 && Cuts::pt < 30 && Cuts::abseta < 1.0);
      declare(fs, "fs");
      // final-state particles in forward rap
      const FinalState fsf(Cuts::pt > 0.2 && Cuts::pt < 30 && Cuts::eta < 4.0 && Cuts::eta > 2.5);
      declare(fsf, "fsf");

      mytxtfile.open("fwd_jets.txt");
    }

    /// Perform the per-event analysis
    void analyze(const Event &event)
    {

      // Get the final state particles and start clustering jets
      Particles fsParticles = applyProjection<FinalState>(event, "fs").particles();
      Particles fsfParticles = applyProjection<FinalState>(event, "fsf").particles();

      // jet selectors
      fastjet::Selector fwd_selector = fastjet::SelectorPtMin(5.0)*fastjet::SelectorRapRange(2.9,3.6);
      fastjet::Selector mid_selector = fastjet::SelectorPtMin(5.0)*fastjet::SelectorRapRange(-0.6,0.6);
      fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, 0.4);

      // cluster foward jets
      PseudoJets fwd_parts;
      for (const Particle &p : fsfParticles)
      {
        PseudoJet pseudojet(p.px(), p.py(), p.pz(), p.E());
        fwd_parts.push_back(pseudojet);
      }
      fastjet::ClusterSequence fwd_cs(fwd_parts, jet_def);
      vector<fastjet::PseudoJet> fwd_jets = sorted_by_pt(fwd_selector(fwd_cs.inclusive_jets()));

      // cluster midrap jets
      PseudoJets mid_parts;
      for (const Particle &p : fsParticles)
      {
        PseudoJet pseudojet(p.px(), p.py(), p.pz(), p.E());
        mid_parts.push_back(pseudojet);
      }
      fastjet::ClusterSequence mid_cs(mid_parts, jet_def);
      vector<fastjet::PseudoJet> mid_jets = sorted_by_pt(mid_selector(mid_cs.inclusive_jets()));

      // midrap neutral particles
      PseudoJets mid_parts_neutral;
      PseudoJets mid_parts_ch;

      for (const Particle &p : fsParticles)
      {
	      PseudoJet pseudojet(p.px(), p.py(), p.pz(), p.E());
      	      if (!p.isCharged()) {
     		mid_parts_neutral.push_back(pseudojet);
	      }
	      else {
	        mid_parts_ch.push_back(pseudojet);
	      }
      }
      vector<fastjet::PseudoJet> mid_parts_neutral_sorted = sorted_by_pt(mid_parts_neutral);
      
      mytxtfile << fwd_jets.size() << ", " << mid_jets.size() << ", " << mid_parts_ch.size() << ", " << mid_parts_neutral_sorted.size() << ", ";
      if (mid_parts_neutral_sorted.size() > 0) {
	   	mytxtfile << mid_parts_neutral_sorted[0].perp();
      }
      else {
		mytxtfile << "0";
      }
      if (mid_jets.size() > 0) {
	      mytxtfile << ", " << mid_jets[0].perp();
      }
      else {
	      mytxtfile << ", 0";
      }
      mytxtfile << "\n";
      
    }
    
    void finalize()
    {
      mytxtfile.close();
    }

  private:
    std::ofstream mytxtfile;
  };

  RIVET_DECLARE_PLUGIN(FORWARD_JETS);
}
