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
  class PTHAT5_TEST : public Analysis
  {
  public:
    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(PTHAT5_TEST);

    void init()
    {
      // Initialise and register projections

      // final-state particles in mid rap
      const FinalState fs(Cuts::pt > 0.2 && Cuts::pt < 30 && Cuts::abseta < 5.0);
      declare(fs, "fs");

      mytxtfile.open("pthat5_jets_all_eta_all_pt.txt");
    }

    /// Perform the per-event analysis
    void analyze(const Event &event)
    {

      // Get the final state particles and start clustering jets
      Particles fsParticles = applyProjection<FinalState>(event, "fs").particles();

      // jet selectors
      //fastjet::Selector mid_selector = fastjet::SelectorRapRange(-4.6,4.6);
      fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, 0.4);

      // cluster midrap jets
      PseudoJets mid_parts;
      for (const Particle &p : fsParticles)
      {
        PseudoJet pseudojet(p.px(), p.py(), p.pz(), p.E());
        mid_parts.push_back(pseudojet);
      }
      fastjet::ClusterSequence mid_cs(mid_parts, jet_def);
      vector<fastjet::PseudoJet> mid_jets = sorted_by_pt(mid_cs.inclusive_jets());
  
      mytxtfile << mid_jets.size() << ", ";
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

  RIVET_DECLARE_PLUGIN(PTHAT5_TEST);
}
