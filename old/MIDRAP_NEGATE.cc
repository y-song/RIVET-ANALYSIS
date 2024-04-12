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
  class MIDRAP_NEGATE : public Analysis
  {
  public:
    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(MIDRAP_NEGATE);

    void init()
    {
      // Initialise and register projections

      // final-state particles in mid rap
      const FinalState fs(Cuts::pt > 0.2 && Cuts::pt < 30 && Cuts::abseta < 1.0);
      declare(fs, "fs");
      // final-state particles in forward rap
      const FinalState fsf(Cuts::pt > 0.2 && Cuts::pt < 30 && Cuts::eta < 4.0 && Cuts::eta > 2.5);
      declare(fsf, "fsf");

      mytxtfile.open("high_midrap_mult_no_high_pt_neutral.txt");
    }

    /// Perform the per-event analysis
    void analyze(const Event &event)
    {

      // Get the final state particles and start clustering jets
      Particles fsParticles = applyProjection<FinalState>(event, "fs").particles();
      Particles fsfParticles = applyProjection<FinalState>(event, "fsf").particles();

      // jet selectors
      fastjet::Selector fwd_selector = fastjet::SelectorPtMin(5.0)*fastjet::SelectorRapRange(2.9,3.6);
      //fastjet::Selector mid_selector = fastjet::SelectorPtMin(4.0)*fastjet::SelectorRapRange(-1.0,1.0);
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
      if (fwd_jets.size() > 0){
	      cout << fwd_jets.size() << " " << mid_parts_ch.size() << " " << mid_parts_neutral_sorted.size();
      if (mid_parts_neutral_sorted.size() > 0) {
	      cout << " " << mid_parts_neutral_sorted[0].perp();
      }
      cout << "\n";
      }
      /*     //cout << mid_parts_neutral_sorted.size() << " " << mid_parts_ch.size() << endl;
 
      if (fwd_jets.size() > 0 && mid_parts_ch.size() > 20 && mid_parts_neutral_sorted.size() > 0 && mid_parts_neutral_sorted[0].perp() < 4.0) {
	      // match: forward jet, high midrap multiplicity, no high pT midrap neutral parts
	      mytxtfile << "1, " << fwd_jets[0].perp() << ", " << fwd_jets[0].constituents().size() << ", " << mid_parts_ch.size() << ", " << mid_parts_neutral_sorted.size() << ", " << mid_parts_neutral_sorted[0].perp() << "\n";
       cout << "1-1 " << mid_parts_neutral_sorted.size() << " " << mid_parts_ch.size() << endl;
 
}
      else if (fwd_jets.size() > 0 && mid_parts_ch.size() > 20 && mid_parts_neutral_sorted.size() == 0) {
	      mytxtfile << "1, " << fwd_jets[0].perp() << ", " << fwd_jets[0].constituents().size() << ", " << mid_parts_ch.size() << ", " << mid_parts_neutral_sorted.size() << ", -999\n";
      cout << "1-2 " << mid_parts_neutral_sorted.size() << " " << mid_parts_ch.size() << endl;
 
}

      else if (fwd_jets.size() > 0 && mid_parts_neutral_sorted.size() > 0) {
	      // miss: forward jet, but did not pass midrap selection
	      mytxtfile << "2, " << fwd_jets[0].perp() << ", " << fwd_jets[0].constituents().size() << ", " << mid_parts_ch.size() << ", " << mid_parts_neutral_sorted.size() << ", " << mid_parts_neutral_sorted[0].perp() << "\n";
       cout << "2-1 " << mid_parts_neutral_sorted.size() << " " << mid_parts_ch.size() << endl;
 
}
      else if (fwd_jets.size() > 0 && mid_parts_neutral_sorted.size() == 0) {
	      mytxtfile << "2, " << fwd_jets[0].perp() << ", " << fwd_jets[0].constituents().size() << ", " << mid_parts_ch.size() << ", " << mid_parts_neutral_sorted.size() << ", -999\n";
       cout << "2-2 " << mid_parts_neutral_sorted.size() << " " << mid_parts_ch.size() << endl;
 
}

      else if (mid_parts_ch.size() > 20 && mid_parts_neutral_sorted.size() > 0 && mid_parts_neutral_sorted[0].perp() < 4.0) {
	      // fake: high midrap mult, no high pT midrap neutral parts, but no forward jet
	      mytxtfile << "3, -999, -999, " << mid_parts_ch.size() << ", " << mid_parts_neutral_sorted.size() << ", " << mid_parts_neutral_sorted[0].perp() << "\n";
       cout << "3-1 " << mid_parts_neutral_sorted.size() << " " << mid_parts_ch.size() << endl;
 
}
      else if (mid_parts_ch.size() > 20 && mid_parts_neutral_sorted.size() == 0) {
	      mytxtfile << "3, -999, -999, " << mid_parts_ch.size() << ", " << mid_parts_neutral_sorted.size() << ", -999" << "\n";
       cout << "3-2 " << mid_parts_neutral_sorted.size() << " " << mid_parts_ch.size() << endl;
 
}

      else if (fwd_jets.size() == 0 && mid_parts_ch.size() <= 20) {
	      // neither
       //cout << mid_parts_neutral_sorted.size() << " " << mid_parts_ch.size() << endl;
 
}
      else if (fwd_jets.size() == 0 && mid_parts_neutral_sorted.size() > 0 && mid_parts_neutral_sorted[0].perp() > 4.0) {
	      // neither
       //cout << mid_parts_neutral_sorted.size() << " " << mid_parts_ch.size() << endl;
 
}
      else {
	      cout << "error: " << endl;
	      cout << fwd_jets[0].perp() << ", " << fwd_jets[0].constituents().size() << ", " << mid_parts_ch.size() << ", " << mid_parts_neutral_sorted.size() << "\n";
       cout << mid_parts_neutral_sorted.size() << " " << mid_parts_ch.size() << endl;
 
}*/
    }
    
    void finalize()
    {
      mytxtfile.close();
    }

  private:
    std::ofstream mytxtfile;
  };

  RIVET_DECLARE_PLUGIN(MIDRAP_NEGATE);
}
