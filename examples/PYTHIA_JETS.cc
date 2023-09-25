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

#include "cmath"
#include "iostream"
#include "fstream"
#include "random"
#include "vector"
#include "stdio.h"
#include "stdlib.h"

using namespace std;

namespace fastjet
{
  class ParticleInfo : public PseudoJet::UserInfoBase
  {
  public:
    ParticleInfo(const double charge, const int &pid) : _charge(charge), _pid(pid) {}
    double charge() const { return _charge; }
    int pid() const { return _pid; }

  protected:
    double _charge;
    int _pid;
  };
}

namespace Rivet
{
  class PYTHIA_JETS : public Analysis
  {
  public:
    RIVET_DEFAULT_ANALYSIS_CTOR(PYTHIA_JETS);

    void init()
    {
      const FinalState fs(Cuts::pt > 0.2 && Cuts::pt < 30 && Cuts::abseta < 1.0);
      declare(fs, "fs");

      // open output file
      mytxtfile.open("pythia_20pthat25.txt");
    }

    /// Perform the per-event analysis
    void analyze(const Event &event)
    {
      
      double xsecweight = handler().nominalCrossSection();
      HepMC::GenEvent *theEvent = (HepMC::GenEvent *)event.genEvent();
      int evid = theEvent->event_number();

      //===========================================
      // process particles
      Particles fsParticles = applyProjection<FinalState>(event, "fs").particles();
      PseudoJets parts;
      for (const Particle &p : fsParticles)
      {
        if (p.pid() == 12 || p.pid() == 14 || p.pid() == 16)
          continue;
        PseudoJet pseudojet(p.px(), p.py(), p.pz(), p.E());
        double charge = p.charge();
        int pid = p.pid();
        pseudojet.set_user_info(new fastjet::ParticleInfo(charge, pid));
        parts.push_back(pseudojet);
      }

      //===========================================
      // process jets
      fastjet::Selector selector = fastjet::SelectorPtMin(5.0) * fastjet::SelectorEtaRange(-0.6, 0.6);
      fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, 0.4);
      fastjet::contrib::SoftDrop sd(0, 0.1);
      fastjet::ClusterSequence cs(parts, jet_def);
      vector<fastjet::PseudoJet> jets = sorted_by_pt(selector(cs.inclusive_jets()));
      vector<fastjet::PseudoJet> sdjets = sd(jets);

      for (unsigned int i = 0; i < jets.size(); i++)
      {
        mytxtfile << evid << ", " << xsecweight << ", " << jets[i].perp() << ", " << jets[i].m() << ", " << jets[i].constituents().size() << ", " << sdjets[i].m() << ", " << sdjets[i].structure_of<fastjet::contrib::SoftDrop>().delta_R() << ", " << sdjets[i].structure_of<fastjet::contrib::SoftDrop>().symmetry() << ", " << sdjets[i].constituents().size() << ", " << i << ", " << theEvent->event_scale() << "\n";
      }
    }

    /// Normalise histograms etc., after the run
    void finalize()
    {
      double norm = crossSection() / sumOfWeights();
      std::cout << norm << std::endl;
      mytxtfile.close();
    }

  private:
    std::ofstream mytxtfile;
  };

  RIVET_DECLARE_PLUGIN(PYTHIA_JETS);
}